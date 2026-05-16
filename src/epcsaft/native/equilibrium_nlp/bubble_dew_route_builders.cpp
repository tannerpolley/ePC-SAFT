#include "route_builders.h"

#include "eos_phase_block.h"
#include "epcsaft_electrolyte.h"
#include "nlp_problem.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace epcsaft::native::equilibrium_nlp {
namespace {

constexpr double kGasConstant = 8.31446261815324;

void require_size(const std::vector<double>& values, std::size_t expected, const std::string& label) {
    if (values.size() == expected) {
        return;
    }
    throw ValueError(label + " size does not match the fixed-temperature pressure route.");
}

void require_positive_finite(double value, const std::string& label) {
    if (std::isfinite(value) && value > 0.0) {
        return;
    }
    throw ValueError(label + " must be positive and finite.");
}

std::vector<double> normalized_positive_values(const std::vector<double>& values, const std::string& label) {
    if (values.empty()) {
        throw ValueError(label + " requires at least one value.");
    }
    double total = 0.0;
    for (double value : values) {
        require_positive_finite(value, label + " value");
        total += value;
    }
    require_positive_finite(total, label + " total");
    std::vector<double> normalized;
    normalized.reserve(values.size());
    for (double value : values) {
        normalized.push_back(value / total);
    }
    return normalized;
}

std::vector<double> shifted_composition(const std::vector<double>& composition) {
    if (composition.size() <= 1) {
        return composition;
    }
    const double triangular_sum = 0.5 * static_cast<double>(composition.size() * (composition.size() + 1));
    std::vector<double> shifted;
    shifted.reserve(composition.size());
    for (std::size_t index = 0; index < composition.size(); ++index) {
        const double triangular = static_cast<double>(index + 1) / triangular_sum;
        shifted.push_back(0.8 * composition[index] + 0.2 * triangular);
    }
    return normalized_positive_values(shifted, "shifted composition");
}

class NeutralFixedTemperaturePressureProblem final : public NlpProblem {
public:
    NeutralFixedTemperaturePressureProblem(
        add_args args,
        double temperature,
        std::vector<double> fixed_composition,
        int fixed_phase_index,
        std::string problem_name
    )
        : args_(std::move(args)),
          temperature_(temperature),
          fixed_composition_(normalized_positive_values(fixed_composition, problem_name + " composition")),
          fixed_phase_index_(fixed_phase_index),
          problem_name_(std::move(problem_name)) {
        require_positive_finite(temperature_, problem_name_ + " temperature");
        if (fixed_phase_index_ < 0 || fixed_phase_index_ >= phase_count()) {
            throw ValueError(problem_name_ + " fixed phase index is out of range.");
        }
        species_count_ = static_cast<int>(fixed_composition_.size());
    }

    std::string name() const override {
        return problem_name_;
    }

    int variable_count() const override {
        return phase_count() * local_variable_count() + 1;
    }

    int constraint_count() const override {
        return composition_constraint_count() + 2 * phase_count();
    }

    int jacobian_nonzero_count() const override {
        return variable_count() * constraint_count();
    }

    NlpBounds bounds() const override {
        NlpBounds out;
        out.variable_lower.reserve(static_cast<std::size_t>(variable_count()));
        out.variable_upper.reserve(static_cast<std::size_t>(variable_count()));
        for (int phase = 0; phase < phase_count(); ++phase) {
            for (int species = 0; species < species_count_; ++species) {
                out.variable_lower.push_back(1.0e-14);
                out.variable_upper.push_back(10.0);
            }
            out.variable_lower.push_back(1.0e-14);
            out.variable_upper.push_back(1.0e6);
        }
        out.variable_lower.push_back(1.0);
        out.variable_upper.push_back(1.0e9);
        out.constraint_lower.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        return out;
    }

    std::vector<double> initial_point() const override {
        const std::vector<double> shifted = shifted_composition(fixed_composition_);
        const double density = std::max(initial_pressure_ / (kGasConstant * temperature_), 1.0e-12);
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(variable_count()));
        for (int phase = 0; phase < phase_count(); ++phase) {
            const std::vector<double>& composition = phase == fixed_phase_index_ ? fixed_composition_ : shifted;
            out.insert(out.end(), composition.begin(), composition.end());
            out.push_back(1.0 / density);
        }
        out.push_back(initial_pressure_);
        return out;
    }

    double objective(const std::vector<double>& variables) const override {
        double value = 0.0;
        for (const EosPhaseBlockResult& block : phase_blocks(variables)) {
            value += block.objective;
        }
        return value;
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(variable_count()));
        double pressure_derivative = 0.0;
        for (const EosPhaseBlockResult& block : phase_blocks(variables)) {
            out.insert(out.end(), block.gradient.begin(), block.gradient.end());
            pressure_derivative += block.volume / block.gas_constant_temperature;
        }
        out.push_back(pressure_derivative);
        return out;
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        const auto amounts = phase_amounts_from_variables(variables);
        std::vector<double> out(static_cast<std::size_t>(constraint_count()), 0.0);
        int row = 0;
        const auto& fixed_amounts = amounts[static_cast<std::size_t>(fixed_phase_index_)];
        const double fixed_total = std::accumulate(fixed_amounts.begin(), fixed_amounts.end(), 0.0);
        for (int species = 0; species < composition_constraint_count(); ++species) {
            out[static_cast<std::size_t>(row++)] =
                fixed_amounts[static_cast<std::size_t>(species)]
                - fixed_composition_[static_cast<std::size_t>(species)] * fixed_total;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const auto& phase_amounts = amounts[static_cast<std::size_t>(phase)];
            out[static_cast<std::size_t>(row++)] =
                std::accumulate(phase_amounts.begin(), phase_amounts.end(), 0.0) - 1.0;
        }
        for (const EosPhaseBlockResult& block : phase_blocks(variables)) {
            out[static_cast<std::size_t>(row++)] = block.pressure_consistency_residual;
        }
        return out;
    }

    NlpJacobianStructure jacobian_structure() const override {
        NlpJacobianStructure out;
        out.rows.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        out.cols.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < constraint_count(); ++row) {
            for (int col = 0; col < variable_count(); ++col) {
                out.rows.push_back(row);
                out.cols.push_back(col);
            }
        }
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        const auto blocks = phase_blocks(variables);
        std::vector<double> out(
            static_cast<std::size_t>(constraint_count() * variable_count()),
            0.0
        );
        int row = 0;
        for (int species = 0; species < composition_constraint_count(); ++species) {
            for (int col = 0; col < species_count_; ++col) {
                out[static_cast<std::size_t>(row * variable_count() + fixed_col(col))] =
                    (col == species ? 1.0 : 0.0) - fixed_composition_[static_cast<std::size_t>(species)];
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const int offset = phase * local_variable_count();
            for (int species = 0; species < species_count_; ++species) {
                out[static_cast<std::size_t>(row * variable_count() + offset + species)] = 1.0;
            }
            ++row;
        }
        for (int phase = 0; phase < phase_count(); ++phase) {
            const EosPhaseBlockResult& block = blocks[static_cast<std::size_t>(phase)];
            if (block.pressure_jacobian_row_major.size() != static_cast<std::size_t>(local_variable_count())) {
                throw ValueError(problem_name_ + " pressure Jacobian size did not match variables.");
            }
            const int offset = phase * local_variable_count();
            for (int col = 0; col < local_variable_count(); ++col) {
                out[static_cast<std::size_t>(row * variable_count() + offset + col)] =
                    block.pressure_jacobian_row_major[static_cast<std::size_t>(col)];
            }
            out[static_cast<std::size_t>(row * variable_count() + variable_count() - 1)] = -1.0;
            ++row;
        }
        return out;
    }

    NlpScaling scaling() const override {
        NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(static_cast<std::size_t>(variable_count()), 1.0);
        out.variables.back() = 1.0e-5;
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0);
        return out;
    }

    int species_count() const {
        return species_count_;
    }

    int phase_count() const {
        return 2;
    }

private:
    int local_variable_count() const {
        return species_count_ + 1;
    }

    int composition_constraint_count() const {
        return std::max(0, species_count_ - 1);
    }

    int fixed_col(int species) const {
        return fixed_phase_index_ * local_variable_count() + species;
    }

    std::vector<std::vector<double>> phase_amounts_from_variables(const std::vector<double>& variables) const {
        require_size(variables, static_cast<std::size_t>(variable_count()), problem_name_ + " variable");
        std::vector<std::vector<double>> amounts(
            static_cast<std::size_t>(phase_count()),
            std::vector<double>(static_cast<std::size_t>(species_count_), 0.0)
        );
        for (int phase = 0; phase < phase_count(); ++phase) {
            const std::size_t offset = static_cast<std::size_t>(phase * local_variable_count());
            for (int species = 0; species < species_count_; ++species) {
                amounts[static_cast<std::size_t>(phase)][static_cast<std::size_t>(species)] =
                    variables[offset + static_cast<std::size_t>(species)];
            }
        }
        return amounts;
    }

    std::vector<double> volumes_from_variables(const std::vector<double>& variables) const {
        require_size(variables, static_cast<std::size_t>(variable_count()), problem_name_ + " variable");
        std::vector<double> volumes(static_cast<std::size_t>(phase_count()), 0.0);
        for (int phase = 0; phase < phase_count(); ++phase) {
            volumes[static_cast<std::size_t>(phase)] =
                variables[static_cast<std::size_t>(phase * local_variable_count() + species_count_)];
        }
        return volumes;
    }

    std::vector<EosPhaseBlockResult> phase_blocks(const std::vector<double>& variables) const {
        require_size(variables, static_cast<std::size_t>(variable_count()), problem_name_ + " variable");
        const double pressure = variables.back();
        require_positive_finite(pressure, problem_name_ + " pressure");
        const auto amounts = phase_amounts_from_variables(variables);
        const auto volumes = volumes_from_variables(variables);
        std::vector<EosPhaseBlockResult> blocks;
        blocks.reserve(static_cast<std::size_t>(phase_count()));
        for (int phase = 0; phase < phase_count(); ++phase) {
            blocks.push_back(
                evaluate_eos_phase_block(
                    args_,
                    temperature_,
                    pressure,
                    amounts[static_cast<std::size_t>(phase)],
                    volumes[static_cast<std::size_t>(phase)]
                )
            );
        }
        return blocks;
    }

    add_args args_;
    double temperature_ = 0.0;
    double initial_pressure_ = 1.0e5;
    std::vector<double> fixed_composition_;
    int fixed_phase_index_ = 0;
    std::string problem_name_;
    int species_count_ = 0;
};

NeutralTwoPhaseEosNlpContract make_contract(const NeutralFixedTemperaturePressureProblem& problem) {
    validate_nlp_problem_shape(problem);

    const std::vector<double> initial = problem.initial_point();
    const NlpBounds bounds = problem.bounds();
    const NlpJacobianStructure structure = problem.jacobian_structure();

    NeutralTwoPhaseEosNlpContract out;
    out.problem_name = problem.name();
    out.derivative_backend = "analytic_cppad";
    out.phase_count = problem.phase_count();
    out.species_count = problem.species_count();
    out.variable_count = problem.variable_count();
    out.constraint_count = problem.constraint_count();
    out.jacobian_nonzero_count = problem.jacobian_nonzero_count();
    out.initial_point = initial;
    out.variable_lower_bounds = bounds.variable_lower;
    out.variable_upper_bounds = bounds.variable_upper;
    out.constraint_lower_bounds = bounds.constraint_lower;
    out.constraint_upper_bounds = bounds.constraint_upper;
    out.objective_at_initial = problem.objective(initial);
    out.gradient_at_initial = problem.objective_gradient(initial);
    out.constraints_at_initial = problem.constraints(initial);
    out.jacobian_rows = structure.rows;
    out.jacobian_cols = structure.cols;
    out.jacobian_values_at_initial = problem.jacobian_values(initial);
    return out;
}

}  // namespace

NeutralTwoPhaseEosNlpContract evaluate_neutral_bubble_p_eos_nlp_contract(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition
) {
    NeutralFixedTemperaturePressureProblem problem(
        args,
        temperature,
        liquid_composition,
        0,
        "neutral_bubble_p_eos"
    );
    return make_contract(problem);
}

NeutralTwoPhaseEosNlpContract evaluate_neutral_dew_p_eos_nlp_contract(
    const add_args& args,
    double temperature,
    const std::vector<double>& vapor_composition
) {
    NeutralFixedTemperaturePressureProblem problem(
        args,
        temperature,
        vapor_composition,
        1,
        "neutral_dew_p_eos"
    );
    return make_contract(problem);
}

}  // namespace epcsaft::native::equilibrium_nlp
