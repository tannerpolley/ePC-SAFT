#include "route_builders.h"

#include "eos_phase_block.h"
#include "epcsaft_electrolyte.h"
#include "ipopt_adapter.h"
#include "nlp_problem.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

namespace epcsaft::native::equilibrium_nlp {

namespace {

void require_size(const std::vector<double>& values, std::size_t expected, const std::string& label) {
    if (values.size() == expected) {
        return;
    }
    throw ValueError(label + " size does not match the neutral two-phase EOS NLP contract.");
}

void require_positive_finite(double value, const std::string& label) {
    if (std::isfinite(value) && value > 0.0) {
        return;
    }
    throw ValueError(label + " must be positive and finite.");
}

class NeutralTwoPhaseEosProblem final : public NlpProblem {
public:
    NeutralTwoPhaseEosProblem(
        add_args args,
        double temperature,
        double target_pressure,
        std::vector<std::vector<double>> phase_amounts,
        std::vector<double> volumes,
        std::vector<double> feed_amounts
    )
        : args_(std::move(args)),
          temperature_(temperature),
          target_pressure_(target_pressure),
          initial_phase_amounts_(std::move(phase_amounts)),
          initial_volumes_(std::move(volumes)),
          feed_amounts_(std::move(feed_amounts)) {
        if (initial_phase_amounts_.size() != 2) {
            throw ValueError("Neutral EOS route builder currently requires exactly two phases.");
        }
        if (!std::isfinite(temperature_) || temperature_ <= 0.0 || !std::isfinite(target_pressure_)) {
            throw ValueError("Neutral EOS route builder received invalid T/P specifications.");
        }
        species_count_ = static_cast<int>(feed_amounts_.size());
        if (species_count_ <= 0) {
            throw ValueError("Neutral EOS route builder requires at least one feed species.");
        }
        require_size(initial_volumes_, initial_phase_amounts_.size(), "Neutral EOS route volume");
        for (const auto& amounts : initial_phase_amounts_) {
            require_size(amounts, static_cast<std::size_t>(species_count_), "Neutral EOS route phase amount");
            for (double amount : amounts) {
                require_positive_finite(amount, "Neutral EOS route phase amount");
            }
        }
        for (double volume : initial_volumes_) {
            require_positive_finite(volume, "Neutral EOS route phase volume");
        }
        for (double amount : feed_amounts_) {
            if (!std::isfinite(amount) || amount < 0.0) {
                throw ValueError("Neutral EOS route feed amounts must be finite and non-negative.");
            }
        }
    }

    std::string name() const override {
        return "neutral_two_phase_eos";
    }

    int variable_count() const override {
        return phase_count() * local_variable_count();
    }

    int constraint_count() const override {
        return species_count_ + phase_count();
    }

    int jacobian_nonzero_count() const override {
        return variable_count() * constraint_count();
    }

    NlpBounds bounds() const override {
        NlpBounds out;
        out.variable_lower.assign(static_cast<std::size_t>(variable_count()), 1.0e-14);
        out.variable_upper.assign(static_cast<std::size_t>(variable_count()), std::numeric_limits<double>::infinity());
        out.constraint_lower.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        return out;
    }

    std::vector<double> initial_point() const override {
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(variable_count()));
        for (std::size_t phase = 0; phase < initial_phase_amounts_.size(); ++phase) {
            out.insert(out.end(), initial_phase_amounts_[phase].begin(), initial_phase_amounts_[phase].end());
            out.push_back(initial_volumes_[phase]);
        }
        return out;
    }

    double objective(const std::vector<double>& variables) const override {
        return phase_system(variables).objective;
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        return phase_system(variables).gradient;
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        return phase_system(variables).constraints;
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
        return phase_system(variables).constraint_jacobian_row_major;
    }

    NlpScaling scaling() const override {
        const double total_feed = std::accumulate(feed_amounts_.begin(), feed_amounts_.end(), 0.0);
        const double amount_scale = std::max(1.0, total_feed);
        NlpScaling out;
        out.objective = 1.0 / amount_scale;
        out.variables.assign(static_cast<std::size_t>(variable_count()), 1.0 / amount_scale);
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0 / amount_scale);
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

    std::vector<std::vector<double>> phase_amounts_from_variables(const std::vector<double>& variables) const {
        require_size(variables, static_cast<std::size_t>(variable_count()), "Neutral EOS route variable");
        std::vector<std::vector<double>> phase_amounts(
            static_cast<std::size_t>(phase_count()),
            std::vector<double>(static_cast<std::size_t>(species_count_), 0.0)
        );
        for (int phase = 0; phase < phase_count(); ++phase) {
            const std::size_t offset = static_cast<std::size_t>(phase * local_variable_count());
            for (int species = 0; species < species_count_; ++species) {
                phase_amounts[static_cast<std::size_t>(phase)][static_cast<std::size_t>(species)] =
                    variables[offset + static_cast<std::size_t>(species)];
            }
        }
        return phase_amounts;
    }

    std::vector<double> volumes_from_variables(const std::vector<double>& variables) const {
        require_size(variables, static_cast<std::size_t>(variable_count()), "Neutral EOS route variable");
        std::vector<double> volumes(static_cast<std::size_t>(phase_count()), 0.0);
        for (int phase = 0; phase < phase_count(); ++phase) {
            const std::size_t volume_index = static_cast<std::size_t>(
                phase * local_variable_count() + species_count_
            );
            volumes[static_cast<std::size_t>(phase)] = variables[volume_index];
        }
        return volumes;
    }

    EosPhaseSystemResult phase_system(const std::vector<double>& variables) const {
        return evaluate_eos_phase_system(
            args_,
            temperature_,
            target_pressure_,
            phase_amounts_from_variables(variables),
            volumes_from_variables(variables),
            feed_amounts_,
            {}
        );
    }

    add_args args_;
    double temperature_ = 0.0;
    double target_pressure_ = 0.0;
    std::vector<std::vector<double>> initial_phase_amounts_;
    std::vector<double> initial_volumes_;
    std::vector<double> feed_amounts_;
    int species_count_ = 0;
};

}  // namespace

NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts
) {
    NeutralTwoPhaseEosProblem problem(args, temperature, target_pressure, phase_amounts, volumes, feed_amounts);
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

IpoptSolveResult solve_neutral_two_phase_eos_ipopt(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    const IpoptSolveOptions& options
) {
    NeutralTwoPhaseEosProblem problem(args, temperature, target_pressure, phase_amounts, volumes, feed_amounts);
    return solve_ipopt_nlp(problem, options);
}

}  // namespace epcsaft::native::equilibrium_nlp
