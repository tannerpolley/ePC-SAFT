#include "route_builders.h"

#include "eos_phase_block.h"
#include "epcsaft_core_internal.h"
#include "ipopt_adapter.h"
#include "nlp_problem.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace epcsaft::native::equilibrium_nlp {

namespace {

constexpr double kGasConstant = 8.31446261815324;
constexpr double kContractPhaseDistance = 1.0e-8;

struct NeutralTwoPhaseEosInitialPoint {
    std::vector<std::vector<double>> phase_amounts;
    std::vector<double> volumes;
};

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

double positive_sum(const std::vector<double>& values, const std::string& label) {
    if (values.empty()) {
        throw ValueError(label + " requires at least one value.");
    }
    double total = 0.0;
    for (double value : values) {
        require_positive_finite(value, label + " value");
        total += value;
    }
    require_positive_finite(total, label + " total");
    return total;
}

std::vector<double> normalized_positive_values(const std::vector<double>& values, const std::string& label) {
    const double total = positive_sum(values, label);
    std::vector<double> normalized;
    normalized.reserve(values.size());
    for (double value : values) {
        normalized.push_back(value / total);
    }
    return normalized;
}

std::vector<double> deterministic_composition_shift(
    const std::vector<double>& composition,
    const std::vector<double>& charges,
    const std::string& route_label
) {
    std::vector<double> shifted = composition;
    if (!charges.empty()) {
        require_size(charges, composition.size(), route_label + " charge");
    }
    if (composition.size() <= 1) {
        if (!charges.empty()) {
            throw ValueError(route_label + " requires at least two species.");
        }
        return shifted;
    }

    std::vector<double> positions;
    positions.reserve(composition.size());
    const double denominator = static_cast<double>(composition.size() - 1);
    for (std::size_t index = 0; index < composition.size(); ++index) {
        positions.push_back(-1.0 + 2.0 * static_cast<double>(index) / denominator);
    }
    double weighted_position = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        weighted_position += composition[index] * positions[index];
    }
    std::vector<double> direction;
    direction.reserve(composition.size());
    double max_abs_direction = 0.0;
    for (double position : positions) {
        const double value = position - weighted_position;
        direction.push_back(value);
        max_abs_direction = std::max(max_abs_direction, std::abs(value));
    }

    if (!charges.empty()) {
        double composition_charge = 0.0;
        double charge_square_weight = 0.0;
        double charge_direction = 0.0;
        for (std::size_t index = 0; index < composition.size(); ++index) {
            composition_charge += composition[index] * charges[index];
            charge_square_weight += composition[index] * charges[index] * charges[index];
            charge_direction += composition[index] * charges[index] * direction[index];
        }
        if (charge_square_weight <= 0.0) {
            throw ValueError(route_label + " requires at least one charged species.");
        }
        if (std::abs(composition_charge) > 1.0e-10) {
            throw ValueError(route_label + " feed must be charge neutral.");
        }
        const double charge_projection = charge_direction / charge_square_weight;
        max_abs_direction = 0.0;
        for (std::size_t index = 0; index < direction.size(); ++index) {
            direction[index] -= charge_projection * charges[index];
            max_abs_direction = std::max(max_abs_direction, std::abs(direction[index]));
        }
        if (max_abs_direction <= 0.0) {
            throw ValueError(route_label + " could not construct a charge-neutral initial direction.");
        }
    }

    double shifted_sum = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        const double scaled_direction = max_abs_direction > 0.0 ? direction[index] / max_abs_direction : 0.0;
        shifted[index] = composition[index] * (1.0 + 0.2 * scaled_direction);
        shifted_sum += shifted[index];
    }
    require_positive_finite(shifted_sum, "shifted composition sum");
    for (double& value : shifted) {
        value /= shifted_sum;
    }
    return shifted;
}

NeutralTwoPhaseEosInitialPoint build_two_phase_eos_initial_point(
    const std::vector<double>& feed_amounts,
    const std::vector<double>& first_composition,
    double temperature,
    double target_pressure,
    const std::string& route_label
) {
    require_positive_finite(temperature, route_label + " temperature");
    require_positive_finite(target_pressure, route_label + " pressure");
    const double total_feed = positive_sum(feed_amounts, route_label + " feed amount");
    require_size(first_composition, feed_amounts.size(), route_label + " first phase composition");

    NeutralTwoPhaseEosInitialPoint out;
    out.phase_amounts.assign(2, std::vector<double>(feed_amounts.size(), 0.0));
    for (std::size_t index = 0; index < feed_amounts.size(); ++index) {
        out.phase_amounts[0][index] = 0.5 * total_feed * first_composition[index];
        out.phase_amounts[1][index] = feed_amounts[index] - out.phase_amounts[0][index];
        require_positive_finite(out.phase_amounts[0][index], route_label + " first phase amount");
        require_positive_finite(out.phase_amounts[1][index], route_label + " second phase amount");
    }

    const double density = std::max(target_pressure / (kGasConstant * temperature), 1.0e-12);
    out.volumes.reserve(2);
    for (const auto& amounts : out.phase_amounts) {
        const double phase_total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
        require_positive_finite(phase_total, route_label + " phase amount total");
        out.volumes.push_back(phase_total / density);
    }
    return out;
}

NeutralTwoPhaseEosInitialPoint build_neutral_two_phase_eos_initial_point(
    const std::vector<double>& feed_amounts,
    double temperature,
    double target_pressure,
    const std::string& route_label
) {
    const std::vector<double> composition = normalized_positive_values(feed_amounts, route_label + " feed amount");
    const std::vector<double> first_composition =
        deterministic_composition_shift(composition, {}, route_label + " composition");
    return build_two_phase_eos_initial_point(feed_amounts, first_composition, temperature, target_pressure, route_label);
}

NeutralTwoPhaseEosInitialPoint build_charge_neutral_two_phase_eos_initial_point(
    const std::vector<double>& feed_amounts,
    const std::vector<double>& charges,
    double temperature,
    double target_pressure,
    const std::string& route_label
) {
    const std::vector<double> composition = normalized_positive_values(feed_amounts, route_label + " feed amount");
    const std::vector<double> first_composition =
        deterministic_composition_shift(composition, charges, route_label + " composition");
    return build_two_phase_eos_initial_point(feed_amounts, first_composition, temperature, target_pressure, route_label);
}

double vector_infinity_norm(const std::vector<double>& values, std::size_t begin, std::size_t end) {
    double norm = 0.0;
    for (std::size_t index = begin; index < end; ++index) {
        norm = std::max(norm, std::abs(values[index]));
    }
    return norm;
}

double phase_distance_inf_norm(const std::vector<double>& first, const std::vector<double>& second) {
    require_size(second, first.size(), "Neutral EOS postsolve phase composition");
    double distance = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        distance = std::max(distance, std::abs(first[index] - second[index]));
    }
    return distance;
}

double phase_charge_inf_norm(
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& charges
) {
    if (charges.empty()) {
        return 0.0;
    }
    double norm = 0.0;
    for (const auto& phase : phase_amounts) {
        require_size(phase, charges.size(), "Electrolyte LLE phase charge");
        double charge = 0.0;
        for (std::size_t species = 0; species < charges.size(); ++species) {
            charge += phase[species] * charges[species];
        }
        norm = std::max(norm, std::abs(charge));
    }
    return norm;
}

double chemical_potential_inf_norm(
    const EosPhaseBlockResult& first,
    const EosPhaseBlockResult& second,
    std::size_t species_count
) {
    if (first.gradient.size() < species_count || second.gradient.size() < species_count) {
        throw ValueError("Neutral EOS postsolve phase gradient sizes did not match species count.");
    }
    double norm = 0.0;
    for (std::size_t species = 0; species < species_count; ++species) {
        norm = std::max(norm, std::abs(first.gradient[species] - second.gradient[species]));
    }
    return norm;
}

std::vector<double> reduced_ln_fugacity_values(
    const add_args& args,
    const EosPhaseBlockResult& block,
    std::size_t species_count
) {
    if (block.composition.size() < species_count) {
        throw ValueError("Neutral EOS postsolve phase composition size did not match species count.");
    }
    FugacityContributionResult fugacity = fugacity_coefficient_result_cpp(
        block.temperature,
        block.density,
        block.composition,
        args
    );
    if (fugacity.lnfugcoef.total.size() < species_count) {
        throw ValueError("Neutral EOS postsolve fugacity payload size did not match species count.");
    }
    std::vector<double> values;
    values.reserve(species_count);
    for (std::size_t species = 0; species < species_count; ++species) {
        require_positive_finite(block.composition[species], "Neutral EOS postsolve phase composition");
        values.push_back(std::log(block.composition[species]) + fugacity.lnfugcoef.total[species]);
    }
    return values;
}

double ln_fugacity_inf_norm(
    const add_args& args,
    const EosPhaseBlockResult& first,
    const EosPhaseBlockResult& second,
    std::size_t species_count
) {
    const std::vector<double> first_values = reduced_ln_fugacity_values(args, first, species_count);
    const std::vector<double> second_values = reduced_ln_fugacity_values(args, second, species_count);
    double norm = 0.0;
    for (std::size_t species = 0; species < species_count; ++species) {
        norm = std::max(norm, std::abs(first_values[species] - second_values[species]));
    }
    return norm;
}

std::vector<std::vector<double>> neutral_phase_amounts_from_route_variables(
    const std::vector<double>& variables,
    std::size_t species_count
) {
    const std::size_t local_variable_count = species_count + 1;
    require_size(variables, 2 * local_variable_count, "Neutral EOS route result variable");
    std::vector<std::vector<double>> phase_amounts(2, std::vector<double>(species_count, 0.0));
    for (std::size_t phase = 0; phase < 2; ++phase) {
        const std::size_t offset = phase * local_variable_count;
        for (std::size_t species = 0; species < species_count; ++species) {
            phase_amounts[phase][species] = variables[offset + species];
        }
    }
    return phase_amounts;
}

std::vector<double> neutral_phase_volumes_from_route_variables(
    const std::vector<double>& variables,
    std::size_t species_count
) {
    const std::size_t local_variable_count = species_count + 1;
    require_size(variables, 2 * local_variable_count, "Neutral EOS route result variable");
    return {variables[species_count], variables[local_variable_count + species_count]};
}

class NeutralTwoPhaseEosProblem final : public NlpProblem {
public:
    NeutralTwoPhaseEosProblem(
        add_args args,
        double temperature,
        double target_pressure,
        std::vector<std::vector<double>> phase_amounts,
        std::vector<double> volumes,
        std::vector<double> feed_amounts,
        std::vector<double> charges = {},
        std::string problem_name = "neutral_two_phase_eos",
        double minimum_phase_distance = 0.0
    )
        : args_(std::move(args)),
          temperature_(temperature),
          target_pressure_(target_pressure),
          initial_phase_amounts_(std::move(phase_amounts)),
          initial_volumes_(std::move(volumes)),
          feed_amounts_(std::move(feed_amounts)),
          charges_(std::move(charges)),
          problem_name_(std::move(problem_name)),
          minimum_phase_distance_(minimum_phase_distance) {
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
        if (!charges_.empty()) {
            require_size(charges_, static_cast<std::size_t>(species_count_), "Two-phase EOS route charge");
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
        if (minimum_phase_distance_ > 0.0) {
            const std::vector<double> first = phase_composition(
                initial_phase_amounts_[0],
                "Neutral EOS route first phase amount total"
            );
            const std::vector<double> second = phase_composition(
                initial_phase_amounts_[1],
                "Neutral EOS route second phase amount total"
            );
            double max_distance = 0.0;
            for (int species = 0; species < species_count_; ++species) {
                const double diff =
                    first[static_cast<std::size_t>(species)] - second[static_cast<std::size_t>(species)];
                if (std::abs(diff) > max_distance) {
                    max_distance = std::abs(diff);
                    separation_species_index_ = species;
                    separation_sign_ = diff >= 0.0 ? 1.0 : -1.0;
                }
            }
            if (max_distance <= 0.0) {
                throw ValueError("Neutral EOS route requires distinct initial phases for phase-separation gating.");
            }
        }
    }

    std::string name() const override {
        return problem_name_;
    }

    int variable_count() const override {
        return phase_count() * local_variable_count();
    }

    int constraint_count() const override {
        return species_count_ + phase_count() + (charges_.empty() ? 0 : phase_count())
            + separation_constraint_count();
    }

    int jacobian_nonzero_count() const override {
        return variable_count() * constraint_count();
    }

    NlpBounds bounds() const override {
        NlpBounds out;
        const double total_feed = std::accumulate(feed_amounts_.begin(), feed_amounts_.end(), 0.0);
        const double amount_upper = std::max(1.0, 10.0 * total_feed);
        const double volume_upper = std::max(1.0, 1.0e6 * total_feed);
        out.variable_lower.assign(static_cast<std::size_t>(variable_count()), 1.0e-14);
        out.variable_upper.reserve(static_cast<std::size_t>(variable_count()));
        for (int phase = 0; phase < phase_count(); ++phase) {
            for (int species = 0; species < species_count_; ++species) {
                out.variable_upper.push_back(amount_upper);
            }
            out.variable_upper.push_back(volume_upper);
        }
        out.constraint_lower.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        if (separation_constraint_count() > 0) {
            out.constraint_lower.back() = minimum_phase_distance_;
            out.constraint_upper.back() = 1.0e12;
        }
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
        std::vector<double> out = phase_system(variables).constraints;
        if (separation_constraint_count() > 0) {
            out.push_back(phase_separation(phase_amounts_from_variables(variables)));
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
        std::vector<double> out = phase_system(variables).constraint_jacobian_row_major;
        if (separation_constraint_count() > 0) {
            const std::vector<double> row = phase_separation_jacobian(phase_amounts_from_variables(variables));
            out.insert(out.end(), row.begin(), row.end());
        }
        return out;
    }

    NlpScaling scaling() const override {
        const double total_feed = std::accumulate(feed_amounts_.begin(), feed_amounts_.end(), 0.0);
        const double amount_scale = std::max(1.0, total_feed);
        NlpScaling out;
        out.objective = 1.0 / amount_scale;
        out.variables.assign(static_cast<std::size_t>(variable_count()), 1.0 / amount_scale);
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0 / amount_scale);
        if (separation_constraint_count() > 0) {
            out.constraints.back() = 1.0;
        }
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

    int separation_constraint_count() const {
        return minimum_phase_distance_ > 0.0 ? 1 : 0;
    }

    std::vector<double> phase_composition(
        const std::vector<double>& amounts,
        const std::string& total_label
    ) const {
        require_size(amounts, static_cast<std::size_t>(species_count_), "Neutral EOS route phase amount");
        const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
        require_positive_finite(total, total_label);
        std::vector<double> composition;
        composition.reserve(amounts.size());
        for (double amount : amounts) {
            composition.push_back(amount / total);
        }
        return composition;
    }

    double phase_separation(const std::vector<std::vector<double>>& phase_amounts) const {
        const std::vector<double> first = phase_composition(
            phase_amounts[0],
            "Neutral EOS route first phase amount total"
        );
        const std::vector<double> second = phase_composition(
            phase_amounts[1],
            "Neutral EOS route second phase amount total"
        );
        const std::size_t species = static_cast<std::size_t>(separation_species_index_);
        return separation_sign_ * (first[species] - second[species]);
    }

    std::vector<double> phase_separation_jacobian(
        const std::vector<std::vector<double>>& phase_amounts
    ) const {
        const std::vector<double> first = phase_composition(
            phase_amounts[0],
            "Neutral EOS route first phase amount total"
        );
        const std::vector<double> second = phase_composition(
            phase_amounts[1],
            "Neutral EOS route second phase amount total"
        );
        const double first_total = std::accumulate(phase_amounts[0].begin(), phase_amounts[0].end(), 0.0);
        const double second_total = std::accumulate(phase_amounts[1].begin(), phase_amounts[1].end(), 0.0);
        require_positive_finite(first_total, "Neutral EOS route first phase amount total");
        require_positive_finite(second_total, "Neutral EOS route second phase amount total");

        std::vector<double> row(static_cast<std::size_t>(variable_count()), 0.0);
        const std::size_t selected = static_cast<std::size_t>(separation_species_index_);
        for (int species = 0; species < species_count_; ++species) {
            const std::size_t index = static_cast<std::size_t>(species);
            const double first_indicator = index == selected ? 1.0 : 0.0;
            row[index] = separation_sign_ * (first_indicator - first[selected]) / first_total;

            const std::size_t second_offset = static_cast<std::size_t>(local_variable_count() + species);
            const double second_indicator = index == selected ? 1.0 : 0.0;
            row[second_offset] = -separation_sign_ * (second_indicator - second[selected]) / second_total;
        }
        return row;
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
            charges_
        );
    }

    add_args args_;
    double temperature_ = 0.0;
    double target_pressure_ = 0.0;
    std::vector<std::vector<double>> initial_phase_amounts_;
    std::vector<double> initial_volumes_;
    std::vector<double> feed_amounts_;
    std::vector<double> charges_;
    std::string problem_name_;
    double minimum_phase_distance_ = 0.0;
    int separation_species_index_ = 0;
    double separation_sign_ = 1.0;
    int species_count_ = 0;
};

NeutralTwoPhaseEosNlpContract make_nlp_contract(
    const NlpProblem& problem,
    int phase_count,
    int species_count
) {
    validate_nlp_problem_shape(problem);

    const std::vector<double> initial = problem.initial_point();
    const NlpBounds bounds = problem.bounds();
    const NlpJacobianStructure structure = problem.jacobian_structure();

    NeutralTwoPhaseEosNlpContract out;
    out.problem_name = problem.name();
    out.derivative_backend = "analytic_cppad";
    out.phase_count = phase_count;
    out.species_count = species_count;
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

NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts
) {
    NeutralTwoPhaseEosProblem problem(args, temperature, target_pressure, phase_amounts, volumes, feed_amounts);
    return make_nlp_contract(problem, problem.phase_count(), problem.species_count());
}

NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_tp_flash_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts
) {
    const NeutralTwoPhaseEosInitialPoint initial =
        build_neutral_two_phase_eos_initial_point(feed_amounts, temperature, target_pressure, "Neutral TP flash route");
    NeutralTwoPhaseEosProblem problem(
        args,
        temperature,
        target_pressure,
        initial.phase_amounts,
        initial.volumes,
        feed_amounts,
        {},
        "neutral_two_phase_eos",
        kContractPhaseDistance
    );
    return make_nlp_contract(problem, problem.phase_count(), problem.species_count());
}

NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_lle_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts
) {
    const NeutralTwoPhaseEosInitialPoint initial =
        build_neutral_two_phase_eos_initial_point(feed_amounts, temperature, target_pressure, "Neutral LLE route");
    NeutralTwoPhaseEosProblem problem(
        args,
        temperature,
        target_pressure,
        initial.phase_amounts,
        initial.volumes,
        feed_amounts,
        {},
        "neutral_two_phase_eos",
        kContractPhaseDistance
    );
    return make_nlp_contract(problem, problem.phase_count(), problem.species_count());
}

NeutralTwoPhaseEosNlpContract evaluate_electrolyte_lle_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts
) {
    const NeutralTwoPhaseEosInitialPoint initial = build_charge_neutral_two_phase_eos_initial_point(
        feed_amounts,
        args.z,
        temperature,
        target_pressure,
        "Electrolyte LLE route"
    );
    NeutralTwoPhaseEosProblem problem(
        args,
        temperature,
        target_pressure,
        initial.phase_amounts,
        initial.volumes,
        feed_amounts,
        args.z,
        "electrolyte_lle_eos",
        kContractPhaseDistance
    );
    return make_nlp_contract(problem, problem.phase_count(), problem.species_count());
}

IpoptSolveResult solve_neutral_two_phase_eos_ipopt(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    const IpoptSolveOptions& options,
    const std::vector<double>& charges,
    const std::string& problem_name,
    double minimum_phase_distance
) {
    NeutralTwoPhaseEosProblem problem(
        args,
        temperature,
        target_pressure,
        phase_amounts,
        volumes,
        feed_amounts,
        charges,
        problem_name,
        minimum_phase_distance
    );
    return solve_ipopt_nlp(problem, options);
}

NeutralTwoPhaseEosPostsolve evaluate_neutral_two_phase_eos_postsolve(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    double material_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    require_positive_finite(material_tolerance, "Neutral EOS postsolve material tolerance");
    require_positive_finite(pressure_tolerance, "Neutral EOS postsolve pressure tolerance");
    require_positive_finite(
        chemical_potential_tolerance,
        "Neutral EOS postsolve chemical-potential tolerance"
    );
    require_positive_finite(phase_distance_tolerance, "Neutral EOS postsolve phase-distance tolerance");

    const EosPhaseSystemResult system = evaluate_eos_phase_system(
        args,
        temperature,
        target_pressure,
        phase_amounts,
        volumes,
        feed_amounts,
        {}
    );
    if (system.phase_count != 2) {
        throw ValueError("Neutral EOS postsolve currently requires exactly two phases.");
    }

    NeutralTwoPhaseEosPostsolve out;
    out.derivative_backend = system.derivative_backend;
    out.phase_count = system.phase_count;
    out.species_count = system.species_count;
    out.objective = system.objective;
    out.constraints = system.constraints;
    out.phase_volumes = volumes;
    out.phase_amount_totals.reserve(system.phase_blocks.size());
    out.phase_compositions.reserve(system.phase_blocks.size());
    for (const EosPhaseBlockResult& block : system.phase_blocks) {
        out.phase_amount_totals.push_back(block.total_amount);
        out.phase_compositions.push_back(block.composition);
    }

    const std::size_t species_count = static_cast<std::size_t>(system.species_count);
    out.material_balance_norm = vector_infinity_norm(system.constraints, 0, species_count);
    out.pressure_consistency_norm = vector_infinity_norm(
        system.constraints,
        species_count,
        system.constraints.size()
    );
    out.chemical_potential_consistency_norm = chemical_potential_inf_norm(
        system.phase_blocks[0],
        system.phase_blocks[1],
        species_count
    );
    out.ln_fugacity_consistency_norm = ln_fugacity_inf_norm(
        args,
        system.phase_blocks[0],
        system.phase_blocks[1],
        species_count
    );
    out.phase_distance = phase_distance_inf_norm(out.phase_compositions[0], out.phase_compositions[1]);

    out.accepted = out.material_balance_norm <= material_tolerance
        && out.pressure_consistency_norm <= pressure_tolerance
        && out.chemical_potential_consistency_norm <= chemical_potential_tolerance
        && out.ln_fugacity_consistency_norm <= chemical_potential_tolerance
        && out.phase_distance >= phase_distance_tolerance;
    if (out.accepted) {
        out.rejection_reason = "accepted";
    } else if (out.material_balance_norm > material_tolerance) {
        out.rejection_reason = "material_balance";
    } else if (out.pressure_consistency_norm > pressure_tolerance) {
        out.rejection_reason = "pressure_consistency";
    } else if (out.chemical_potential_consistency_norm > chemical_potential_tolerance) {
        out.rejection_reason = "chemical_potential_consistency";
    } else if (out.ln_fugacity_consistency_norm > chemical_potential_tolerance) {
        out.rejection_reason = "ln_fugacity_consistency";
    } else {
        out.rejection_reason = "phase_distance";
    }
    return out;
}

NeutralTwoPhaseEosRouteResult solve_neutral_two_phase_eos_route(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    const IpoptSolveOptions& options,
    double material_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance,
    double minimum_phase_distance,
    const std::vector<double>& charges,
    const std::string& problem_name,
    double charge_tolerance
) {
    if (!charges.empty()) {
        require_positive_finite(charge_tolerance, problem_name + " charge tolerance");
    }
    const IpoptAdapterInfo adapter = native_ipopt_adapter_info();
    NeutralTwoPhaseEosRouteResult out;
    out.compiled = adapter.compiled;
    out.adapter_available = adapter.adapter_available;
    out.adapter_kind = adapter.adapter_kind;
    out.hessian_strategy = adapter.hessian_strategy;
    out.exact_gradient_required = adapter.exact_gradient_required;
    out.exact_jacobian_required = adapter.exact_jacobian_required;
    out.problem_name = problem_name;
    if (!adapter.compiled) {
        out.status = "ipopt_dependency_required";
        return out;
    }

    const IpoptSolveResult solve = solve_neutral_two_phase_eos_ipopt(
        args,
        temperature,
        target_pressure,
        phase_amounts,
        volumes,
        feed_amounts,
        options,
        charges,
        problem_name,
        minimum_phase_distance
    );
    out.ran = solve.solver_ran;
    out.solver_accepted = solve.accepted;
    out.solver_status = solve.solver_status;
    out.application_status = solve.application_status;
    out.hessian_strategy = solve.hessian_strategy;
    out.objective = solve.objective;
    out.variables = solve.variables;
    out.constraints = solve.constraints;
    if (!solve.accepted) {
        out.status = "solver_rejected";
        return out;
    }

    const std::size_t species_count = feed_amounts.size();
    out.phase_amounts = neutral_phase_amounts_from_route_variables(solve.variables, species_count);
    out.phase_volumes = neutral_phase_volumes_from_route_variables(solve.variables, species_count);
    out.postsolve = evaluate_neutral_two_phase_eos_postsolve(
        args,
        temperature,
        target_pressure,
        out.phase_amounts,
        out.phase_volumes,
        feed_amounts,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
    if (!charges.empty()) {
        out.postsolve.charge_balance_norm = phase_charge_inf_norm(out.phase_amounts, charges);
        out.postsolve.accepted = out.postsolve.accepted && out.postsolve.charge_balance_norm <= charge_tolerance;
        if (out.postsolve.charge_balance_norm > charge_tolerance) {
            out.postsolve.rejection_reason = "charge_balance";
        }
    }
    out.accepted = out.postsolve.accepted;
    out.status = out.accepted ? "accepted" : "postsolve_rejected";
    return out;
}

NeutralTwoPhaseEosRouteResult solve_neutral_two_phase_eos_tp_flash_route(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts,
    const IpoptSolveOptions& options,
    double material_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    const NeutralTwoPhaseEosInitialPoint initial =
        build_neutral_two_phase_eos_initial_point(feed_amounts, temperature, target_pressure, "Neutral TP flash route");
    return solve_neutral_two_phase_eos_route(
        args,
        temperature,
        target_pressure,
        initial.phase_amounts,
        initial.volumes,
        feed_amounts,
        options,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
        phase_distance_tolerance
    );
}

NeutralTwoPhaseEosRouteResult solve_neutral_two_phase_eos_lle_route(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts,
    const IpoptSolveOptions& options,
    double material_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    const NeutralTwoPhaseEosInitialPoint initial =
        build_neutral_two_phase_eos_initial_point(feed_amounts, temperature, target_pressure, "Neutral LLE route");
    return solve_neutral_two_phase_eos_route(
        args,
        temperature,
        target_pressure,
        initial.phase_amounts,
        initial.volumes,
        feed_amounts,
        options,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
        phase_distance_tolerance
    );
}

NeutralTwoPhaseEosRouteResult solve_electrolyte_lle_eos_route(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts,
    const IpoptSolveOptions& options,
    double material_tolerance,
    double pressure_tolerance,
    double charge_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    const NeutralTwoPhaseEosInitialPoint initial = build_charge_neutral_two_phase_eos_initial_point(
        feed_amounts,
        args.z,
        temperature,
        target_pressure,
        "Electrolyte LLE route"
    );
    return solve_neutral_two_phase_eos_route(
        args,
        temperature,
        target_pressure,
        initial.phase_amounts,
        initial.volumes,
        feed_amounts,
        options,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance,
        phase_distance_tolerance,
        args.z,
        "electrolyte_lle_eos",
        charge_tolerance
    );
}

}  // namespace epcsaft::native::equilibrium_nlp
