#include "stability_route_builders.h"

#include "epcsaft_core_internal.h"
#include "nlp_problem.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace epcsaft::native::equilibrium_nlp {
namespace {

std::string solve_string(const IpoptSolveResult& solve, const std::string& key, const std::string& default_value) {
    const auto item = solve.diagnostics_string.find(key);
    return item == solve.diagnostics_string.end() ? default_value : item->second;
}

int solve_int(const IpoptSolveResult& solve, const std::string& key, int default_value = 0) {
    const auto item = solve.diagnostics_int.find(key);
    return item == solve.diagnostics_int.end() ? default_value : item->second;
}

double solve_double(const IpoptSolveResult& solve, const std::string& key, double default_value = 0.0) {
    const auto item = solve.diagnostics_double.find(key);
    return item == solve.diagnostics_double.end() ? default_value : item->second;
}

bool solve_bool(const IpoptSolveResult& solve, const std::string& key, bool default_value = false) {
    const auto item = solve.diagnostics_bool.find(key);
    return item == solve.diagnostics_bool.end() ? default_value : item->second;
}

void apply_stability_ipopt_metadata(StabilityRouteResult& out, const IpoptSolveResult& solve) {
    out.gradient_approximation = solve_string(solve, "gradient_approximation", "exact");
    out.jacobian_approximation = solve_string(solve, "jacobian_approximation", "exact");
    out.hessian_approximation = solve_string(solve, "hessian_approximation", out.hessian_approximation);
    out.hessian_backend = solve_string(solve, "hessian_backend", out.hessian_backend);
    out.scaling_method = solve_string(solve, "scaling_method", out.scaling_method);
    out.linear_solver_requested = solve_string(solve, "linear_solver_requested", out.linear_solver_requested);
    out.linear_solver_selected = solve_string(solve, "linear_solver_selected", out.linear_solver_selected);
    out.iteration_count = solve_int(solve, "iteration_count");
    out.iteration_history_limit = solve_int(solve, "iteration_history_limit");
    out.iteration_history_size = solve_int(solve, "iteration_history_size");
    out.variable_scaling_count = solve_int(solve, "variable_scaling_count");
    out.constraint_scaling_count = solve_int(solve, "constraint_scaling_count");
    out.eval_h_calls = solve_int(solve, "eval_h_calls");
    out.objective_scaling = solve_double(solve, "objective_scaling", out.objective_scaling);
    out.acceptable_tolerance = solve_double(solve, "acceptable_tolerance", out.acceptable_tolerance);
    out.constraint_violation_tolerance =
        solve_double(solve, "constraint_violation_tolerance", out.constraint_violation_tolerance);
    out.dual_infeasibility_tolerance =
        solve_double(solve, "dual_infeasibility_tolerance", out.dual_infeasibility_tolerance);
    out.complementarity_tolerance =
        solve_double(solve, "complementarity_tolerance", out.complementarity_tolerance);
    out.variable_scaling_min = solve_double(solve, "variable_scaling_min", out.variable_scaling_min);
    out.variable_scaling_max = solve_double(solve, "variable_scaling_max", out.variable_scaling_max);
    out.constraint_scaling_min = solve_double(solve, "constraint_scaling_min", out.constraint_scaling_min);
    out.constraint_scaling_max = solve_double(solve, "constraint_scaling_max", out.constraint_scaling_max);
    out.exact_hessian_available = solve_bool(solve, "exact_hessian_available");
    out.warm_start_requested = solve_bool(solve, "warm_start_requested");
    out.warm_start_used = solve_bool(solve, "warm_start_used");
    out.bound_lower_multipliers = solve.bound_lower_multipliers;
    out.bound_upper_multipliers = solve.bound_upper_multipliers;
    out.constraint_multipliers = solve.constraint_multipliers;
    out.iteration_history = solve.iteration_history;
}

void require_size(const std::vector<double>& values, std::size_t expected, const std::string& label) {
    if (values.size() == expected) {
        return;
    }
    throw ValueError(label + " size does not match the stability route.");
}

void require_positive_finite(double value, const std::string& label) {
    if (std::isfinite(value) && value > 0.0) {
        return;
    }
    throw ValueError(label + " must be positive and finite.");
}

std::string phase_label(int phase) {
    if (phase == 0) {
        return "liq";
    }
    if (phase == 1) {
        return "vap";
    }
    throw ValueError("stability route phase must be 0/liquid or 1/vapor.");
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

std::vector<double> shifted_composition(const std::vector<double>& composition, double shift_sign = 1.0) {
    if (composition.size() <= 1) {
        std::vector<double> out;
        out.reserve(composition.size());
        for (double value : composition) {
            out.push_back(value);
        }
        return out;
    }
    const double triangular_sum = 0.5 * static_cast<double>(composition.size() * (composition.size() + 1));
    std::vector<double> shifted;
    shifted.reserve(composition.size());
    double total = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        const double direction = static_cast<double>(index + 1) / triangular_sum - 1.0 / static_cast<double>(composition.size());
        const double value = composition[index] * (1.0 + 0.2 * shift_sign * direction);
        require_positive_finite(value, "stability route shifted composition");
        shifted.push_back(value);
        total += value;
    }
    require_positive_finite(total, "stability route shifted composition total");
    for (double& value : shifted) {
        value /= total;
    }
    return shifted;
}

std::vector<double> charge_neutral_shifted_composition(
    const std::vector<double>& composition,
    const std::vector<double>& charges,
    const std::string& label,
    double shift_sign = 1.0
) {
    require_size(charges, composition.size(), label + " charge");
    if (composition.size() <= 1) {
        throw ValueError(label + " requires at least two species.");
    }
    double composition_charge = 0.0;
    double charge_square_weight = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        if (!std::isfinite(charges[index])) {
            throw ValueError(label + " charge values must be finite.");
        }
        composition_charge += composition[index] * charges[index];
        charge_square_weight += composition[index] * charges[index] * charges[index];
    }
    if (charge_square_weight <= 0.0) {
        throw ValueError(label + " requires at least one charged species.");
    }
    if (std::abs(composition_charge) > 1.0e-10) {
        throw ValueError(label + " fixed composition must be charge neutral.");
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
    double charge_direction = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        const double value = positions[index] - weighted_position;
        direction.push_back(value);
        charge_direction += composition[index] * charges[index] * value;
    }
    const double charge_projection = charge_direction / charge_square_weight;
    double max_abs_direction = 0.0;
    for (std::size_t index = 0; index < direction.size(); ++index) {
        direction[index] -= charge_projection * charges[index];
        max_abs_direction = std::max(max_abs_direction, std::abs(direction[index]));
    }
    if (max_abs_direction <= 0.0) {
        throw ValueError(label + " could not construct a charge-neutral initial direction.");
    }

    std::vector<double> shifted;
    shifted.reserve(composition.size());
    double total = 0.0;
    for (std::size_t index = 0; index < composition.size(); ++index) {
        const double value =
            composition[index] * (1.0 + 0.2 * shift_sign * direction[index] / max_abs_direction);
        require_positive_finite(value, label + " shifted composition");
        shifted.push_back(value);
        total += value;
    }
    require_positive_finite(total, label + " shifted composition total");
    for (double& value : shifted) {
        value /= total;
    }
    return shifted;
}

struct NamedInitialComposition {
    std::string seed_name;
    std::vector<double> composition;
};

std::vector<NamedInitialComposition> stability_seed_candidates(
    const std::vector<double>& feed_composition,
    const std::vector<double>& charges,
    const std::vector<double>& trial_initial_composition
) {
    std::vector<NamedInitialComposition> out;
    if (!trial_initial_composition.empty()) {
        out.push_back({"provided_trial_initial_composition", trial_initial_composition});
    }
    if (charges.empty()) {
        out.push_back({"canonical_shifted_feed", shifted_composition(feed_composition, 1.0)});
        out.push_back({"mirrored_shifted_feed", shifted_composition(feed_composition, -1.0)});
        return out;
    }
    out.push_back({
        "canonical_charge_neutral_feed",
        charge_neutral_shifted_composition(
            feed_composition,
            charges,
            "stability trial initial composition",
            1.0
        )
    });
    out.push_back({
        "mirrored_charge_neutral_feed",
        charge_neutral_shifted_composition(
            feed_composition,
            charges,
            "stability trial initial composition",
            -1.0
        )
    });
    return out;
}

RouteSeedAttempt stability_seed_attempt_from_result(const StabilityRouteResult& result) {
    RouteSeedAttempt out;
    out.seed_name = result.seed_name;
    out.status = result.status;
    out.solver_status = result.solver_status;
    out.application_status = result.application_status;
    out.solver_accepted = result.solver_accepted;
    out.accepted = result.accepted;
    out.stable = result.stable;
    out.iteration_count = result.iteration_count;
    out.objective = result.objective;
    out.min_tpd = result.min_tpd;
    return out;
}

bool stability_attempt_better(const StabilityRouteResult& candidate, const StabilityRouteResult& current) {
    if (candidate.accepted != current.accepted) {
        return candidate.accepted;
    }
    if (candidate.accepted && current.accepted) {
        return candidate.min_tpd < current.min_tpd;
    }
    if (candidate.solver_accepted != current.solver_accepted) {
        return candidate.solver_accepted;
    }
    if (candidate.ran != current.ran) {
        return candidate.ran;
    }
    return false;
}

PhaseStateCompositionSensitivityResult phase_state_sensitivity(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& composition,
    int phase,
    const std::string& label
) {
    PhaseStateCompositionSensitivityResult result = phase_state_ln_fugacity_composition_sensitivity_cpp(
        temperature,
        pressure,
        composition,
        phase,
        args
    );
    if (!result.supported) {
        const std::string message = result.message.empty()
            ? "phase-state fugacity composition sensitivity was not available."
            : result.message;
        throw ValueError(label + " " + message);
    }
    require_size(result.ln_fugacity, composition.size(), label + " ln fugacity");
    require_size(
        result.jacobian_row_major,
        composition.size() * composition.size(),
        label + " ln fugacity composition Jacobian"
    );
    return result;
}

std::vector<double> reduced_potential(
    const std::vector<double>& composition,
    const std::vector<double>& ln_fugacity
) {
    require_size(ln_fugacity, composition.size(), "stability reduced potential");
    std::vector<double> out;
    out.reserve(composition.size());
    for (std::size_t index = 0; index < composition.size(); ++index) {
        require_positive_finite(composition[index], "stability composition");
        out.push_back(std::log(composition[index]) + ln_fugacity[index]);
    }
    return out;
}

class StabilityTpdProblem final : public NlpProblem {
public:
    StabilityTpdProblem(
        add_args args,
        double temperature,
        double pressure,
        std::vector<double> feed_composition,
        int parent_phase,
        int trial_phase,
        std::string problem_name,
        std::vector<double> charges = {},
        bool require_charge_constraint = false,
        std::vector<double> initial_composition = {}
    )
        : args_(std::move(args)),
          temperature_(temperature),
          pressure_(pressure),
          feed_composition_(normalized_positive_values(feed_composition, "stability feed")),
          parent_phase_(parent_phase),
          trial_phase_(trial_phase),
          parent_phase_label_(phase_label(parent_phase)),
          trial_phase_label_(phase_label(trial_phase)),
          problem_name_(std::move(problem_name)),
          charges_(std::move(charges)),
          initial_composition_(
              initial_composition.empty()
                  ? std::vector<double>{}
                  : normalized_positive_values(initial_composition, "stability trial initial composition")
          ) {
        require_positive_finite(temperature_, "stability temperature");
        require_positive_finite(pressure_, "stability pressure");
        species_count_ = static_cast<int>(feed_composition_.size());
        if (!initial_composition_.empty()) {
            require_size(
                initial_composition_,
                static_cast<std::size_t>(species_count_),
                "stability trial initial composition"
            );
        }
        if (charges_.empty()) {
            if (require_charge_constraint) {
                throw ValueError("electrolyte stability route requires charge data.");
            }
        } else {
            require_size(charges_, static_cast<std::size_t>(species_count_), "stability charge");
            bool has_charged_species = false;
            double feed_charge = 0.0;
            for (int index = 0; index < species_count_; ++index) {
                const double charge = charges_[static_cast<std::size_t>(index)];
                if (!std::isfinite(charge)) {
                    throw ValueError("stability charge values must be finite.");
                }
                has_charged_species = has_charged_species || std::abs(charge) > charge_epsilon_;
                feed_charge += feed_composition_[static_cast<std::size_t>(index)] * charge;
            }
            if (!has_charged_species) {
                throw ValueError("electrolyte stability route requires charged species.");
            }
            if (std::abs(feed_charge) > charge_balance_tolerance_) {
                throw ValueError("electrolyte stability feed must be charge neutral.");
            }
            if (!initial_composition_.empty()) {
                double initial_charge = 0.0;
                for (int index = 0; index < species_count_; ++index) {
                    initial_charge += initial_composition_[static_cast<std::size_t>(index)]
                        * charges_[static_cast<std::size_t>(index)];
                }
                if (std::abs(initial_charge) > charge_balance_tolerance_) {
                    throw ValueError("electrolyte stability trial initial composition must be charge neutral.");
                }
            }
        }
        parent_state_ = phase_state_sensitivity(
            args_,
            temperature_,
            pressure_,
            feed_composition_,
            parent_phase_,
            "parent stability state"
        );
        parent_reduced_potential_ = reduced_potential(feed_composition_, parent_state_.ln_fugacity);
    }

    std::string name() const override {
        return problem_name_;
    }

    int variable_count() const override {
        return species_count_;
    }

    int constraint_count() const override {
        return has_charge_constraint() ? 2 : 1;
    }

    int jacobian_nonzero_count() const override {
        return species_count_ * constraint_count();
    }

    NlpBounds bounds() const override {
        NlpBounds out;
        out.variable_lower.assign(static_cast<std::size_t>(species_count_), minimum_composition_);
        out.variable_upper.assign(static_cast<std::size_t>(species_count_), 1.0);
        out.constraint_lower.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        return out;
    }

    std::vector<double> initial_point() const override {
        if (!initial_composition_.empty()) {
            return initial_composition_;
        }
        if (has_charge_constraint()) {
            return feed_composition_;
        }
        return shifted_composition(feed_composition_);
    }

    double objective(const std::vector<double>& variables) const override {
        const PhaseStateCompositionSensitivityResult trial = trial_state(variables);
        double value = 0.0;
        for (int index = 0; index < species_count_; ++index) {
            const double xi = variables[static_cast<std::size_t>(index)];
            value += xi * (
                std::log(xi)
                + trial.ln_fugacity[static_cast<std::size_t>(index)]
                - parent_reduced_potential_[static_cast<std::size_t>(index)]
            );
        }
        return value;
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        const PhaseStateCompositionSensitivityResult trial = trial_state(variables);
        std::vector<double> gradient(static_cast<std::size_t>(species_count_), 0.0);
        for (int col = 0; col < species_count_; ++col) {
            const double xj = variables[static_cast<std::size_t>(col)];
            gradient[static_cast<std::size_t>(col)] =
                std::log(xj)
                + 1.0
                + trial.ln_fugacity[static_cast<std::size_t>(col)]
                - parent_reduced_potential_[static_cast<std::size_t>(col)];
            for (int row = 0; row < species_count_; ++row) {
                gradient[static_cast<std::size_t>(col)] +=
                    variables[static_cast<std::size_t>(row)]
                    * trial.jacobian_row_major[static_cast<std::size_t>(row * species_count_ + col)];
            }
        }
        return gradient;
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        require_trial_variables(variables);
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(constraint_count()));
        out.push_back(std::accumulate(variables.begin(), variables.end(), 0.0) - 1.0);
        if (has_charge_constraint()) {
            double charge_balance = 0.0;
            for (int index = 0; index < species_count_; ++index) {
                charge_balance +=
                    variables[static_cast<std::size_t>(index)] * charges_[static_cast<std::size_t>(index)];
            }
            out.push_back(charge_balance);
        }
        return out;
    }

    NlpJacobianStructure jacobian_structure() const override {
        NlpJacobianStructure out;
        out.rows.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        out.cols.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < constraint_count(); ++row) {
            for (int col = 0; col < species_count_; ++col) {
                out.rows.push_back(row);
                out.cols.push_back(col);
            }
        }
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        require_trial_variables(variables);
        std::vector<double> out(static_cast<std::size_t>(species_count_), 1.0);
        if (has_charge_constraint()) {
            out.insert(out.end(), charges_.begin(), charges_.end());
        }
        return out;
    }

    NlpScaling scaling() const override {
        NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(static_cast<std::size_t>(species_count_), 1.0);
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0);
        return out;
    }

    int species_count() const {
        return species_count_;
    }

    const std::string& parent_phase_label() const {
        return parent_phase_label_;
    }

    const std::string& trial_phase_label() const {
        return trial_phase_label_;
    }

    const std::vector<double>& feed_composition() const {
        return feed_composition_;
    }

    const std::vector<double>& parent_reduced_potential() const {
        return parent_reduced_potential_;
    }

private:
    bool has_charge_constraint() const {
        return !charges_.empty();
    }

    void require_trial_variables(const std::vector<double>& variables) const {
        require_size(variables, static_cast<std::size_t>(species_count_), "stability trial variable");
        for (double value : variables) {
            require_positive_finite(value, "stability trial composition");
        }
    }

    PhaseStateCompositionSensitivityResult trial_state(const std::vector<double>& variables) const {
        require_trial_variables(variables);
        return phase_state_sensitivity(
            args_,
            temperature_,
            pressure_,
            variables,
            trial_phase_,
            "trial stability state"
        );
    }

    add_args args_;
    double temperature_ = 0.0;
    double pressure_ = 0.0;
    double minimum_composition_ = 1.0e-14;
    double charge_epsilon_ = 1.0e-12;
    double charge_balance_tolerance_ = 1.0e-10;
    std::vector<double> feed_composition_;
    int parent_phase_ = 0;
    int trial_phase_ = 0;
    std::string parent_phase_label_;
    std::string trial_phase_label_;
    std::string problem_name_;
    std::vector<double> charges_;
    std::vector<double> initial_composition_;
    PhaseStateCompositionSensitivityResult parent_state_;
    std::vector<double> parent_reduced_potential_;
    int species_count_ = 0;
};

StabilityNlpContract make_contract(const StabilityTpdProblem& problem) {
    validate_nlp_problem_shape(problem);

    const std::vector<double> initial = problem.initial_point();
    const NlpBounds bounds = problem.bounds();
    const NlpJacobianStructure structure = problem.jacobian_structure();

    StabilityNlpContract out;
    out.problem_name = problem.name();
    out.derivative_backend = "cppad_implicit";
    out.species_count = problem.species_count();
    out.variable_count = problem.variable_count();
    out.constraint_count = problem.constraint_count();
    out.jacobian_nonzero_count = problem.jacobian_nonzero_count();
    out.parent_phase = problem.parent_phase_label();
    out.trial_phase = problem.trial_phase_label();
    out.feed_composition = problem.feed_composition();
    out.parent_reduced_potential = problem.parent_reduced_potential();
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

StabilityRouteResult solve_stability_tpd_route(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition,
    int parent_phase,
    int trial_phase,
    const std::string& problem_name,
    const std::string& seed_name,
    std::vector<double> charges,
    bool require_charge_constraint,
    const IpoptSolveOptions& options,
    double stability_tolerance,
    const std::vector<double>& trial_initial_composition
) {
    const IpoptAdapterInfo adapter = native_ipopt_adapter_info();
    StabilityRouteResult out;
    out.compiled = adapter.compiled;
    out.adapter_available = adapter.adapter_available;
    out.adapter_kind = adapter.adapter_kind;
    out.problem_name = problem_name;
    out.seed_name = seed_name;
    out.exact_gradient_required = adapter.exact_gradient_required;
    out.exact_jacobian_required = adapter.exact_jacobian_required;
    out.parent_phase = phase_label(parent_phase);
    out.trial_phase = phase_label(trial_phase);
    if (!adapter.compiled) {
        out.status = "ipopt_dependency_required";
        return out;
    }

    StabilityTpdProblem problem(
        args,
        temperature,
        pressure,
        feed_composition,
        parent_phase,
        trial_phase,
        problem_name,
        std::move(charges),
        require_charge_constraint,
        trial_initial_composition
    );
    out.parent_reduced_potential = problem.parent_reduced_potential();
    out.initial_composition = problem.initial_point();
    const IpoptSolveResult solve = solve_ipopt_nlp(problem, options);
    out.ran = solve.solver_ran;
    out.solver_accepted = solve.accepted;
    out.solver_status = solve.solver_status;
    out.application_status = solve.application_status;
    apply_stability_ipopt_metadata(out, solve);
    out.objective = solve.objective;
    out.min_tpd = solve.objective;
    out.variables = solve.variables;
    out.constraints = solve.constraints;
    if (!solve.accepted) {
        out.status = "solver_rejected";
        return out;
    }
    out.accepted = true;
    out.stable = solve.objective >= -std::abs(stability_tolerance);
    out.trial_composition = solve.variables;
    out.status = "accepted";
    return out;
}

}  // namespace

StabilityNlpContract evaluate_neutral_stability_tpd_nlp_contract(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition,
    int parent_phase,
    int trial_phase
) {
    StabilityTpdProblem problem(
        args,
        temperature,
        pressure,
        feed_composition,
        parent_phase,
        trial_phase,
        "neutral_stability_tpd"
    );
    return make_contract(problem);
}

StabilityNlpContract evaluate_electrolyte_stability_tpd_nlp_contract(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition
) {
    StabilityTpdProblem problem(
        args,
        temperature,
        pressure,
        feed_composition,
        0,
        0,
        "electrolyte_stability_tpd",
        args.z,
        true
    );
    return make_contract(problem);
}

StabilityRouteResult solve_neutral_stability_tpd_route(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition,
    int parent_phase,
    int trial_phase,
    const IpoptSolveOptions& options,
    double stability_tolerance,
    const std::vector<double>& trial_initial_composition
) {
    const std::vector<NamedInitialComposition> seeds =
        stability_seed_candidates(feed_composition, {}, trial_initial_composition);
    StabilityRouteResult best;
    bool have_best = false;
    std::vector<RouteSeedAttempt> attempts;
    attempts.reserve(seeds.size() + (options.initial_variables.empty() ? 0 : 1));

    auto run_attempt = [&](
        const std::string& seed_name,
        const std::vector<double>& composition,
        const IpoptSolveOptions& attempt_options
    ) {
        StabilityRouteResult result = solve_stability_tpd_route(
            args,
            temperature,
            pressure,
            feed_composition,
            parent_phase,
            trial_phase,
            "neutral_stability_tpd",
            seed_name,
            {},
            false,
            attempt_options,
            stability_tolerance,
            composition
        );
        result.initial_point_strategy = "deterministic_multistart";
        attempts.push_back(stability_seed_attempt_from_result(result));
        if (!have_best || stability_attempt_better(result, best)) {
            best = result;
            have_best = true;
        }
        return result;
    };

    if (!options.initial_variables.empty()) {
        const StabilityRouteResult continuation =
            run_attempt("continuation_state", trial_initial_composition, options);
        if (continuation.accepted && !continuation.stable) {
            best.seed_attempts = attempts;
            return best;
        }
    }

    for (const auto& seed : seeds) {
        IpoptSolveOptions attempt_options = options;
        attempt_options.initial_variables.clear();
        attempt_options.initial_bound_lower_multipliers.clear();
        attempt_options.initial_bound_upper_multipliers.clear();
        attempt_options.initial_constraint_multipliers.clear();
        const StabilityRouteResult attempt = run_attempt(seed.seed_name, seed.composition, attempt_options);
        if (attempt.accepted && !attempt.stable) {
            break;
        }
    }

    best.initial_point_strategy = "deterministic_multistart";
    best.seed_attempts = attempts;
    return best;
}

StabilityRouteResult solve_electrolyte_stability_tpd_route(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition,
    const IpoptSolveOptions& options,
    double stability_tolerance,
    const std::vector<double>& trial_initial_composition
) {
    const std::vector<NamedInitialComposition> seeds =
        stability_seed_candidates(feed_composition, args.z, trial_initial_composition);
    StabilityRouteResult best;
    bool have_best = false;
    std::vector<RouteSeedAttempt> attempts;
    attempts.reserve(seeds.size() + (options.initial_variables.empty() ? 0 : 1));

    auto run_attempt = [&](
        const std::string& seed_name,
        const std::vector<double>& composition,
        const IpoptSolveOptions& attempt_options
    ) {
        StabilityRouteResult result = solve_stability_tpd_route(
            args,
            temperature,
            pressure,
            feed_composition,
            0,
            0,
            "electrolyte_stability_tpd",
            seed_name,
            args.z,
            true,
            attempt_options,
            stability_tolerance,
            composition
        );
        result.initial_point_strategy = "deterministic_multistart";
        attempts.push_back(stability_seed_attempt_from_result(result));
        if (!have_best || stability_attempt_better(result, best)) {
            best = result;
            have_best = true;
        }
        return result;
    };

    if (!options.initial_variables.empty()) {
        const StabilityRouteResult continuation =
            run_attempt("continuation_state", trial_initial_composition, options);
        if (continuation.accepted && !continuation.stable) {
            best.seed_attempts = attempts;
            return best;
        }
    }

    for (const auto& seed : seeds) {
        IpoptSolveOptions attempt_options = options;
        attempt_options.initial_variables.clear();
        attempt_options.initial_bound_lower_multipliers.clear();
        attempt_options.initial_bound_upper_multipliers.clear();
        attempt_options.initial_constraint_multipliers.clear();
        const StabilityRouteResult attempt = run_attempt(seed.seed_name, seed.composition, attempt_options);
        if (attempt.accepted && !attempt.stable) {
            break;
        }
    }

    best.initial_point_strategy = "deterministic_multistart";
    best.seed_attempts = attempts;
    return best;
}

}  // namespace epcsaft::native::equilibrium_nlp
