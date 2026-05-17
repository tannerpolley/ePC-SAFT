#include "stability_route_builders.h"

#include "epcsaft_core_internal.h"
#include "nlp_problem.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>

namespace epcsaft::native::equilibrium_nlp {
namespace {

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

std::vector<double> shifted_composition(const std::vector<double>& composition) {
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
        const double value = composition[index] * (1.0 + 0.2 * direction);
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

class NeutralStabilityTpdProblem final : public NlpProblem {
public:
    NeutralStabilityTpdProblem(
        add_args args,
        double temperature,
        double pressure,
        std::vector<double> feed_composition,
        int parent_phase,
        int trial_phase
    )
        : args_(std::move(args)),
          temperature_(temperature),
          pressure_(pressure),
          feed_composition_(normalized_positive_values(feed_composition, "stability feed")),
          parent_phase_(parent_phase),
          trial_phase_(trial_phase),
          parent_phase_label_(phase_label(parent_phase)),
          trial_phase_label_(phase_label(trial_phase)) {
        require_positive_finite(temperature_, "stability temperature");
        require_positive_finite(pressure_, "stability pressure");
        species_count_ = static_cast<int>(feed_composition_.size());
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
        return "neutral_stability_tpd";
    }

    int variable_count() const override {
        return species_count_;
    }

    int constraint_count() const override {
        return 1;
    }

    int jacobian_nonzero_count() const override {
        return species_count_;
    }

    NlpBounds bounds() const override {
        NlpBounds out;
        out.variable_lower.assign(static_cast<std::size_t>(species_count_), minimum_composition_);
        out.variable_upper.assign(static_cast<std::size_t>(species_count_), 1.0);
        out.constraint_lower = {0.0};
        out.constraint_upper = {0.0};
        return out;
    }

    std::vector<double> initial_point() const override {
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
        return {std::accumulate(variables.begin(), variables.end(), 0.0) - 1.0};
    }

    NlpJacobianStructure jacobian_structure() const override {
        NlpJacobianStructure out;
        out.rows.assign(static_cast<std::size_t>(species_count_), 0);
        out.cols.reserve(static_cast<std::size_t>(species_count_));
        for (int col = 0; col < species_count_; ++col) {
            out.cols.push_back(col);
        }
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        require_trial_variables(variables);
        return std::vector<double>(static_cast<std::size_t>(species_count_), 1.0);
    }

    NlpScaling scaling() const override {
        NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(static_cast<std::size_t>(species_count_), 1.0);
        out.constraints = {1.0};
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
    std::vector<double> feed_composition_;
    int parent_phase_ = 0;
    int trial_phase_ = 0;
    std::string parent_phase_label_;
    std::string trial_phase_label_;
    PhaseStateCompositionSensitivityResult parent_state_;
    std::vector<double> parent_reduced_potential_;
    int species_count_ = 0;
};

StabilityNlpContract make_contract(const NeutralStabilityTpdProblem& problem) {
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

}  // namespace

StabilityNlpContract evaluate_neutral_stability_tpd_nlp_contract(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition,
    int parent_phase,
    int trial_phase
) {
    NeutralStabilityTpdProblem problem(args, temperature, pressure, feed_composition, parent_phase, trial_phase);
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
    double stability_tolerance
) {
    const IpoptAdapterInfo adapter = native_ipopt_adapter_info();
    StabilityRouteResult out;
    out.compiled = adapter.compiled;
    out.adapter_available = adapter.adapter_available;
    out.adapter_kind = adapter.adapter_kind;
    out.exact_gradient_required = adapter.exact_gradient_required;
    out.exact_jacobian_required = adapter.exact_jacobian_required;
    out.parent_phase = phase_label(parent_phase);
    out.trial_phase = phase_label(trial_phase);
    if (!adapter.compiled) {
        out.status = "ipopt_dependency_required";
        return out;
    }

    NeutralStabilityTpdProblem problem(args, temperature, pressure, feed_composition, parent_phase, trial_phase);
    out.parent_reduced_potential = problem.parent_reduced_potential();
    const IpoptSolveResult solve = solve_ipopt_nlp(problem, options);
    out.ran = solve.solver_ran;
    out.solver_accepted = solve.accepted;
    out.solver_status = solve.solver_status;
    out.application_status = solve.application_status;
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

}  // namespace epcsaft::native::equilibrium_nlp
