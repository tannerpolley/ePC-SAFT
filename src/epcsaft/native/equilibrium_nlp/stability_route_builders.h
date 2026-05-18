#pragma once

#include "ipopt_adapter.h"

#include <string>
#include <vector>

struct add_args;

namespace epcsaft::native::equilibrium_nlp {

struct StabilityNlpContract {
    std::string problem_name;
    std::string derivative_backend;
    int species_count = 0;
    int variable_count = 0;
    int constraint_count = 0;
    int jacobian_nonzero_count = 0;
    std::string parent_phase;
    std::string trial_phase;
    std::vector<double> feed_composition;
    std::vector<double> parent_reduced_potential;
    std::vector<double> initial_point;
    std::vector<double> variable_lower_bounds;
    std::vector<double> variable_upper_bounds;
    std::vector<double> constraint_lower_bounds;
    std::vector<double> constraint_upper_bounds;
    double objective_at_initial = 0.0;
    std::vector<double> gradient_at_initial;
    std::vector<double> constraints_at_initial;
    std::vector<int> jacobian_rows;
    std::vector<int> jacobian_cols;
    std::vector<double> jacobian_values_at_initial;
};

struct StabilityRouteResult {
    bool compiled = false;
    bool adapter_available = false;
    bool ran = false;
    bool solver_accepted = false;
    bool accepted = false;
    bool stable = false;
    bool exact_gradient_required = true;
    bool exact_jacobian_required = true;
    std::string backend = "ipopt";
    std::string adapter_kind = "native_tnlp_adapter";
    std::string problem_name = "neutral_stability_tpd";
    std::string derivative_backend = "cppad_implicit";
    std::string status;
    std::string solver_status;
    std::string application_status;
    std::string parent_phase;
    std::string trial_phase;
    std::string seed_name = "canonical_shifted_feed";
    double objective = 0.0;
    double min_tpd = 0.0;
    std::vector<double> variables;
    std::vector<double> constraints;
    std::vector<double> trial_composition;
    std::vector<double> initial_composition;
    std::vector<double> parent_reduced_potential;
};

StabilityNlpContract evaluate_neutral_stability_tpd_nlp_contract(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition,
    int parent_phase,
    int trial_phase
);

StabilityNlpContract evaluate_electrolyte_stability_tpd_nlp_contract(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition
);

StabilityRouteResult solve_neutral_stability_tpd_route(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition,
    int parent_phase,
    int trial_phase,
    const IpoptSolveOptions& options,
    double stability_tolerance,
    const std::vector<double>& trial_initial_composition = {}
);

StabilityRouteResult solve_electrolyte_stability_tpd_route(
    const add_args& args,
    double temperature,
    double pressure,
    const std::vector<double>& feed_composition,
    const IpoptSolveOptions& options,
    double stability_tolerance,
    const std::vector<double>& trial_initial_composition = {}
);

}  // namespace epcsaft::native::equilibrium_nlp
