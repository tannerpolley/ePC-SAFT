#pragma once

#include "ipopt_adapter.h"

#include <string>
#include <vector>

struct add_args;

namespace epcsaft::native::equilibrium_nlp {

struct NeutralTwoPhaseEosNlpContract {
    std::string problem_name;
    std::string derivative_backend;
    int phase_count = 0;
    int species_count = 0;
    int variable_count = 0;
    int constraint_count = 0;
    int jacobian_nonzero_count = 0;
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

struct NeutralTwoPhaseEosPostsolve {
    bool accepted = false;
    std::string rejection_reason;
    std::string derivative_backend;
    int phase_count = 0;
    int species_count = 0;
    double material_balance_norm = 0.0;
    double pressure_consistency_norm = 0.0;
    double chemical_potential_consistency_norm = 0.0;
    double phase_distance = 0.0;
    double objective = 0.0;
    std::vector<double> constraints;
    std::vector<double> phase_amount_totals;
    std::vector<double> phase_volumes;
    std::vector<std::vector<double>> phase_compositions;
};

struct NeutralTwoPhaseEosRouteResult {
    bool compiled = false;
    bool adapter_available = false;
    bool ran = false;
    bool solver_accepted = false;
    bool accepted = false;
    bool exact_gradient_required = true;
    bool exact_jacobian_required = true;
    std::string backend = "ipopt";
    std::string adapter_kind = "native_tnlp_adapter";
    std::string problem_name = "neutral_two_phase_eos";
    std::string derivative_backend = "analytic_cppad";
    std::string status = "not_started";
    std::string solver_status = "not_started";
    std::string application_status = "not_started";
    std::string hessian_strategy = "limited_memory";
    double objective = 0.0;
    std::vector<double> variables;
    std::vector<double> constraints;
    std::vector<std::vector<double>> phase_amounts;
    std::vector<double> phase_volumes;
    NeutralTwoPhaseEosPostsolve postsolve;
};

NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts
);

IpoptSolveResult solve_neutral_two_phase_eos_ipopt(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    const IpoptSolveOptions& options
);

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
);

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
    double phase_distance_tolerance
);

}  // namespace epcsaft::native::equilibrium_nlp
