#pragma once

#include "ipopt_adapter.h"

#include <string>
#include <vector>

struct add_args;

namespace epcsaft::native::equilibrium_nlp {

struct NeutralTwoPhaseEosNlpContract {
    std::string problem_name;
    std::string derivative_backend;
    std::string variable_model;
    std::string density_backend;
    int phase_count = 0;
    int species_count = 0;
    int balance_row_count = 0;
    int reaction_count = 0;
    int variable_count = 0;
    int constraint_count = 0;
    int jacobian_nonzero_count = 0;
    std::vector<double> standard_mu_rt;
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
    double ln_fugacity_consistency_norm = 0.0;
    double charge_balance_norm = 0.0;
    double fixed_composition_norm = 0.0;
    double phase_amount_total_norm = 0.0;
    double phase_distance = 0.0;
    double objective = 0.0;
    double gibbs_feed = 0.0;
    double gibbs_split = 0.0;
    double gibbs_delta = 0.0;
    double minimum_phase_fraction = 0.0;
    std::string density_backend;
    std::vector<double> constraints;
    std::vector<double> phase_amount_totals;
    std::vector<double> phase_volumes;
    std::vector<double> phase_densities;
    std::vector<std::vector<double>> phase_compositions;
    std::vector<std::vector<double>> phase_ln_fugacity_coefficients;
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
    std::string gradient_approximation = "exact";
    std::string jacobian_approximation = "exact";
    std::string hessian_approximation = "limited-memory";
    std::string hessian_backend = "limited-memory";
    std::string scaling_method = "user-scaling";
    std::string status;
    std::string solver_status;
    std::string application_status;
    std::string last_callback_exception;
    int iteration_count = 0;
    int iteration_history_limit = 0;
    int iteration_history_size = 0;
    int variable_scaling_count = 0;
    int constraint_scaling_count = 0;
    int eval_h_calls = 0;
    double objective_scaling = 1.0;
    double variable_scaling_min = 0.0;
    double variable_scaling_max = 0.0;
    double constraint_scaling_min = 0.0;
    double constraint_scaling_max = 0.0;
    bool exact_hessian_available = false;
    bool warm_start_requested = false;
    bool warm_start_used = false;
    double objective = 0.0;
    std::vector<double> variables;
    std::vector<double> constraints;
    std::vector<double> bound_lower_multipliers;
    std::vector<double> bound_upper_multipliers;
    std::vector<double> constraint_multipliers;
    std::vector<IpoptIterationRecord> iteration_history;
    std::vector<std::vector<double>> phase_amounts;
    std::vector<double> phase_volumes;
    NeutralTwoPhaseEosPostsolve postsolve;
};

struct ReactiveTwoPhaseEosPostsolve {
    bool accepted = false;
    std::string rejection_reason;
    std::string derivative_backend = "analytic_cppad";
    std::string density_backend;
    int phase_count = 0;
    int species_count = 0;
    int balance_row_count = 0;
    int reaction_count = 0;
    double conserved_balance_norm = 0.0;
    double charge_balance_norm = 0.0;
    double pressure_consistency_norm = 0.0;
    double phase_equilibrium_norm = 0.0;
    double reaction_stationarity_norm = 0.0;
    double phase_distance = 0.0;
    double objective = 0.0;
    std::vector<double> standard_mu_rt;
    std::vector<double> constraints;
    std::vector<double> reaction_stationarity_residuals;
    std::vector<double> phase_equilibrium_residuals;
    std::vector<double> phase_charge_residuals;
    std::vector<double> phase_amount_totals;
    std::vector<double> phase_volumes;
    std::vector<double> phase_densities;
    std::vector<std::vector<double>> phase_compositions;
    std::vector<std::vector<double>> phase_ln_fugacity_coefficients;
};

struct ReactiveTwoPhaseEosRouteResult {
    bool compiled = false;
    bool adapter_available = false;
    bool ran = false;
    bool solver_accepted = false;
    bool accepted = false;
    bool exact_gradient_required = true;
    bool exact_jacobian_required = true;
    std::string backend = "ipopt";
    std::string adapter_kind = "native_tnlp_adapter";
    std::string problem_name = "reactive_two_phase_eos";
    std::string derivative_backend = "analytic_cppad";
    std::string gradient_approximation = "exact";
    std::string jacobian_approximation = "exact";
    std::string hessian_approximation = "limited-memory";
    std::string hessian_backend = "limited-memory";
    std::string scaling_method = "user-scaling";
    std::string status;
    std::string solver_status;
    std::string application_status;
    int iteration_count = 0;
    int iteration_history_limit = 0;
    int iteration_history_size = 0;
    int variable_scaling_count = 0;
    int constraint_scaling_count = 0;
    int eval_h_calls = 0;
    double objective_scaling = 1.0;
    double variable_scaling_min = 0.0;
    double variable_scaling_max = 0.0;
    double constraint_scaling_min = 0.0;
    double constraint_scaling_max = 0.0;
    bool exact_hessian_available = false;
    bool warm_start_requested = false;
    bool warm_start_used = false;
    int phase_count = 0;
    int species_count = 0;
    int balance_row_count = 0;
    int reaction_count = 0;
    double objective = 0.0;
    std::vector<double> standard_mu_rt;
    std::vector<double> variables;
    std::vector<double> constraints;
    std::vector<double> bound_lower_multipliers;
    std::vector<double> bound_upper_multipliers;
    std::vector<double> constraint_multipliers;
    std::vector<IpoptIterationRecord> iteration_history;
    std::vector<std::vector<double>> phase_amounts;
    std::vector<double> phase_volumes;
    ReactiveTwoPhaseEosPostsolve postsolve;
};

void apply_ipopt_solve_metadata(NeutralTwoPhaseEosRouteResult& out, const IpoptSolveResult& solve);

void apply_ipopt_solve_metadata(ReactiveTwoPhaseEosRouteResult& out, const IpoptSolveResult& solve);

NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts
);

NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_tp_flash_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts
);

NeutralTwoPhaseEosNlpContract evaluate_neutral_two_phase_eos_lle_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts
);

NeutralTwoPhaseEosNlpContract evaluate_electrolyte_lle_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts
);

NeutralTwoPhaseEosNlpContract evaluate_reactive_two_phase_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    int balance_rows,
    const std::vector<double>& balance_matrix_row_major,
    const std::vector<double>& total_vector,
    int reaction_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants
);

ReactiveTwoPhaseEosRouteResult solve_reactive_two_phase_eos_route(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    int balance_rows,
    const std::vector<double>& balance_matrix_row_major,
    const std::vector<double>& total_vector,
    int reaction_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants,
    const IpoptSolveOptions& options,
    double conserved_balance_tolerance,
    double pressure_tolerance,
    double reaction_stationarity_tolerance,
    double phase_distance_tolerance,
    const std::vector<double>& charges
);

NeutralTwoPhaseEosNlpContract evaluate_reactive_lle_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts,
    int balance_rows,
    const std::vector<double>& balance_matrix_row_major,
    const std::vector<double>& total_vector,
    int reaction_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants
);

NeutralTwoPhaseEosNlpContract evaluate_reactive_electrolyte_lle_eos_nlp_contract(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts,
    int balance_rows,
    const std::vector<double>& balance_matrix_row_major,
    const std::vector<double>& total_vector,
    int reaction_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants
);

ReactiveTwoPhaseEosRouteResult solve_reactive_lle_eos_route(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts,
    int balance_rows,
    const std::vector<double>& balance_matrix_row_major,
    const std::vector<double>& total_vector,
    int reaction_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants,
    const IpoptSolveOptions& options,
    double conserved_balance_tolerance,
    double pressure_tolerance,
    double reaction_stationarity_tolerance,
    double phase_distance_tolerance
);

ReactiveTwoPhaseEosRouteResult solve_reactive_electrolyte_lle_eos_route(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& feed_amounts,
    int balance_rows,
    const std::vector<double>& balance_matrix_row_major,
    const std::vector<double>& total_vector,
    int reaction_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants,
    const IpoptSolveOptions& options,
    double conserved_balance_tolerance,
    double pressure_tolerance,
    double reaction_stationarity_tolerance,
    double phase_distance_tolerance
);

ReactiveTwoPhaseEosPostsolve evaluate_reactive_two_phase_eos_postsolve(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    int balance_rows,
    const std::vector<double>& balance_matrix_row_major,
    const std::vector<double>& total_vector,
    int reaction_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants,
    double conserved_balance_tolerance,
    double pressure_tolerance,
    double reaction_stationarity_tolerance,
    double phase_distance_tolerance,
    const std::vector<double>& charges
);

NeutralTwoPhaseEosNlpContract evaluate_neutral_bubble_p_eos_nlp_contract(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition
);

NeutralTwoPhaseEosNlpContract evaluate_electrolyte_bubble_p_eos_nlp_contract(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition
);

NeutralTwoPhaseEosNlpContract evaluate_neutral_dew_p_eos_nlp_contract(
    const add_args& args,
    double temperature,
    const std::vector<double>& vapor_composition
);

NeutralTwoPhaseEosNlpContract evaluate_neutral_bubble_t_eos_nlp_contract(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& liquid_composition
);

NeutralTwoPhaseEosNlpContract evaluate_neutral_dew_t_eos_nlp_contract(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& vapor_composition
);

IpoptSolveResult solve_neutral_two_phase_eos_ipopt(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    const IpoptSolveOptions& options,
    const std::vector<double>& charges = {},
    const std::string& problem_name = "neutral_two_phase_eos",
    double minimum_phase_distance = 0.0
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
    double phase_distance_tolerance,
    const std::vector<double>& charges = {}
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
    double phase_distance_tolerance,
    double minimum_phase_distance = 0.0,
    const std::vector<double>& charges = {},
    const std::string& problem_name = "neutral_two_phase_eos",
    double charge_tolerance = 0.0
);

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
);

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
);

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
);

NeutralTwoPhaseEosRouteResult solve_neutral_bubble_p_eos_route(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
);

NeutralTwoPhaseEosRouteResult solve_electrolyte_bubble_p_eos_route(
    const add_args& args,
    double temperature,
    const std::vector<double>& liquid_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double charge_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
);

NeutralTwoPhaseEosRouteResult solve_neutral_dew_p_eos_route(
    const add_args& args,
    double temperature,
    const std::vector<double>& vapor_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
);

NeutralTwoPhaseEosRouteResult solve_neutral_bubble_t_eos_route(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& liquid_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
);

NeutralTwoPhaseEosRouteResult solve_neutral_dew_t_eos_route(
    const add_args& args,
    double target_pressure,
    const std::vector<double>& vapor_composition,
    const IpoptSolveOptions& options,
    double phase_total_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
);

}  // namespace epcsaft::native::equilibrium_nlp
