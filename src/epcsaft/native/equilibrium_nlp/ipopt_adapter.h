#pragma once

#include "nlp_problem.h"

#include <map>
#include <string>
#include <vector>

namespace epcsaft::native::equilibrium_nlp {

struct IpoptAdapterInfo {
    bool compiled = false;
    bool adapter_available = false;
    std::string backend = "ipopt";
    std::string status = "disabled";
    std::string adapter_kind = "native_tnlp_adapter";
    bool exact_gradient_required = true;
    bool exact_jacobian_required = true;
};

struct IpoptSolveOptions {
    int max_iterations = 100;
    int print_level = 0;
    double tolerance = 1.0e-8;
    double acceptable_tolerance = 0.0;
    double constraint_violation_tolerance = 0.0;
    double dual_infeasibility_tolerance = 0.0;
    double complementarity_tolerance = 0.0;
    double max_wall_time_seconds = 0.0;
    bool limited_memory_hessian = true;
    std::string hessian_mode = "auto";
    int iteration_history_limit = 20;
    std::string linear_solver = "auto";
    std::vector<double> initial_variables;
    std::vector<double> initial_bound_lower_multipliers;
    std::vector<double> initial_bound_upper_multipliers;
    std::vector<double> initial_constraint_multipliers;
};

struct IpoptIterationRecord {
    int iteration = 0;
    double objective = 0.0;
    double primal_infeasibility = 0.0;
    double dual_infeasibility = 0.0;
    double barrier_parameter = 0.0;
    double step_size_primal = 0.0;
    double step_size_dual = 0.0;
    double regularization_size = 0.0;
    double complementarity = 0.0;
    int step_trial_count = 0;
    bool restoration_phase = false;
};

struct RouteSeedAttempt {
    std::string seed_name;
    std::string status;
    std::string solver_status;
    std::string application_status;
    bool solver_accepted = false;
    bool accepted = false;
    bool stable = false;
    int iteration_count = 0;
    double objective = 0.0;
    double phase_distance = 0.0;
    double material_balance_norm = 0.0;
    double conserved_balance_norm = 0.0;
    double charge_balance_norm = 0.0;
    double pressure_consistency_norm = 0.0;
    double chemical_potential_consistency_norm = 0.0;
    double phase_equilibrium_norm = 0.0;
    double reaction_stationarity_norm = 0.0;
    double min_tpd = 0.0;
};

struct IpoptSolveResult {
    bool solver_ran = false;
    bool solved = false;
    bool acceptable = false;
    bool accepted = false;
    std::string backend = "ipopt";
    std::string adapter_kind = "native_tnlp_adapter";
    std::string solver_status;
    std::string application_status;
    double objective = 0.0;
    std::vector<double> variables;
    std::vector<double> constraints;
    std::vector<double> bound_lower_multipliers;
    std::vector<double> bound_upper_multipliers;
    std::vector<double> constraint_multipliers;
    std::vector<IpoptIterationRecord> iteration_history;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
};

IpoptAdapterInfo native_ipopt_adapter_info();

IpoptSolveResult solve_ipopt_nlp(
    const NlpProblem& problem,
    const IpoptSolveOptions& options
);

IpoptSolveResult solve_ipopt_quadratic_smoke();

IpoptSolveResult solve_ipopt_quadratic_smoke(const IpoptSolveOptions& options);

}  // namespace epcsaft::native::equilibrium_nlp
