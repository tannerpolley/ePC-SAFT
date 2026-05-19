#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "ipopt_adapter.h"

class ePCSAFTMixtureNative;

struct ChemicalEquilibriumOptionsNative {
    int max_iterations = 50;
    double tolerance = 1.0e-8;
    double min_mole_fraction = 1.0e-14;
    std::string jacobian_backend = "auto";
    std::string solver_backend = "auto";
    std::string hessian_mode = "auto";
    int iteration_history_limit = 20;
    std::string linear_solver = "auto";
    double acceptable_tolerance = 0.0;
    double constraint_violation_tolerance = 0.0;
    double dual_infeasibility_tolerance = 0.0;
    double complementarity_tolerance = 0.0;
    std::vector<double> initial_variables;
    std::vector<double> initial_bound_lower_multipliers;
    std::vector<double> initial_bound_upper_multipliers;
    std::vector<double> initial_constraint_multipliers;
    std::string phase = "liq";
    std::string activity_output = "auto";
};

struct ChemicalEquilibriumResultNative {
    bool success = false;
    std::string message;
    std::vector<double> composition;
    std::vector<double> activity_coefficients;
    std::vector<double> mass_balance_residuals;
    double charge_residual = 0.0;
    std::vector<double> reaction_residuals;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
    std::map<std::string, std::vector<double>> diagnostics_vector;
    std::vector<double> continuation_variables;
    std::vector<double> continuation_bound_lower_multipliers;
    std::vector<double> continuation_bound_upper_multipliers;
    std::vector<double> continuation_constraint_multipliers;
    std::vector<epcsaft::native::equilibrium_nlp::IpoptIterationRecord> iteration_history;
};

struct ChemicalResidualEvaluationNative {
    std::string variable_model = "log_species_amounts";
    std::vector<double> variables;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> residual;
    std::vector<double> jacobian_row_major;
    int jacobian_rows = 0;
    int jacobian_cols = 0;
    std::vector<double> gradient;
    double objective = 0.0;
    std::vector<double> composition;
    std::vector<double> activity_coefficients;
    std::vector<double> mass_balance_residuals;
    double charge_residual = 0.0;
    std::vector<double> reaction_residuals;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
    std::map<std::string, std::vector<double>> diagnostics_vector;
};

ChemicalEquilibriumResultNative chemical_equilibrium_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& initial_x,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options
);

ChemicalResidualEvaluationNative evaluate_chemical_equilibrium_residual_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& initial_x,
    const std::vector<double>& variables,
    bool has_variables,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options
);
