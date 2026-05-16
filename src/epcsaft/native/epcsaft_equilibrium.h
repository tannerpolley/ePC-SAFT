#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "epcsaft_electrolyte.h"

struct EquilibriumOptionsNative {
    int max_iterations = 80;
    double tolerance = 1.0e-6;
    double damping = 0.5;
    double min_composition = 1.0e-12;
    bool include_phase_diagnostics = false;
    bool stability_precheck = true;
    std::string density_diagnostics = "auto";
    bool experimental_coupled_density_lle = false;
    std::string jacobian_backend = "auto";
    double timeout_seconds = 0.0;
    int max_seed_attempts = 0;
    int max_density_failures = 0;
    int max_total_objective_evaluations = 0;
};

struct EquilibriumPhaseNative {
    std::string label;
    std::vector<double> composition;
    double density = 0.0;
    double temperature = 0.0;
    double pressure = 0.0;
    double phase_fraction = 0.0;
    std::vector<double> ln_fugacity_coefficient;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, std::string> diagnostics_string;
};

struct EquilibriumAttemptDiagnosticsNative {
    std::string seed_name;
    std::string rejection_reason;
    double beta_org = 0.0;
    double phase_distance = 0.0;
    double solver_residual_norm = 0.0;
    double material_balance_error = 0.0;
    double charge_balance_error = 0.0;
    double gibbs_delta = 0.0;
    int iterations = 0;
};

struct EquilibriumResultNative {
    std::string backend;
    std::string problem_kind;
    std::vector<EquilibriumPhaseNative> phases;
    bool stable = false;
    bool split_detected = false;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
    std::map<std::string, std::vector<double>> diagnostics_vector;
    std::vector<EquilibriumAttemptDiagnosticsNative> attempt_diagnostics;
    std::vector<DensitySolveDiagnostics> density_diagnostics;
};

struct ElectrolyteLLEResidualEvaluationNative {
    std::string variable_model = "ascani_transformed_salt_pairs";
    std::vector<double> variables;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> residual;
    std::vector<double> jacobian_row_major;
    int jacobian_rows = 0;
    int jacobian_cols = 0;
    std::vector<double> gradient;
    double objective = 0.0;
    std::vector<double> aq_composition;
    std::vector<double> org_composition;
    std::vector<double> aq_ln_fugacity_coefficient;
    std::vector<double> org_ln_fugacity_coefficient;
    double aq_density = 0.0;
    double org_density = 0.0;
    double phase_fraction_org = 0.0;
    double material_balance_error = 0.0;
    double charge_balance_error = 0.0;
    double phase_distance = 0.0;
    double gibbs_delta = 0.0;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
    std::map<std::string, std::vector<double>> diagnostics_vector;
};

struct ReactivePhaseResidualEvaluationNative {
    std::string variable_model = "log_phase_species_amounts";
    std::vector<double> variables;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> residual;
    std::vector<double> jacobian_row_major;
    int jacobian_rows = 0;
    int jacobian_cols = 0;
    std::vector<double> gradient;
    double objective = 0.0;
    std::vector<double> phase1_composition;
    std::vector<double> phase2_composition;
    std::vector<double> phase1_amounts;
    std::vector<double> phase2_amounts;
    std::vector<double> phase1_ln_fugacity_coefficient;
    std::vector<double> phase2_ln_fugacity_coefficient;
    double phase1_density = 0.0;
    double phase2_density = 0.0;
    double phase_fraction_phase2 = 0.0;
    std::vector<double> element_balance_residuals;
    std::vector<double> reaction_residuals_phase1;
    std::vector<double> reaction_residuals_phase2;
    std::vector<double> reaction_residuals_cross_phase;
    std::vector<double> neutral_phase_equilibrium_residuals;
    std::vector<double> ionic_equilibrium_residuals;
    std::vector<double> phase_charge_residuals;
    double phase_distance = 0.0;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
    std::map<std::string, std::vector<double>> diagnostics_vector;
};

ElectrolyteLLEResidualEvaluationNative evaluate_electrolyte_lle_residual_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& species,
    const std::vector<double>& variables = {},
    bool has_variables = false,
    const std::vector<double>& initial_aq = {},
    const std::vector<double>& initial_org = {},
    double initial_beta_org = 0.5,
    bool has_initial_phases = false
);

ReactivePhaseResidualEvaluationNative evaluate_reactive_phase_equilibrium_residual_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major = {},
    const std::vector<double>& variables = {},
    bool has_variables = false,
    const std::vector<double>& initial_phase1 = {},
    const std::vector<double>& initial_phase2 = {},
    double initial_phase_fraction_phase2 = 0.5,
    bool has_initial_phases = false
);
EquilibriumResultNative reactive_phase_equilibrium_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major = {},
    const std::vector<double>& initial_phase1 = {},
    const std::vector<double>& initial_phase2 = {},
    double initial_phase_fraction_phase2 = 0.5,
    bool has_initial_phases = false
);
