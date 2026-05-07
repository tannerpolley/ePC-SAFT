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
};

struct ElectrolyteBubbleOptionsNative {
    double initial_pressure = 1.0e5;
    double min_pressure = 1.0;
    double max_pressure = 1.0e8;
    int max_iterations = 80;
    int max_vapor_iterations = 30;
    int max_bracket_expansions = 40;
    double tolerance = 1.0e-6;
    double vapor_tolerance = 1.0e-10;
    double pressure_factor = 2.0;
    double min_composition = 1.0e-14;
    double charge_tolerance = 1.0e-8;
    bool return_best_effort = false;
    std::vector<double> initial_y_vap;
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

struct StabilityTrialNative {
    std::string parent_phase;
    std::string trial_phase;
    std::string seed_name;
    std::vector<double> composition;
    double tpd = 0.0;
    int iterations = 0;
    bool converged = false;
    bool unstable = false;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, std::string> diagnostics_string;
};

struct StabilityResultNative {
    std::string backend = "neutral_tpd";
    std::string problem_kind = "stability";
    bool stable = true;
    double min_tpd = 0.0;
    std::string parent_phase = "liq";
    std::string trial_phase = "liq";
    std::vector<double> trial_composition;
    std::vector<StabilityTrialNative> trials;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
    std::map<std::string, std::vector<double>> diagnostics_vector;
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

StabilityResultNative neutral_stability_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& parent_phases,
    const std::vector<std::string>& trial_phases
);

StabilityResultNative electrolyte_stability_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& species
);

EquilibriumResultNative tp_flash_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options
);

EquilibriumResultNative lle_flash_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& initial_liq1 = {},
    const std::vector<double>& initial_liq2 = {},
    double initial_beta = 0.5,
    bool has_initial_phases = false
);

EquilibriumResultNative electrolyte_lle_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& species,
    const std::vector<double>& initial_aq = {},
    const std::vector<double>& initial_org = {},
    double initial_beta_org = 0.5,
    bool has_initial_phases = false
);

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

EquilibriumResultNative electrolyte_bubble_pressure_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    const std::vector<double>& x_liq,
    const ElectrolyteBubbleOptionsNative& options,
    const std::vector<std::string>& species,
    const std::vector<std::string>& vapor_species
);
