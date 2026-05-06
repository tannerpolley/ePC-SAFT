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
