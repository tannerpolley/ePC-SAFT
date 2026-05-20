#pragma once

#include <map>
#include <memory>
#include <string>
#include <vector>

#include "epcsaft_electrolyte.h"
#include "equilibrium_nlp/route_builders.h"
#include "equilibrium_nlp/route_metadata.h"

struct EquilibriumOptionsNative {
    double min_composition = 1.0e-12;
    std::string jacobian_backend = "auto";
};

struct ElectrolyteLLEResidualEvaluationNative {
    std::string variable_model = "ascani_transformed_salt_pairs";
    std::string density_backend;
    std::vector<std::string> residual_families;
    std::vector<std::string> constraint_families;
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
    std::string density_backend;
    std::vector<std::string> residual_families;
    std::vector<std::string> constraint_families;
    std::vector<double> variables;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<double> residual;
    std::vector<double> jacobian_row_major;
    std::vector<double> residual_hessian_tensor_row_major;
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

namespace epcsaft::native {

using NativeRouteMetadata = epcsaft::native::equilibrium_nlp::RouteMetadata;

void apply_route_metadata(
    ::ElectrolyteLLEResidualEvaluationNative& out,
    const NativeRouteMetadata& metadata
);

void apply_route_metadata(
    ::ReactivePhaseResidualEvaluationNative& out,
    const NativeRouteMetadata& metadata
);

void apply_route_metadata(
    epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosNlpContract& out,
    const NativeRouteMetadata& metadata
);

void apply_route_metadata(
    epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosPostsolve& out,
    const NativeRouteMetadata& metadata
);

void apply_route_metadata(
    epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosRouteResult& out,
    const NativeRouteMetadata& metadata
);

void apply_route_metadata(
    epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosPostsolve& out,
    const NativeRouteMetadata& metadata
);

void apply_route_metadata(
    epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult& out,
    const NativeRouteMetadata& metadata
);

}  // namespace epcsaft::native

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

epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosNlpContract evaluate_electrolyte_lle_liquid_root_nlp_contract_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& species,
    double phase_distance_tolerance
);

epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosRouteResult solve_electrolyte_lle_liquid_root_route_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& equilibrium_options,
    const std::vector<std::string>& species,
    const epcsaft::native::equilibrium_nlp::IpoptSolveOptions& solve_options,
    double material_tolerance,
    double charge_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
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
    bool has_initial_phases = false,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase1_mixture = nullptr,
    const std::vector<int>& phase1_global_indices = {},
    const std::shared_ptr<ePCSAFTMixtureNative>& phase2_mixture = nullptr,
    const std::vector<int>& phase2_global_indices = {}
);

epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosNlpContract evaluate_reactive_phase_liquid_root_nlp_contract_native(
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
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    double phase_distance_tolerance
);

epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult solve_reactive_phase_liquid_root_route_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& equilibrium_options,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    const epcsaft::native::equilibrium_nlp::IpoptSolveOptions& solve_options,
    double conserved_balance_tolerance,
    double reaction_tolerance,
    double phase_equilibrium_tolerance,
    double phase_distance_tolerance
);

epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult solve_reactive_phase_liquid_root_phase_model_route_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase1_mixture,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase2_mixture,
    const std::vector<int>& phase1_global_indices,
    const std::vector<int>& phase2_global_indices,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& equilibrium_options,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    const epcsaft::native::equilibrium_nlp::IpoptSolveOptions& solve_options,
    double conserved_balance_tolerance,
    double reaction_tolerance,
    double phase_equilibrium_tolerance,
    double phase_distance_tolerance
);
