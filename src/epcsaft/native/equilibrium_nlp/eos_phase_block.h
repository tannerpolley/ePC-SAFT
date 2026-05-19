#pragma once

#include <string>
#include <vector>

#include "electrolyte_block.h"

struct add_args;

namespace epcsaft::native::equilibrium_nlp {

struct EosPhaseBlockResult {
    std::string block;
    std::string derivative_backend;
    std::vector<std::string> variable_names;
    std::vector<std::string> constraint_names;
    double temperature = 0.0;
    double target_pressure = 0.0;
    double gas_constant_temperature = 0.0;
    double total_amount = 0.0;
    double volume = 0.0;
    double density = 0.0;
    std::vector<double> composition;
    double residual_helmholtz = 0.0;
    double eos_pressure = 0.0;
    double pressure_density_derivative = 0.0;
    double compressibility_factor = 0.0;
    double ideal_helmholtz = 0.0;
    double residual_helmholtz_term = 0.0;
    double pressure_work = 0.0;
    ElectrolyteContributionBlockResult electrolyte_contribution;
    double objective = 0.0;
    std::vector<double> gradient;
    std::string objective_curvature_backend;
    int objective_curvature_rows = 0;
    int objective_curvature_cols = 0;
    std::vector<double> objective_curvature_row_major;
    std::string objective_third_derivative_backend;
    int objective_third_derivative_rows = 0;
    int objective_third_derivative_cols = 0;
    std::vector<double> objective_third_derivative_tensor_row_major;
    std::string pressure_hessian_backend;
    int pressure_hessian_rows = 0;
    int pressure_hessian_cols = 0;
    std::vector<double> pressure_hessian_row_major;
    double pressure_consistency_residual = 0.0;
    std::string constraint_jacobian_backend;
    int constraint_jacobian_rows = 0;
    int constraint_jacobian_cols = 0;
    std::vector<double> constraint_jacobian_row_major;
    std::string pressure_jacobian_backend;
    int pressure_jacobian_rows = 0;
    int pressure_jacobian_cols = 0;
    std::vector<double> pressure_jacobian_row_major;
};

struct EosPhaseSystemResult {
    std::string block;
    std::string derivative_backend;
    int phase_count = 0;
    int species_count = 0;
    std::vector<std::string> variable_names;
    std::vector<std::string> constraint_names;
    double temperature = 0.0;
    double target_pressure = 0.0;
    std::vector<double> feed_amounts;
    std::vector<EosPhaseBlockResult> phase_blocks;
    double objective = 0.0;
    double association_objective = 0.0;
    std::vector<double> phase_association_objectives;
    std::vector<double> gradient;
    std::string objective_hessian_backend;
    int objective_hessian_rows = 0;
    int objective_hessian_cols = 0;
    std::vector<double> objective_hessian_row_major;
    std::vector<double> constraints;
    std::vector<double> phase_charge_residuals;
    std::vector<double> phase_association_residuals;
    std::string constraint_jacobian_backend;
    int constraint_jacobian_rows = 0;
    int constraint_jacobian_cols = 0;
    std::vector<double> constraint_jacobian_row_major;
    std::string constraint_hessian_backend;
    int constraint_hessian_rows = 0;
    int constraint_hessian_cols = 0;
    std::vector<double> constraint_hessian_tensor_row_major;
    std::vector<bool> constraint_has_hessian;
};

EosPhaseBlockResult evaluate_eos_phase_block(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& amounts,
    double volume
);

EosPhaseSystemResult evaluate_eos_phase_system(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    const std::vector<double>& charges,
    const std::vector<std::vector<double>>& association_site_fractions = {},
    const std::vector<double>& association_delta_row_major = {}
);

}  // namespace epcsaft::native::equilibrium_nlp
