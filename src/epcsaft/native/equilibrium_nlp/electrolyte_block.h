#pragma once

#include <string>
#include <vector>

struct add_args;

namespace epcsaft::native::equilibrium_nlp {

struct PhaseChargeBlockResult {
    std::string block;
    std::string derivative_backend;
    int phase_count = 0;
    int species_count = 0;
    int local_variable_count = 0;
    std::vector<std::string> constraint_names;
    std::vector<double> residuals;
    int jacobian_rows = 0;
    int jacobian_cols = 0;
    std::vector<double> jacobian_row_major;
};

PhaseChargeBlockResult evaluate_phase_charge_block(
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& charges,
    int local_variable_count
);

struct ElectrolyteContributionBlockResult {
    std::string block;
    std::string value_backend;
    std::string term_basis;
    bool active = false;
    double temperature = 0.0;
    double density = 0.0;
    std::vector<double> composition;
    std::vector<double> amounts;
    std::vector<double> charges;
    double phase_charge_residual = 0.0;
    double ion_residual_helmholtz = 0.0;
    double born_residual_helmholtz = 0.0;
    double electrolyte_residual_helmholtz = 0.0;
    double total_residual_helmholtz = 0.0;
};

ElectrolyteContributionBlockResult evaluate_electrolyte_contribution_block(
    const add_args& args,
    double temperature,
    double density,
    const std::vector<double>& composition,
    const std::vector<double>& amounts
);

}  // namespace epcsaft::native::equilibrium_nlp
