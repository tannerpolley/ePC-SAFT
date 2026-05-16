#pragma once

#include <string>
#include <vector>

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

}  // namespace epcsaft::native::equilibrium_nlp
