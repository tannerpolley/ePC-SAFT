#pragma once

#include <string>
#include <vector>

namespace epcsaft::native::equilibrium_nlp {

struct AssociationMassActionBlockResult {
    std::string block;
    std::string derivative_backend;
    int site_count = 0;
    std::vector<std::string> constraint_names;
    std::vector<double> residuals;
    std::vector<double> density_derivative;
    int jacobian_rows = 0;
    int jacobian_cols = 0;
    std::vector<double> site_fraction_jacobian_row_major;
    std::vector<double> site_composition_jacobian_row_major;
};

AssociationMassActionBlockResult evaluate_association_mass_action_block(
    double density,
    const std::vector<double>& site_fractions,
    const std::vector<double>& site_composition,
    const std::vector<double>& delta_row_major
);

}  // namespace epcsaft::native::equilibrium_nlp
