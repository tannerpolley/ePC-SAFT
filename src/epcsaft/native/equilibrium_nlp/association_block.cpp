#include "association_block.h"

#include "epcsaft_electrolyte.h"

#include <cmath>

namespace epcsaft::native::equilibrium_nlp {

namespace {

void require_finite(double value, const std::string& label) {
    if (std::isfinite(value)) {
        return;
    }
    throw ValueError(label + " must be finite.");
}

void require_positive_finite(double value, const std::string& label) {
    if (std::isfinite(value) && value > 0.0) {
        return;
    }
    throw ValueError(label + " must be positive and finite.");
}

}  // namespace

AssociationMassActionBlockResult evaluate_association_mass_action_block(
    double density,
    const std::vector<double>& site_fractions,
    const std::vector<double>& site_composition,
    const std::vector<double>& delta_row_major
) {
    require_finite(density, "Association block density");
    if (density < 0.0) {
        throw ValueError("Association block density must be non-negative.");
    }
    const std::size_t site_count = site_fractions.size();
    if (site_count == 0) {
        throw ValueError("Association block requires at least one site.");
    }
    if (site_composition.size() != site_count) {
        throw ValueError("Association block site composition size must match site fraction size.");
    }
    if (delta_row_major.size() != site_count * site_count) {
        throw ValueError("Association block delta matrix must be square with one row per site.");
    }

    AssociationMassActionBlockResult result;
    result.block = "association_mass_action";
    result.derivative_backend = "analytic";
    result.site_count = static_cast<int>(site_count);
    result.constraint_names.reserve(site_count);
    result.residuals.assign(site_count, 0.0);
    result.density_derivative.assign(site_count, 0.0);
    result.jacobian_rows = static_cast<int>(site_count);
    result.jacobian_cols = static_cast<int>(site_count);
    result.site_fraction_jacobian_row_major.assign(site_count * site_count, 0.0);
    result.site_composition_jacobian_row_major.assign(site_count * site_count, 0.0);

    for (std::size_t site = 0; site < site_count; ++site) {
        require_positive_finite(site_fractions[site], "Association site fraction");
        require_finite(site_composition[site], "Association site composition");
        if (site_composition[site] < 0.0) {
            throw ValueError("Association site composition must be non-negative.");
        }
        result.constraint_names.push_back("association_site_" + std::to_string(site));
        result.residuals[site] = 1.0 / site_fractions[site] - 1.0;
        result.site_fraction_jacobian_row_major[site * site_count + site] =
            -1.0 / (site_fractions[site] * site_fractions[site]);
    }

    for (std::size_t row = 0; row < site_count; ++row) {
        for (std::size_t col = 0; col < site_count; ++col) {
            const std::size_t index = row * site_count + col;
            require_finite(delta_row_major[index], "Association delta matrix entry");
            result.residuals[row] -= density
                * site_composition[col]
                * site_fractions[col]
                * delta_row_major[index];
            result.density_derivative[row] -= site_composition[col]
                * site_fractions[col]
                * delta_row_major[index];
            result.site_fraction_jacobian_row_major[index] -= density
                * site_composition[col]
                * delta_row_major[index];
            result.site_composition_jacobian_row_major[index] = -density
                * site_fractions[col]
                * delta_row_major[index];
        }
    }
    return result;
}

}  // namespace epcsaft::native::equilibrium_nlp
