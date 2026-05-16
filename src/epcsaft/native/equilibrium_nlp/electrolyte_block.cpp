#include "electrolyte_block.h"

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

}  // namespace

PhaseChargeBlockResult evaluate_phase_charge_block(
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& charges,
    int local_variable_count
) {
    if (phase_amounts.empty()) {
        throw ValueError("Phase charge block requires at least one phase.");
    }
    if (charges.empty()) {
        throw ValueError("Phase charge block requires at least one charge.");
    }
    if (local_variable_count <= 0) {
        throw ValueError("Phase charge block requires a positive local variable count.");
    }
    const std::size_t phase_count = phase_amounts.size();
    const std::size_t species_count = charges.size();
    if (static_cast<std::size_t>(local_variable_count) < species_count) {
        throw ValueError("Phase charge block local variable count cannot be smaller than species count.");
    }
    for (double charge : charges) {
        require_finite(charge, "Phase charge");
    }

    PhaseChargeBlockResult result;
    result.block = "phase_charge";
    result.derivative_backend = "analytic";
    result.phase_count = static_cast<int>(phase_count);
    result.species_count = static_cast<int>(species_count);
    result.local_variable_count = local_variable_count;
    result.constraint_names.reserve(phase_count);
    result.residuals.assign(phase_count, 0.0);
    const std::size_t variable_count = phase_count * static_cast<std::size_t>(local_variable_count);
    result.jacobian_rows = static_cast<int>(phase_count);
    result.jacobian_cols = static_cast<int>(variable_count);
    result.jacobian_row_major.assign(phase_count * variable_count, 0.0);

    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        if (phase_amounts[phase].size() != species_count) {
            throw ValueError("Phase charge block amount sizes must match charge count.");
        }
        result.constraint_names.push_back("phase_" + std::to_string(phase) + ".charge_balance");
        const std::size_t column_offset = phase * static_cast<std::size_t>(local_variable_count);
        for (std::size_t species = 0; species < species_count; ++species) {
            require_finite(phase_amounts[phase][species], "Phase charge amount");
            result.residuals[phase] += charges[species] * phase_amounts[phase][species];
            result.jacobian_row_major[phase * variable_count + column_offset + species] = charges[species];
        }
    }
    return result;
}

}  // namespace epcsaft::native::equilibrium_nlp
