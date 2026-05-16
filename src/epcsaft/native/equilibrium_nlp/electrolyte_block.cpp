#include "electrolyte_block.h"

#include "epcsaft_core_internal.h"

#include <algorithm>
#include <cmath>
#include <numeric>

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

ElectrolyteContributionBlockResult evaluate_electrolyte_contribution_block(
    const add_args& args,
    double temperature,
    double density,
    const std::vector<double>& composition,
    const std::vector<double>& amounts
) {
    require_finite(temperature, "Electrolyte contribution temperature");
    require_finite(density, "Electrolyte contribution density");
    if (temperature <= 0.0) {
        throw ValueError("Electrolyte contribution temperature must be positive.");
    }
    if (density <= 0.0) {
        throw ValueError("Electrolyte contribution density must be positive.");
    }
    if (composition.empty()) {
        throw ValueError("Electrolyte contribution block requires at least one composition entry.");
    }
    if (!amounts.empty() && amounts.size() != composition.size()) {
        throw ValueError("Electrolyte contribution amount count must match composition size.");
    }
    for (double value : composition) {
        require_finite(value, "Electrolyte contribution composition");
        if (value < 0.0) {
            throw ValueError("Electrolyte contribution composition entries must be non-negative.");
        }
    }
    const double composition_sum = std::accumulate(composition.begin(), composition.end(), 0.0);
    if (!(std::isfinite(composition_sum) && composition_sum > 0.0)) {
        throw ValueError("Electrolyte contribution composition sum must be positive.");
    }
    if (!args.z.empty() && args.z.size() != composition.size()) {
        throw ValueError("Electrolyte contribution charge count must match composition size.");
    }

    ElectrolyteContributionBlockResult result;
    result.block = "electrolyte_contribution";
    result.value_backend = "native_eos";
    result.term_basis = "dimensionless_residual_helmholtz";
    result.temperature = temperature;
    result.density = density;
    result.composition = composition;
    result.amounts = amounts;
    result.charges = args.z;

    if (!args.z.empty()) {
        const std::vector<double>& charge_basis = amounts.empty() ? composition : amounts;
        for (std::size_t index = 0; index < args.z.size(); ++index) {
            result.phase_charge_residual += args.z[index] * charge_basis[index];
        }
    }

    const ScalarContributionTerms terms = residual_helmholtz_result_cpp(
        temperature,
        density,
        composition,
        args
    );
    result.ion_residual_helmholtz = terms.ion;
    result.born_residual_helmholtz = terms.born;
    result.electrolyte_residual_helmholtz = terms.ion + terms.born;
    result.total_residual_helmholtz = terms.total;
    result.active = std::any_of(args.z.begin(), args.z.end(), [](double charge) {
        return std::isfinite(charge) && std::abs(charge) > 0.0;
    }) || std::abs(result.ion_residual_helmholtz) > 0.0
        || std::abs(result.born_residual_helmholtz) > 0.0;
    return result;
}

}  // namespace epcsaft::native::equilibrium_nlp
