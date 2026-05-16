#include "eos_phase_block.h"

#include "electrolyte_block.h"
#include "epcsaft_core_internal.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>

namespace epcsaft::native::equilibrium_nlp {

namespace {

void require_positive_finite(double value, const std::string& label) {
    if (std::isfinite(value) && value > 0.0) {
        return;
    }
    throw ValueError(label + " must be positive and finite.");
}

void validate_amounts(const std::vector<double>& amounts) {
    if (amounts.empty()) {
        throw ValueError("EOS phase block requires at least one species amount.");
    }
    for (double amount : amounts) {
        require_positive_finite(amount, "EOS phase species amount");
    }
}

std::vector<std::string> amount_variable_names(std::size_t count) {
    std::vector<std::string> names;
    names.reserve(count + 1);
    for (std::size_t index = 0; index < count; ++index) {
        names.push_back("n_" + std::to_string(index));
    }
    names.push_back("volume");
    return names;
}

std::vector<std::string> phase_system_variable_names(std::size_t phase_count, std::size_t species_count) {
    std::vector<std::string> names;
    names.reserve(phase_count * (species_count + 1));
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        for (std::size_t species = 0; species < species_count; ++species) {
            names.push_back("phase_" + std::to_string(phase) + ".n_" + std::to_string(species));
        }
        names.push_back("phase_" + std::to_string(phase) + ".volume");
    }
    return names;
}

std::vector<std::string> phase_system_constraint_names(std::size_t phase_count, std::size_t species_count) {
    std::vector<std::string> names;
    names.reserve(species_count + phase_count);
    for (std::size_t species = 0; species < species_count; ++species) {
        names.push_back("material_balance_" + std::to_string(species));
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        names.push_back("phase_" + std::to_string(phase) + ".pressure_consistency");
    }
    return names;
}

std::vector<double> composition_from_amounts(const std::vector<double>& amounts, double total_amount) {
    std::vector<double> composition;
    composition.reserve(amounts.size());
    for (double amount : amounts) {
        composition.push_back(amount / total_amount);
    }
    return composition;
}

double ideal_helmholtz_amount_volume_term(const std::vector<double>& amounts, double volume) {
    double value = 0.0;
    for (double amount : amounts) {
        value += amount * (std::log(amount / volume) - 1.0);
    }
    return value;
}

void require_finite_nonnegative(double value, const std::string& label) {
    if (std::isfinite(value) && value >= 0.0) {
        return;
    }
    throw ValueError(label + " must be finite and non-negative.");
}

}  // namespace

EosPhaseBlockResult evaluate_eos_phase_block(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<double>& amounts,
    double volume
) {
    require_positive_finite(temperature, "EOS phase temperature");
    if (!std::isfinite(target_pressure)) {
        throw ValueError("EOS phase target pressure must be finite.");
    }
    require_positive_finite(volume, "EOS phase volume");
    validate_amounts(amounts);

    const double total_amount = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    require_positive_finite(total_amount, "EOS phase total amount");
    const std::vector<double> composition = composition_from_amounts(amounts, total_amount);
    if (!args.m.empty() && args.m.size() != composition.size()) {
        std::ostringstream msg;
        msg << "EOS phase block composition size " << composition.size()
            << " does not match parameter size " << args.m.size() << ".";
        throw ValueError(msg.str());
    }

    const double density = total_amount / volume;
    require_positive_finite(density, "EOS phase molar density");
    const double rt = kb * N_AV * temperature;
    const ScalarContributionTerms helmholtz = residual_helmholtz_result_cpp(
        temperature,
        density,
        composition,
        args
    );
    const CompressibilityFactorResult z = compressibility_factor_result_cpp(
        temperature,
        density,
        composition,
        args
    );
    const ResidualChemicalPotentialResult mu = residual_chemical_potential_result_cpp(
        temperature,
        density,
        composition,
        args
    );
    if (mu.mu.total.size() != amounts.size()) {
        throw ValueError("EOS phase residual chemical potential size did not match species count.");
    }
    const double eos_pressure = p_cpp(temperature, density, composition, args);

    EosPhaseBlockResult result;
    result.block = "eos_phase";
    result.derivative_backend = "analytic";
    result.variable_names = amount_variable_names(amounts.size());
    result.constraint_names = {"pressure_consistency"};
    result.temperature = temperature;
    result.target_pressure = target_pressure;
    result.gas_constant_temperature = rt;
    result.total_amount = total_amount;
    result.volume = volume;
    result.density = density;
    result.composition = composition;
    result.residual_helmholtz = helmholtz.total;
    result.eos_pressure = eos_pressure;
    result.compressibility_factor = z.terms.total;
    result.ideal_helmholtz = ideal_helmholtz_amount_volume_term(amounts, volume);
    result.residual_helmholtz_term = total_amount * helmholtz.total;
    result.pressure_work = target_pressure * volume / rt;
    result.objective = result.ideal_helmholtz + result.residual_helmholtz_term + result.pressure_work;
    result.pressure_consistency_residual = eos_pressure - target_pressure;
    result.gradient.reserve(amounts.size() + 1);
    for (std::size_t index = 0; index < amounts.size(); ++index) {
        result.gradient.push_back(std::log(amounts[index] / volume) + mu.mu.total[index]);
    }
    result.gradient.push_back((target_pressure - eos_pressure) / rt);
    double cppad_objective = 0.0;
    std::vector<double> cppad_gradient;
    std::vector<double> cppad_hessian;
    eos_phase_objective_derivatives_cpp(
        temperature,
        target_pressure,
        amounts,
        volume,
        args,
        &cppad_objective,
        &cppad_gradient,
        &cppad_hessian
    );
    const int nvars = static_cast<int>(amounts.size()) + 1;
    if (cppad_gradient.size() != result.gradient.size()
        || cppad_hessian.size() != static_cast<std::size_t>(nvars * nvars)) {
        throw ValueError("EOS phase objective CppAD derivative shape did not match the phase variable model.");
    }
    const double objective_scale = std::max(1.0, std::abs(result.objective));
    if (std::abs(cppad_objective - result.objective) > 1.0e-8 * objective_scale) {
        throw ValueError("EOS phase objective CppAD value did not match the analytical block value.");
    }
    result.objective_curvature_backend = "cppad";
    result.objective_curvature_rows = nvars;
    result.objective_curvature_cols = nvars;
    result.objective_curvature_row_major = std::move(cppad_hessian);
    result.constraint_jacobian_backend = "cppad";
    result.constraint_jacobian_rows = 1;
    result.constraint_jacobian_cols = nvars;
    result.constraint_jacobian_row_major.reserve(static_cast<std::size_t>(nvars));
    const int volume_row = nvars - 1;
    for (int col = 0; col < nvars; ++col) {
        result.constraint_jacobian_row_major.push_back(
            -rt * result.objective_curvature_row_major[
                static_cast<std::size_t>(volume_row * nvars + col)
            ]
        );
    }
    result.pressure_jacobian_backend = result.constraint_jacobian_backend;
    result.pressure_jacobian_rows = result.constraint_jacobian_rows;
    result.pressure_jacobian_cols = result.constraint_jacobian_cols;
    result.pressure_jacobian_row_major = result.constraint_jacobian_row_major;
    result.pressure_density_derivative =
        -result.pressure_jacobian_row_major[static_cast<std::size_t>(volume_row)] * volume / density;
    return result;
}

EosPhaseSystemResult evaluate_eos_phase_system(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    const std::vector<double>& charges
) {
    if (phase_amounts.empty()) {
        throw ValueError("EOS phase system requires at least one phase.");
    }
    const std::size_t phase_count = phase_amounts.size();
    const std::size_t species_count = feed_amounts.size();
    if (species_count == 0) {
        throw ValueError("EOS phase system requires at least one feed species amount.");
    }
    if (volumes.size() != phase_count) {
        throw ValueError("EOS phase system volume count must match phase count.");
    }
    for (double amount : feed_amounts) {
        require_finite_nonnegative(amount, "EOS phase system feed amount");
    }

    EosPhaseSystemResult result;
    result.block = "eos_phase_system";
    result.derivative_backend = "analytic_cppad";
    result.phase_count = static_cast<int>(phase_count);
    result.species_count = static_cast<int>(species_count);
    result.variable_names = phase_system_variable_names(phase_count, species_count);
    result.constraint_names = phase_system_constraint_names(phase_count, species_count);
    result.temperature = temperature;
    result.target_pressure = target_pressure;
    result.feed_amounts = feed_amounts;
    result.phase_blocks.reserve(phase_count);

    const std::size_t local_variable_count = species_count + 1;
    result.gradient.reserve(phase_count * local_variable_count);
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        if (phase_amounts[phase].size() != species_count) {
            throw ValueError("EOS phase system phase amount sizes must match feed species count.");
        }
        result.phase_blocks.push_back(
            evaluate_eos_phase_block(args, temperature, target_pressure, phase_amounts[phase], volumes[phase])
        );
        const EosPhaseBlockResult& block = result.phase_blocks.back();
        if (block.gradient.size() != local_variable_count
            || block.pressure_jacobian_row_major.size() != local_variable_count) {
            throw ValueError("EOS phase block derivative size did not match the phase system variables.");
        }
        result.objective += block.objective;
        result.gradient.insert(result.gradient.end(), block.gradient.begin(), block.gradient.end());
    }

    if (!charges.empty() && charges.size() != species_count) {
        throw ValueError("EOS phase system charge count must match feed species count.");
    }
    const PhaseChargeBlockResult charge_block = charges.empty()
        ? PhaseChargeBlockResult{}
        : evaluate_phase_charge_block(phase_amounts, charges, static_cast<int>(local_variable_count));
    const std::size_t charge_constraint_count = charge_block.residuals.size();
    const std::size_t constraint_count = species_count + phase_count + charge_constraint_count;
    const std::size_t variable_count = phase_count * local_variable_count;
    result.constraints.assign(constraint_count, 0.0);
    for (std::size_t species = 0; species < species_count; ++species) {
        result.constraints[species] = -feed_amounts[species];
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            result.constraints[species] += phase_amounts[phase][species];
        }
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        result.constraints[species_count + phase] = result.phase_blocks[phase].pressure_consistency_residual;
    }
    for (std::size_t row = 0; row < charge_constraint_count; ++row) {
        result.constraint_names.push_back(charge_block.constraint_names[row]);
        result.constraints[species_count + phase_count + row] = charge_block.residuals[row];
    }
    result.phase_charge_residuals = charge_block.residuals;

    result.constraint_jacobian_backend = "analytic_cppad";
    result.constraint_jacobian_rows = static_cast<int>(constraint_count);
    result.constraint_jacobian_cols = static_cast<int>(variable_count);
    result.constraint_jacobian_row_major.assign(constraint_count * variable_count, 0.0);
    for (std::size_t species = 0; species < species_count; ++species) {
        for (std::size_t phase = 0; phase < phase_count; ++phase) {
            result.constraint_jacobian_row_major[
                species * variable_count + phase * local_variable_count + species
            ] = 1.0;
        }
    }
    for (std::size_t phase = 0; phase < phase_count; ++phase) {
        const EosPhaseBlockResult& block = result.phase_blocks[phase];
        if (block.pressure_jacobian_row_major.size() != local_variable_count) {
            throw ValueError("EOS phase pressure Jacobian size did not match the phase variable model.");
        }
        const std::size_t row = species_count + phase;
        const std::size_t col_offset = phase * local_variable_count;
        for (std::size_t col = 0; col < local_variable_count; ++col) {
            result.constraint_jacobian_row_major[row * variable_count + col_offset + col] =
                block.pressure_jacobian_row_major[col];
        }
    }
    for (std::size_t row = 0; row < charge_constraint_count; ++row) {
        const std::size_t target_row = species_count + phase_count + row;
        for (std::size_t col = 0; col < variable_count; ++col) {
            result.constraint_jacobian_row_major[target_row * variable_count + col] =
                charge_block.jacobian_row_major[row * variable_count + col];
        }
    }
    return result;
}

}  // namespace epcsaft::native::equilibrium_nlp
