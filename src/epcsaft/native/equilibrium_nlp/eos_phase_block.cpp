#include "eos_phase_block.h"

#include "epcsaft_core_internal.h"

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
    return result;
}

}  // namespace epcsaft::native::equilibrium_nlp
