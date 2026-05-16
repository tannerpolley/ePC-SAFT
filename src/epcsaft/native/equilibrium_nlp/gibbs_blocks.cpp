#include "gibbs_blocks.h"

#include "epcsaft_electrolyte.h"

#include <cmath>
#include <numeric>

namespace epcsaft::native::equilibrium_nlp {

namespace {

void validate_amounts_and_standard_state(
    const std::vector<double>& amounts,
    const std::vector<double>& standard_mu_rt
) {
    if (amounts.empty()) {
        throw ValueError("Ideal Gibbs block requires at least one species.");
    }
    if (amounts.size() != standard_mu_rt.size()) {
        throw ValueError("Ideal Gibbs block requires one standard chemical potential per species.");
    }
    for (std::size_t index = 0; index < amounts.size(); ++index) {
        if (!(std::isfinite(amounts[index]) && amounts[index] > 0.0)) {
            throw ValueError("Ideal Gibbs species amounts must be positive and finite.");
        }
        if (!std::isfinite(standard_mu_rt[index])) {
            throw ValueError("Ideal Gibbs standard chemical potentials must be finite.");
        }
    }
}

}  // namespace

IdealReducedGibbsResult evaluate_ideal_reduced_gibbs(
    const std::vector<double>& amounts,
    const std::vector<double>& standard_mu_rt,
    bool include_hessian
) {
    validate_amounts_and_standard_state(amounts, standard_mu_rt);
    const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    if (!(std::isfinite(total) && total > 0.0)) {
        throw ValueError("Ideal Gibbs total amount must be positive and finite.");
    }

    IdealReducedGibbsResult out;
    const std::size_t species = amounts.size();
    out.mole_fractions.resize(species);
    out.gradient.resize(species);
    for (std::size_t index = 0; index < species; ++index) {
        const double x_i = amounts[index] / total;
        out.mole_fractions[index] = x_i;
        out.value += amounts[index] * (std::log(x_i) + standard_mu_rt[index]);
        out.gradient[index] = std::log(x_i) + standard_mu_rt[index];
    }

    if (include_hessian) {
        out.hessian_row_major.assign(species * species, -1.0 / total);
        for (std::size_t index = 0; index < species; ++index) {
            out.hessian_row_major[index * species + index] += 1.0 / amounts[index];
        }
    }
    return out;
}

}  // namespace epcsaft::native::equilibrium_nlp
