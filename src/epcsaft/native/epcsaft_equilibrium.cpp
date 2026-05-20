#include "epcsaft_equilibrium.h"
#include "equilibrium/equilibrium_helpers.h"
#include "equilibrium_nlp/second_order.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cstddef>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

PhaseStateCompositionSensitivityResult phase_state_ln_fugacity_explicit_density_composition_sensitivity_cpp(
    double t,
    double rho,
    std::vector<double> x,
    const add_args &cppargs
);

namespace {

using namespace epcsaft::native::equilibrium;
namespace equilibrium_nlp = epcsaft::native::equilibrium_nlp;

constexpr double kLiquidLleMinimumPhaseDistance = 1.0e-1;
constexpr double kLiquidLleMinimumRetainedPhaseFraction = 1.0e-3;

struct PhaseStateNative {
    std::shared_ptr<ePCSAFTStateNative> state;
    std::vector<double> ln_phi;
    double density = 0.0;
};

struct ElectrolyteSaltPairNative {
    std::string label;
    int cation = -1;
    int anion = -1;
    int cation_stoich = 1;
    int anion_stoich = 1;
    double cation_charge = 1.0;
    double anion_charge = -1.0;
};

struct ElectrolyteBasisNative {
    std::vector<int> neutral_indices;
    std::vector<int> cation_indices;
    std::vector<int> anion_indices;
    std::vector<ElectrolyteSaltPairNative> salt_pairs;
    std::vector<double> formula_feed;
    std::vector<double> species_charges;
    int basis_rank = 0;
    std::string variable_model = "ascani_transformed_salt_pairs";
};

struct ElectrolyteCandidateNative {
    double beta_formula = 0.5;
    double beta_org = 0.5;
    std::vector<double> aq_formula;
    std::vector<double> org_formula;
    std::vector<double> aq_comp;
    std::vector<double> org_comp;
    PhaseStateNative aq_state;
    PhaseStateNative org_state;
    std::vector<double> residual;
    std::vector<double> material_residual;
    bool explicit_density = false;
    double solver_residual_norm = std::numeric_limits<double>::infinity();
    double material_error = std::numeric_limits<double>::infinity();
    double charge_balance_error = std::numeric_limits<double>::infinity();
    double gibbs_feed = 0.0;
    double gibbs_split = 0.0;
    double gibbs_delta = 0.0;
    double phase_distance_value = 0.0;
    int iteration = 0;
    double objective = std::numeric_limits<double>::infinity();
};

PhaseStateNative phase_state(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& composition,
    const std::string& phase,
    const std::string& label = ""
) {
    PhaseStateNative out;
    int phase_int = phase_token_to_int(phase);
    double rho = mixture->solve_density_scoped(t, p, composition, phase_int, label.empty() ? phase : label);
    out.state = mixture->state(t, composition, phase_int, false, 0.0, true, rho);
    out.ln_phi = out.state->ln_fugacity_coefficient();
    out.density = out.state->density();
    return out;
}

std::string ion_stem(const std::string& label, double charge = 0.0) {
    std::string out;
    for (char ch : label) {
        if (ch != '+' && ch != '-') {
            out.push_back(ch);
        }
    }
    int charge_int = static_cast<int>(std::round(std::abs(charge)));
    if (charge_int > 1 && (label.find('+') != std::string::npos || label.find('-') != std::string::npos)) {
        std::string suffix = std::to_string(charge_int);
        if (out.size() >= suffix.size() && out.substr(out.size() - suffix.size()) == suffix) {
            out.erase(out.size() - suffix.size());
        }
    }
    return out;
}

std::pair<int, int> neutral_salt_stoichiometry(double cation_charge, double anion_charge) {
    int cation_z = static_cast<int>(std::round(std::abs(cation_charge)));
    int anion_z = static_cast<int>(std::round(std::abs(anion_charge)));
    if (cation_z <= 0 || anion_z <= 0
        || std::abs(static_cast<double>(cation_z) - std::abs(cation_charge)) > 1.0e-12
        || std::abs(static_cast<double>(anion_z) - std::abs(anion_charge)) > 1.0e-12) {
        throw ValueError("electrolyte salt stoichiometry currently requires integer ion charges.");
    }
    int divisor = std::gcd(cation_z, anion_z);
    return {anion_z / divisor, cation_z / divisor};
}

ElectrolyteBasisNative build_electrolyte_basis_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<double>& feed,
    const std::vector<std::string>& species
) {
    const std::vector<double>& charges = mixture->args().z;
    if (charges.size() != feed.size()) {
        throw ValueError("mixture parameters must include one charge value per species in params['z'].");
    }
    ElectrolyteBasisNative basis;
    basis.species_charges = charges;
    for (std::size_t i = 0; i < charges.size(); ++i) {
        if (std::abs(charges[i]) <= 1.0e-12) {
            basis.neutral_indices.push_back(static_cast<int>(i));
        } else if (charges[i] > 1.0e-12) {
            basis.cation_indices.push_back(static_cast<int>(i));
        } else {
            basis.anion_indices.push_back(static_cast<int>(i));
        }
    }
    if (basis.neutral_indices.size() < 2 || basis.cation_indices.empty() || basis.anion_indices.empty()) {
        throw ValueError("electrolyte_lle requires at least two neutral species plus cations and anions.");
    }
    if (basis.anion_indices.size() == 1) {
        int anion = basis.anion_indices[0];
        for (int cation : basis.cation_indices) {
            auto stoich = neutral_salt_stoichiometry(charges[static_cast<std::size_t>(cation)], charges[static_cast<std::size_t>(anion)]);
            ElectrolyteSaltPairNative pair;
            std::string cat_label = species.size() == feed.size() ? species[static_cast<std::size_t>(cation)] : std::to_string(cation);
            std::string an_label = species.size() == feed.size() ? species[static_cast<std::size_t>(anion)] : std::to_string(anion);
            pair.cation = cation;
            pair.anion = anion;
            pair.cation_stoich = stoich.first;
            pair.anion_stoich = stoich.second;
            pair.cation_charge = charges[static_cast<std::size_t>(cation)];
            pair.anion_charge = charges[static_cast<std::size_t>(anion)];
            pair.label = ion_stem(cat_label, pair.cation_charge)
                + (pair.cation_stoich == 1 ? "" : std::to_string(pair.cation_stoich))
                + ion_stem(an_label, pair.anion_charge)
                + (pair.anion_stoich == 1 ? "" : std::to_string(pair.anion_stoich));
            basis.salt_pairs.push_back(pair);
        }
    } else if (basis.cation_indices.size() == 1) {
        int cation = basis.cation_indices[0];
        for (int anion : basis.anion_indices) {
            auto stoich = neutral_salt_stoichiometry(charges[static_cast<std::size_t>(cation)], charges[static_cast<std::size_t>(anion)]);
            ElectrolyteSaltPairNative pair;
            std::string cat_label = species.size() == feed.size() ? species[static_cast<std::size_t>(cation)] : std::to_string(cation);
            std::string an_label = species.size() == feed.size() ? species[static_cast<std::size_t>(anion)] : std::to_string(anion);
            pair.cation = cation;
            pair.anion = anion;
            pair.cation_stoich = stoich.first;
            pair.anion_stoich = stoich.second;
            pair.cation_charge = charges[static_cast<std::size_t>(cation)];
            pair.anion_charge = charges[static_cast<std::size_t>(anion)];
            pair.label = ion_stem(cat_label, pair.cation_charge)
                + (pair.cation_stoich == 1 ? "" : std::to_string(pair.cation_stoich))
                + ion_stem(an_label, pair.anion_charge)
                + (pair.anion_stoich == 1 ? "" : std::to_string(pair.anion_stoich));
            basis.salt_pairs.push_back(pair);
        }
    } else {
        int anchor_anion = basis.anion_indices[0];
        for (int cation : basis.cation_indices) {
            auto stoich = neutral_salt_stoichiometry(charges[static_cast<std::size_t>(cation)], charges[static_cast<std::size_t>(anchor_anion)]);
            ElectrolyteSaltPairNative pair;
            std::string cat_label = species.size() == feed.size() ? species[static_cast<std::size_t>(cation)] : std::to_string(cation);
            std::string an_label = species.size() == feed.size() ? species[static_cast<std::size_t>(anchor_anion)] : std::to_string(anchor_anion);
            pair.cation = cation;
            pair.anion = anchor_anion;
            pair.cation_stoich = stoich.first;
            pair.anion_stoich = stoich.second;
            pair.cation_charge = charges[static_cast<std::size_t>(cation)];
            pair.anion_charge = charges[static_cast<std::size_t>(anchor_anion)];
            pair.label = ion_stem(cat_label, pair.cation_charge)
                + (pair.cation_stoich == 1 ? "" : std::to_string(pair.cation_stoich))
                + ion_stem(an_label, pair.anion_charge)
                + (pair.anion_stoich == 1 ? "" : std::to_string(pair.anion_stoich));
            basis.salt_pairs.push_back(pair);
        }
        int anchor_cation = basis.cation_indices[0];
        for (std::size_t pos = 1; pos < basis.anion_indices.size(); ++pos) {
            int anion = basis.anion_indices[pos];
            auto stoich = neutral_salt_stoichiometry(charges[static_cast<std::size_t>(anchor_cation)], charges[static_cast<std::size_t>(anion)]);
            ElectrolyteSaltPairNative pair;
            std::string cat_label = species.size() == feed.size() ? species[static_cast<std::size_t>(anchor_cation)] : std::to_string(anchor_cation);
            std::string an_label = species.size() == feed.size() ? species[static_cast<std::size_t>(anion)] : std::to_string(anion);
            pair.cation = anchor_cation;
            pair.anion = anion;
            pair.cation_stoich = stoich.first;
            pair.anion_stoich = stoich.second;
            pair.cation_charge = charges[static_cast<std::size_t>(anchor_cation)];
            pair.anion_charge = charges[static_cast<std::size_t>(anion)];
            pair.label = ion_stem(cat_label, pair.cation_charge)
                + (pair.cation_stoich == 1 ? "" : std::to_string(pair.cation_stoich))
                + ion_stem(an_label, pair.anion_charge)
                + (pair.anion_stoich == 1 ? "" : std::to_string(pair.anion_stoich));
            basis.salt_pairs.push_back(pair);
        }
    }
    std::vector<double> reconstructed_charged(feed.size(), 0.0);
    for (const auto& pair : basis.salt_pairs) {
        double amount = feed[static_cast<std::size_t>(pair.cation)] / static_cast<double>(pair.cation_stoich);
        reconstructed_charged[static_cast<std::size_t>(pair.cation)] += amount * static_cast<double>(pair.cation_stoich);
        reconstructed_charged[static_cast<std::size_t>(pair.anion)] += amount * static_cast<double>(pair.anion_stoich);
    }
    for (int index : basis.cation_indices) {
        if (std::abs(reconstructed_charged[static_cast<std::size_t>(index)] - feed[static_cast<std::size_t>(index)]) > 1.0e-8) {
            throw ValueError("electrolyte_lle feed cannot be represented on the charge-constrained salt basis.");
        }
    }
    for (int index : basis.anion_indices) {
        if (std::abs(reconstructed_charged[static_cast<std::size_t>(index)] - feed[static_cast<std::size_t>(index)]) > 1.0e-8) {
            throw ValueError("electrolyte_lle feed cannot be represented on the charge-constrained salt basis.");
        }
    }
    std::vector<double> formula_moles;
    for (int index : basis.neutral_indices) {
        formula_moles.push_back(feed[static_cast<std::size_t>(index)]);
    }
    for (const auto& pair : basis.salt_pairs) {
        formula_moles.push_back(feed[static_cast<std::size_t>(pair.cation)] / static_cast<double>(pair.cation_stoich));
    }
    basis.formula_feed = clip_normalize(formula_moles, 1.0e-300);
    basis.basis_rank = static_cast<int>(basis.salt_pairs.size());
    return basis;
}

std::pair<std::vector<double>, double> formula_to_explicit(
    const std::vector<double>& formula,
    const ElectrolyteBasisNative& basis,
    std::size_t ncomp
) {
    std::vector<double> explicit_x(ncomp, 0.0);
    for (std::size_t pos = 0; pos < basis.neutral_indices.size(); ++pos) {
        explicit_x[static_cast<std::size_t>(basis.neutral_indices[pos])] += formula[pos];
    }
    std::size_t offset = basis.neutral_indices.size();
    for (std::size_t salt_pos = 0; salt_pos < basis.salt_pairs.size(); ++salt_pos) {
        double amount = formula[offset + salt_pos];
        const auto& pair = basis.salt_pairs[salt_pos];
        explicit_x[static_cast<std::size_t>(pair.cation)] += amount * static_cast<double>(pair.cation_stoich);
        explicit_x[static_cast<std::size_t>(pair.anion)] += amount * static_cast<double>(pair.anion_stoich);
    }
    double total = std::accumulate(explicit_x.begin(), explicit_x.end(), 0.0);
    if (total <= 0.0 || !std::isfinite(total)) {
        throw SolutionError("Formula-basis electrolyte phase expanded to a non-positive explicit composition.");
    }
    for (double& value : explicit_x) {
        value /= total;
    }
    return {explicit_x, total};
}

std::vector<double> explicit_to_formula(const std::vector<double>& composition, const ElectrolyteBasisNative& basis) {
    std::vector<double> values;
    for (int index : basis.neutral_indices) {
        values.push_back(composition[static_cast<std::size_t>(index)]);
    }
    for (const auto& pair : basis.salt_pairs) {
        values.push_back(composition[static_cast<std::size_t>(pair.cation)] / static_cast<double>(pair.cation_stoich));
    }
    return clip_normalize(values, 1.0e-300);
}

std::vector<double> electrolyte_basis_vectors_row_major(const ElectrolyteBasisNative& basis, std::size_t ncomp) {
    std::vector<double> out(basis.salt_pairs.size() * ncomp, 0.0);
    for (std::size_t row = 0; row < basis.salt_pairs.size(); ++row) {
        const auto& pair = basis.salt_pairs[row];
        out[row * ncomp + static_cast<std::size_t>(pair.cation)] = static_cast<double>(pair.cation_stoich);
        out[row * ncomp + static_cast<std::size_t>(pair.anion)] = static_cast<double>(pair.anion_stoich);
    }
    return out;
}

double electrolyte_gibbs_proxy(const std::vector<double>& composition, const PhaseStateNative& state) {
    double out = 0.0;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        out += composition[i] * (std::log(composition[i]) + state.ln_phi[i]);
    }
    return out;
}

std::vector<double> electrolyte_residual_vector(
    const std::vector<double>& aq_comp,
    const std::vector<double>& org_comp,
    const PhaseStateNative& aq_state,
    const PhaseStateNative& org_state,
    const ElectrolyteBasisNative& basis
) {
    std::vector<double> residuals;
    for (int index : basis.neutral_indices) {
        std::size_t i = static_cast<std::size_t>(index);
        residuals.push_back((std::log(org_comp[i]) + org_state.ln_phi[i]) - (std::log(aq_comp[i]) + aq_state.ln_phi[i]));
    }
    for (const auto& pair : basis.salt_pairs) {
        std::size_t c = static_cast<std::size_t>(pair.cation);
        std::size_t a = static_cast<std::size_t>(pair.anion);
        double org_pair = pair.cation_stoich * (std::log(org_comp[c]) + org_state.ln_phi[c])
            + pair.anion_stoich * (std::log(org_comp[a]) + org_state.ln_phi[a]);
        double aq_pair = pair.cation_stoich * (std::log(aq_comp[c]) + aq_state.ln_phi[c])
            + pair.anion_stoich * (std::log(aq_comp[a]) + aq_state.ln_phi[a]);
        residuals.push_back(org_pair - aq_pair);
    }
    return residuals;
}

struct ElectrolyteTransformDerivativesNative {
    std::vector<double> aq_formula_dvar;
    std::vector<double> org_formula_dvar;
    std::vector<double> aq_formula_hessian_row_major;
    std::vector<double> org_formula_hessian_row_major;
    std::vector<double> aq_comp_dvar;
    std::vector<double> org_comp_dvar;
    std::vector<double> aq_comp_hessian_row_major;
    std::vector<double> org_comp_hessian_row_major;
    std::vector<double> beta_org_dvar;
    std::vector<double> beta_org_hessian_row_major;
    int ncomp = 0;
    int nformula = 0;
    int nvar = 0;
};

struct ScalarRouteDerivativeState {
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian_row_major;
};

ScalarRouteDerivativeState make_constant_scalar_state(double value, std::size_t nvar) {
    ScalarRouteDerivativeState out;
    out.value = value;
    out.gradient.assign(nvar, 0.0);
    out.hessian_row_major.assign(nvar * nvar, 0.0);
    return out;
}

PhaseStateNative phase_state_at_density(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double rho,
    const std::vector<double>& composition,
    const std::string& phase
) {
    if (!std::isfinite(rho) || !(rho > 0.0)) {
        throw SolutionError("Explicit-density phase state requires a positive finite density.");
    }
    PhaseStateNative out;
    const int phase_int = phase_token_to_int(phase);
    out.state = mixture->state(t, composition, phase_int, false, 0.0, true, rho);
    out.ln_phi = out.state->ln_fugacity_coefficient();
    out.density = out.state->density();
    return out;
}

ScalarRouteDerivativeState add_scalar_states(
    const ScalarRouteDerivativeState& left,
    const ScalarRouteDerivativeState& right
) {
    ScalarRouteDerivativeState out = left;
    for (std::size_t i = 0; i < out.gradient.size(); ++i) {
        out.gradient[i] += right.gradient[i];
    }
    for (std::size_t i = 0; i < out.hessian_row_major.size(); ++i) {
        out.hessian_row_major[i] += right.hessian_row_major[i];
    }
    out.value += right.value;
    return out;
}

ScalarRouteDerivativeState subtract_scalar_states(
    const ScalarRouteDerivativeState& left,
    const ScalarRouteDerivativeState& right
) {
    ScalarRouteDerivativeState out = left;
    for (std::size_t i = 0; i < out.gradient.size(); ++i) {
        out.gradient[i] -= right.gradient[i];
    }
    for (std::size_t i = 0; i < out.hessian_row_major.size(); ++i) {
        out.hessian_row_major[i] -= right.hessian_row_major[i];
    }
    out.value -= right.value;
    return out;
}

ScalarRouteDerivativeState scale_scalar_state(
    const ScalarRouteDerivativeState& state,
    double factor
) {
    ScalarRouteDerivativeState out = state;
    out.value *= factor;
    for (double& value : out.gradient) {
        value *= factor;
    }
    for (double& value : out.hessian_row_major) {
        value *= factor;
    }
    return out;
}

ScalarRouteDerivativeState multiply_scalar_states(
    const ScalarRouteDerivativeState& left,
    const ScalarRouteDerivativeState& right
) {
    const std::size_t nvar = left.gradient.size();
    ScalarRouteDerivativeState out;
    out.value = left.value * right.value;
    out.gradient.assign(nvar, 0.0);
    out.hessian_row_major.assign(nvar * nvar, 0.0);
    for (std::size_t i = 0; i < nvar; ++i) {
        out.gradient[i] = left.gradient[i] * right.value + left.value * right.gradient[i];
        for (std::size_t j = 0; j < nvar; ++j) {
            out.hessian_row_major[i * nvar + j] =
                left.hessian_row_major[i * nvar + j] * right.value
                + left.gradient[i] * right.gradient[j]
                + left.gradient[j] * right.gradient[i]
                + left.value * right.hessian_row_major[i * nvar + j];
        }
    }
    return out;
}

ScalarRouteDerivativeState divide_scalar_states(
    const ScalarRouteDerivativeState& numerator,
    const ScalarRouteDerivativeState& denominator
) {
    const std::size_t nvar = numerator.gradient.size();
    if (!(std::isfinite(denominator.value) && std::abs(denominator.value) > 0.0)) {
        throw SolutionError("Electrolyte transformed-variable derivative received an invalid denominator.");
    }
    ScalarRouteDerivativeState out;
    out.value = numerator.value / denominator.value;
    out.gradient.assign(nvar, 0.0);
    out.hessian_row_major.assign(nvar * nvar, 0.0);
    const double inv_den = 1.0 / denominator.value;
    const double inv_den2 = inv_den * inv_den;
    const double inv_den3 = inv_den2 * inv_den;
    for (std::size_t i = 0; i < nvar; ++i) {
        out.gradient[i] =
            numerator.gradient[i] * inv_den
            - numerator.value * denominator.gradient[i] * inv_den2;
        for (std::size_t j = 0; j < nvar; ++j) {
            out.hessian_row_major[i * nvar + j] =
                numerator.hessian_row_major[i * nvar + j] * inv_den
                - (
                    numerator.gradient[i] * denominator.gradient[j]
                    + numerator.gradient[j] * denominator.gradient[i]
                    + numerator.value * denominator.hessian_row_major[i * nvar + j]
                ) * inv_den2
                + 2.0 * numerator.value
                    * denominator.gradient[i]
                    * denominator.gradient[j]
                    * inv_den3;
        }
    }
    return out;
}

std::vector<double> formula_expansion_matrix_row_major(const ElectrolyteBasisNative& basis, std::size_t ncomp) {
    const std::size_t nformula = basis.formula_feed.size();
    std::vector<double> matrix(ncomp * nformula, 0.0);
    for (std::size_t pos = 0; pos < basis.neutral_indices.size(); ++pos) {
        matrix[static_cast<std::size_t>(basis.neutral_indices[pos]) * nformula + pos] = 1.0;
    }
    std::size_t offset = basis.neutral_indices.size();
    for (std::size_t salt_pos = 0; salt_pos < basis.salt_pairs.size(); ++salt_pos) {
        const auto& pair = basis.salt_pairs[salt_pos];
        const std::size_t col = offset + salt_pos;
        matrix[static_cast<std::size_t>(pair.cation) * nformula + col] = static_cast<double>(pair.cation_stoich);
        matrix[static_cast<std::size_t>(pair.anion) * nformula + col] = static_cast<double>(pair.anion_stoich);
    }
    return matrix;
}

std::vector<double> formula_expansion_column_sums(const ElectrolyteBasisNative& basis, std::size_t ncomp) {
    const std::size_t nformula = basis.formula_feed.size();
    std::vector<double> sums(nformula, 0.0);
    std::vector<double> matrix = formula_expansion_matrix_row_major(basis, ncomp);
    for (std::size_t col = 0; col < nformula; ++col) {
        for (std::size_t row = 0; row < ncomp; ++row) {
            sums[col] += matrix[row * nformula + col];
        }
    }
    return sums;
}

std::vector<double> formula_to_explicit_jacobian_row_major(
    const std::vector<double>& formula,
    const ElectrolyteBasisNative& basis,
    std::size_t ncomp
) {
    const std::size_t nformula = formula.size();
    std::vector<double> matrix = formula_expansion_matrix_row_major(basis, ncomp);
    std::vector<double> amounts(ncomp, 0.0);
    std::vector<double> col_sums(nformula, 0.0);
    for (std::size_t i = 0; i < ncomp; ++i) {
        for (std::size_t k = 0; k < nformula; ++k) {
            double value = matrix[i * nformula + k];
            amounts[i] += value * formula[k];
            col_sums[k] += value;
        }
    }
    double scale = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    if (!(scale > 0.0) || !std::isfinite(scale)) {
        throw SolutionError("Formula-basis explicit-composition Jacobian received a non-positive expansion scale.");
    }
    std::vector<double> jacobian(ncomp * nformula, 0.0);
    for (std::size_t i = 0; i < ncomp; ++i) {
        for (std::size_t k = 0; k < nformula; ++k) {
            jacobian[i * nformula + k] = (matrix[i * nformula + k] * scale - amounts[i] * col_sums[k])
                / (scale * scale);
        }
    }
    return jacobian;
}

ElectrolyteTransformDerivativesNative electrolyte_transform_derivatives(
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const ElectrolyteCandidateNative& candidate,
    std::size_t total_variable_count = 0
) {
    const std::size_t ncomp = feed.size();
    const std::size_t nformula = basis.formula_feed.size();
    const std::size_t nvar = total_variable_count == 0 ? nformula : total_variable_count;
    ElectrolyteTransformDerivativesNative out;
    out.ncomp = static_cast<int>(ncomp);
    out.nformula = static_cast<int>(nformula);
    out.nvar = static_cast<int>(nvar);
    out.aq_formula_dvar.assign(nformula * nvar, 0.0);
    out.org_formula_dvar.assign(nformula * nvar, 0.0);
    out.aq_formula_hessian_row_major.assign(nformula * nvar * nvar, 0.0);
    out.org_formula_hessian_row_major.assign(nformula * nvar * nvar, 0.0);
    out.aq_comp_dvar.assign(ncomp * nvar, 0.0);
    out.org_comp_dvar.assign(ncomp * nvar, 0.0);
    out.aq_comp_hessian_row_major.assign(ncomp * nvar * nvar, 0.0);
    out.org_comp_hessian_row_major.assign(ncomp * nvar * nvar, 0.0);
    out.beta_org_dvar.assign(nvar, 0.0);
    std::vector<double> column_sums = formula_expansion_column_sums(basis, ncomp);
    std::vector<double> matrix = formula_expansion_matrix_row_major(basis, ncomp);
    out.beta_org_hessian_row_major.assign(nvar * nvar, 0.0);

    ScalarRouteDerivativeState beta_formula = make_constant_scalar_state(candidate.beta_formula, nvar);
    beta_formula.gradient[0] = candidate.beta_formula * (1.0 - candidate.beta_formula);
    beta_formula.hessian_row_major[0] = beta_formula.gradient[0] * (1.0 - 2.0 * candidate.beta_formula);

    std::vector<ScalarRouteDerivativeState> org_formula_states;
    org_formula_states.reserve(nformula);
    for (std::size_t component = 0; component < nformula; ++component) {
        ScalarRouteDerivativeState state = make_constant_scalar_state(candidate.org_formula[component], nvar);
        for (std::size_t var = 1; var < nformula; ++var) {
            const std::size_t alpha = var - 1;
            const double delta_component_alpha = component == alpha ? 1.0 : 0.0;
            state.gradient[var] =
                candidate.org_formula[component] * (delta_component_alpha - candidate.org_formula[alpha]);
        }
        for (std::size_t var_i = 1; var_i < nformula; ++var_i) {
            const std::size_t alpha = var_i - 1;
            const double delta_component_alpha = component == alpha ? 1.0 : 0.0;
            for (std::size_t var_j = 1; var_j < nformula; ++var_j) {
                const std::size_t beta = var_j - 1;
                const double delta_component_beta = component == beta ? 1.0 : 0.0;
                const double delta_alpha_beta = alpha == beta ? 1.0 : 0.0;
                state.hessian_row_major[var_i * nvar + var_j] =
                    candidate.org_formula[component] * (
                        (delta_component_beta - candidate.org_formula[beta])
                        * (delta_component_alpha - candidate.org_formula[alpha])
                        - candidate.org_formula[alpha] * (delta_alpha_beta - candidate.org_formula[beta])
                    );
            }
        }
        org_formula_states.push_back(std::move(state));
    }

    ScalarRouteDerivativeState one_minus_beta =
        subtract_scalar_states(make_constant_scalar_state(1.0, nvar), beta_formula);
    std::vector<ScalarRouteDerivativeState> aq_formula_states;
    aq_formula_states.reserve(nformula);
    for (std::size_t component = 0; component < nformula; ++component) {
        ScalarRouteDerivativeState numerator = subtract_scalar_states(
            make_constant_scalar_state(basis.formula_feed[component], nvar),
            multiply_scalar_states(beta_formula, org_formula_states[component])
        );
        aq_formula_states.push_back(divide_scalar_states(numerator, one_minus_beta));
    }

    auto linear_combination = [&](const std::vector<ScalarRouteDerivativeState>& states, const std::vector<double>& weights) {
        ScalarRouteDerivativeState out_state = make_constant_scalar_state(0.0, nvar);
        for (std::size_t k = 0; k < states.size(); ++k) {
            out_state = add_scalar_states(out_state, scale_scalar_state(states[k], weights[k]));
        }
        return out_state;
    };

    auto fill_formula_outputs = [&](const std::vector<ScalarRouteDerivativeState>& states,
                                    std::vector<double>& grad_out,
                                    std::vector<double>& hess_out) {
        for (std::size_t component = 0; component < states.size(); ++component) {
            for (std::size_t var = 0; var < nvar; ++var) {
                grad_out[component * nvar + var] = states[component].gradient[var];
                for (std::size_t other = 0; other < nvar; ++other) {
                    hess_out[component * nvar * nvar + var * nvar + other] =
                        states[component].hessian_row_major[var * nvar + other];
                }
            }
        }
    };

    fill_formula_outputs(org_formula_states, out.org_formula_dvar, out.org_formula_hessian_row_major);
    fill_formula_outputs(aq_formula_states, out.aq_formula_dvar, out.aq_formula_hessian_row_major);

    auto explicit_states_from_formula = [&](const std::vector<ScalarRouteDerivativeState>& formula_states) {
        std::vector<ScalarRouteDerivativeState> amount_states;
        amount_states.reserve(ncomp);
        for (std::size_t species = 0; species < ncomp; ++species) {
            ScalarRouteDerivativeState amount = make_constant_scalar_state(0.0, nvar);
            for (std::size_t component = 0; component < nformula; ++component) {
                amount = add_scalar_states(
                    amount,
                    scale_scalar_state(
                        formula_states[component],
                        matrix[species * nformula + component]
                    )
                );
            }
            amount_states.push_back(std::move(amount));
        }
        ScalarRouteDerivativeState scale = linear_combination(formula_states, column_sums);
        if (!(scale.value > 0.0) || !std::isfinite(scale.value)) {
            throw SolutionError("Electrolyte transformed-variable derivative produced an invalid expansion scale.");
        }
        std::vector<ScalarRouteDerivativeState> composition_states;
        composition_states.reserve(ncomp);
        for (std::size_t species = 0; species < ncomp; ++species) {
            composition_states.push_back(divide_scalar_states(amount_states[species], scale));
        }
        return std::make_pair(std::move(composition_states), std::move(scale));
    };

    auto [aq_comp_states, aq_scale] = explicit_states_from_formula(aq_formula_states);
    auto [org_comp_states, org_scale] = explicit_states_from_formula(org_formula_states);

    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t var = 0; var < nvar; ++var) {
            out.aq_comp_dvar[species * nvar + var] = aq_comp_states[species].gradient[var];
            out.org_comp_dvar[species * nvar + var] = org_comp_states[species].gradient[var];
            for (std::size_t other = 0; other < nvar; ++other) {
                out.aq_comp_hessian_row_major[species * nvar * nvar + var * nvar + other] =
                    aq_comp_states[species].hessian_row_major[var * nvar + other];
                out.org_comp_hessian_row_major[species * nvar * nvar + var * nvar + other] =
                    org_comp_states[species].hessian_row_major[var * nvar + other];
            }
        }
    }

    ScalarRouteDerivativeState beta_numerator = multiply_scalar_states(beta_formula, org_scale);
    ScalarRouteDerivativeState beta_denominator = add_scalar_states(
        multiply_scalar_states(one_minus_beta, aq_scale),
        multiply_scalar_states(beta_formula, org_scale)
    );
    ScalarRouteDerivativeState beta_org = divide_scalar_states(beta_numerator, beta_denominator);
    out.beta_org_dvar = beta_org.gradient;
    out.beta_org_hessian_row_major = beta_org.hessian_row_major;
    return out;
}

double phase_log_fugacity_derivative_for_species(
    const std::vector<double>& composition,
    const std::vector<double>& dcomp_dvar,
    const std::vector<double>& lnphi_jacobian,
    const std::vector<double>& lnphi_density_derivative,
    double density,
    std::size_t density_var,
    bool has_density_var,
    std::size_t species_index,
    std::size_t var,
    std::size_t ncomp,
    std::size_t nvar
) {
    double value = dcomp_dvar[species_index * nvar + var] / composition[species_index];
    for (std::size_t j = 0; j < ncomp; ++j) {
        value += lnphi_jacobian[species_index * ncomp + j] * dcomp_dvar[j * nvar + var];
    }
    if (has_density_var && var == density_var) {
        value += density * lnphi_density_derivative[species_index];
    }
    return value;
}

double phase_log_fugacity_second_derivative_for_species(
    const std::vector<double>& composition,
    const std::vector<double>& dcomp_dvar,
    const std::vector<double>& d2comp_dvar2,
    const std::vector<double>& lnphi_jacobian,
    const std::vector<double>& lnphi_hessian_tensor,
    const std::vector<double>& lnphi_density_derivative,
    const std::vector<double>& lnphi_density_second_derivative,
    const std::vector<double>& lnphi_density_composition_cross_derivative,
    double density,
    std::size_t density_var,
    bool has_density_var,
    std::size_t species_index,
    std::size_t first_var,
    std::size_t second_var,
    std::size_t ncomp,
    std::size_t nvar
) {
    const double comp = composition[species_index];
    const double dcomp_first = dcomp_dvar[species_index * nvar + first_var];
    const double dcomp_second = dcomp_dvar[species_index * nvar + second_var];
    double value =
        d2comp_dvar2[species_index * nvar * nvar + first_var * nvar + second_var] / comp
        - dcomp_first * dcomp_second / (comp * comp);
    for (std::size_t j = 0; j < ncomp; ++j) {
        value += lnphi_jacobian[species_index * ncomp + j]
            * d2comp_dvar2[j * nvar * nvar + first_var * nvar + second_var];
        for (std::size_t k = 0; k < ncomp; ++k) {
            value += lnphi_hessian_tensor[species_index * ncomp * ncomp + j * ncomp + k]
                * dcomp_dvar[j * nvar + first_var]
                * dcomp_dvar[k * nvar + second_var];
        }
        if (has_density_var) {
            if (first_var == density_var) {
                value += density
                    * lnphi_density_composition_cross_derivative[species_index * ncomp + j]
                    * dcomp_dvar[j * nvar + second_var];
            }
            if (second_var == density_var) {
                value += density
                    * lnphi_density_composition_cross_derivative[species_index * ncomp + j]
                    * dcomp_dvar[j * nvar + first_var];
            }
        }
    }
    if (has_density_var && first_var == density_var && second_var == density_var) {
        value += density * density * lnphi_density_second_derivative[species_index]
            + density * lnphi_density_derivative[species_index];
    }
    return value;
}

struct ElectrolyteResidualSecondOrderNative {
    std::vector<double> jacobian_row_major;
    std::vector<double> residual_hessian_tensor_row_major;
    std::vector<double> phase_separation_hessian_row_major;
    std::vector<double> aqueous_formula_hessian_tensor_row_major;
    std::size_t residual_rows = 0;
    std::size_t nvar = 0;
};

PhaseStateCompositionSensitivityResult electrolyte_phase_sensitivity(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& composition,
    double density,
    bool explicit_density
) {
    (void)p;
    if (!explicit_density) {
        throw SolutionError("Electrolyte LLE residual derivatives require explicit density variables.");
    }
    PhaseStateCompositionSensitivityResult out =
        phase_state_ln_fugacity_explicit_density_composition_sensitivity_cpp(
            t,
            density,
            composition,
            mixture->args()
        );
    if (!out.supported) {
        throw SolutionError("Electrolyte residual derivatives require supported phase-state fugacity sensitivities.");
    }
    return out;
}

std::size_t electrolyte_density_variable_index(
    std::size_t formula_variable_count,
    bool organic_phase
) {
    return formula_variable_count + (organic_phase ? 1 : 0);
}

double electrolyte_pressure_constraint_scale(
    double target_pressure,
    double density,
    const PhaseStateCompositionSensitivityResult& sensitivity
) {
    return std::max({std::abs(target_pressure), std::abs(density * sensitivity.pressure_density_derivative), 1.0});
}

std::vector<double> electrolyte_residual_jacobian_row_major(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const ElectrolyteCandidateNative& candidate
) {
    const std::size_t ncomp = feed.size();
    const std::size_t formula_nvar = basis.formula_feed.size();
    const std::size_t nvar = formula_nvar + (candidate.explicit_density ? 2 : 0);
    ElectrolyteTransformDerivativesNative transform = electrolyte_transform_derivatives(feed, basis, candidate, nvar);
    PhaseStateCompositionSensitivityResult aq_sensitivity = electrolyte_phase_sensitivity(
        mixture,
        t,
        p,
        candidate.aq_comp,
        candidate.aq_state.density,
        candidate.explicit_density
    );
    PhaseStateCompositionSensitivityResult org_sensitivity = electrolyte_phase_sensitivity(
        mixture,
        t,
        p,
        candidate.org_comp,
        candidate.org_state.density,
        candidate.explicit_density
    );
    const std::size_t phase_rows = basis.neutral_indices.size() + basis.salt_pairs.size();
    const std::size_t rows = phase_rows + ncomp;
    std::vector<double> jacobian(rows * nvar, 0.0);
    std::size_t row = 0;
    auto add_species_contribution = [&](std::size_t output_row, std::size_t species_index, double coefficient) {
        for (std::size_t var = 0; var < nvar; ++var) {
            double org_value = phase_log_fugacity_derivative_for_species(
                candidate.org_comp,
                transform.org_comp_dvar,
                org_sensitivity.fixed_density_jacobian_row_major,
                org_sensitivity.ln_fugacity_density_derivative,
                candidate.org_state.density,
                electrolyte_density_variable_index(formula_nvar, true),
                candidate.explicit_density,
                species_index,
                var,
                ncomp,
                nvar
            );
            double aq_value = phase_log_fugacity_derivative_for_species(
                candidate.aq_comp,
                transform.aq_comp_dvar,
                aq_sensitivity.fixed_density_jacobian_row_major,
                aq_sensitivity.ln_fugacity_density_derivative,
                candidate.aq_state.density,
                electrolyte_density_variable_index(formula_nvar, false),
                candidate.explicit_density,
                species_index,
                var,
                ncomp,
                nvar
            );
            jacobian[output_row * nvar + var] += coefficient * (org_value - aq_value);
        }
    };
    for (int index : basis.neutral_indices) {
        add_species_contribution(row, static_cast<std::size_t>(index), 1.0);
        ++row;
    }
    for (const auto& pair : basis.salt_pairs) {
        add_species_contribution(row, static_cast<std::size_t>(pair.cation), static_cast<double>(pair.cation_stoich));
        add_species_contribution(row, static_cast<std::size_t>(pair.anion), static_cast<double>(pair.anion_stoich));
        ++row;
    }
    for (std::size_t i = 0; i < ncomp; ++i) {
        for (std::size_t var = 0; var < nvar; ++var) {
            jacobian[row * nvar + var] =
                (candidate.org_comp[i] - candidate.aq_comp[i]) * transform.beta_org_dvar[var]
                + (1.0 - candidate.beta_org) * transform.aq_comp_dvar[i * nvar + var]
                + candidate.beta_org * transform.org_comp_dvar[i * nvar + var];
        }
        ++row;
    }
    return jacobian;
}

ElectrolyteResidualSecondOrderNative electrolyte_residual_second_order_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const ElectrolyteCandidateNative& candidate,
    int separation_formula_index,
    double separation_sign
) {
    const std::size_t ncomp = feed.size();
    const std::size_t formula_nvar = basis.formula_feed.size();
    const std::size_t nvar = formula_nvar + (candidate.explicit_density ? 2 : 0);
    ElectrolyteTransformDerivativesNative transform = electrolyte_transform_derivatives(feed, basis, candidate, nvar);
    PhaseStateCompositionSensitivityResult aq_sensitivity = electrolyte_phase_sensitivity(
        mixture,
        t,
        p,
        candidate.aq_comp,
        candidate.aq_state.density,
        candidate.explicit_density
    );
    PhaseStateCompositionSensitivityResult org_sensitivity = electrolyte_phase_sensitivity(
        mixture,
        t,
        p,
        candidate.org_comp,
        candidate.org_state.density,
        candidate.explicit_density
    );
    const std::size_t phase_rows = basis.neutral_indices.size() + basis.salt_pairs.size();
    const std::size_t residual_rows = phase_rows + ncomp;
    ElectrolyteResidualSecondOrderNative out;
    out.residual_rows = residual_rows;
    out.nvar = nvar;
    out.jacobian_row_major = electrolyte_residual_jacobian_row_major(mixture, t, p, feed, basis, candidate);
    out.residual_hessian_tensor_row_major.assign(residual_rows * nvar * nvar, 0.0);
    out.phase_separation_hessian_row_major.assign(nvar * nvar, 0.0);
    out.aqueous_formula_hessian_tensor_row_major.assign(basis.formula_feed.size() * nvar * nvar, 0.0);

    auto add_species_contribution = [&](std::size_t output_row, std::size_t species_index, double coefficient) {
        for (std::size_t first = 0; first < nvar; ++first) {
            for (std::size_t second = 0; second < nvar; ++second) {
                const double org_value = phase_log_fugacity_second_derivative_for_species(
                    candidate.org_comp,
                    transform.org_comp_dvar,
                    transform.org_comp_hessian_row_major,
                    org_sensitivity.fixed_density_jacobian_row_major,
                    org_sensitivity.fixed_density_hessian_tensor_row_major,
                    org_sensitivity.ln_fugacity_density_derivative,
                    org_sensitivity.ln_fugacity_density_second_derivative,
                    org_sensitivity.ln_fugacity_density_composition_cross_derivative,
                    candidate.org_state.density,
                    electrolyte_density_variable_index(formula_nvar, true),
                    candidate.explicit_density,
                    species_index,
                    first,
                    second,
                    ncomp,
                    nvar
                );
                const double aq_value = phase_log_fugacity_second_derivative_for_species(
                    candidate.aq_comp,
                    transform.aq_comp_dvar,
                    transform.aq_comp_hessian_row_major,
                    aq_sensitivity.fixed_density_jacobian_row_major,
                    aq_sensitivity.fixed_density_hessian_tensor_row_major,
                    aq_sensitivity.ln_fugacity_density_derivative,
                    aq_sensitivity.ln_fugacity_density_second_derivative,
                    aq_sensitivity.ln_fugacity_density_composition_cross_derivative,
                    candidate.aq_state.density,
                    electrolyte_density_variable_index(formula_nvar, false),
                    candidate.explicit_density,
                    species_index,
                    first,
                    second,
                    ncomp,
                    nvar
                );
                out.residual_hessian_tensor_row_major[
                    output_row * nvar * nvar + first * nvar + second
                ] += coefficient * (org_value - aq_value);
            }
        }
    };

    std::size_t row = 0;
    for (int index : basis.neutral_indices) {
        add_species_contribution(row, static_cast<std::size_t>(index), 1.0);
        ++row;
    }
    for (const auto& pair : basis.salt_pairs) {
        add_species_contribution(row, static_cast<std::size_t>(pair.cation), static_cast<double>(pair.cation_stoich));
        add_species_contribution(row, static_cast<std::size_t>(pair.anion), static_cast<double>(pair.anion_stoich));
        ++row;
    }
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t first = 0; first < nvar; ++first) {
            for (std::size_t second = 0; second < nvar; ++second) {
                double value =
                    -transform.beta_org_hessian_row_major[first * nvar + second] * candidate.aq_comp[species]
                    - transform.beta_org_dvar[first] * transform.aq_comp_dvar[species * nvar + second]
                    - transform.beta_org_dvar[second] * transform.aq_comp_dvar[species * nvar + first]
                    + (1.0 - candidate.beta_org)
                        * transform.aq_comp_hessian_row_major[species * nvar * nvar + first * nvar + second]
                    + transform.beta_org_hessian_row_major[first * nvar + second] * candidate.org_comp[species]
                    + transform.beta_org_dvar[first] * transform.org_comp_dvar[species * nvar + second]
                    + transform.beta_org_dvar[second] * transform.org_comp_dvar[species * nvar + first]
                    + candidate.beta_org
                        * transform.org_comp_hessian_row_major[species * nvar * nvar + first * nvar + second];
                out.residual_hessian_tensor_row_major[
                    row * nvar * nvar + first * nvar + second
                ] = value;
            }
        }
        ++row;
    }

    const std::size_t formula = static_cast<std::size_t>(separation_formula_index);
    for (std::size_t first = 0; first < nvar; ++first) {
        for (std::size_t second = 0; second < nvar; ++second) {
            out.phase_separation_hessian_row_major[first * nvar + second] = separation_sign * (
                transform.org_formula_hessian_row_major[formula * nvar * nvar + first * nvar + second]
                - transform.aq_formula_hessian_row_major[formula * nvar * nvar + first * nvar + second]
            );
        }
    }
    for (std::size_t formula_row = 0; formula_row < basis.formula_feed.size(); ++formula_row) {
        for (std::size_t first = 0; first < nvar; ++first) {
            for (std::size_t second = 0; second < nvar; ++second) {
                out.aqueous_formula_hessian_tensor_row_major[
                    formula_row * nvar * nvar + first * nvar + second
                ] = transform.aq_formula_hessian_row_major[
                    formula_row * nvar * nvar + first * nvar + second
                ];
            }
        }
    }
    return out;
}

std::vector<double> composition_to_logits(const std::vector<double>& composition) {
    std::vector<double> out(composition.size() - 1, 0.0);
    double denom = composition.back();
    for (std::size_t i = 0; i + 1 < composition.size(); ++i) {
        out[i] = std::log(composition[i] / denom);
    }
    return out;
}

std::vector<double> logits_to_composition(const std::vector<double>& logits) {
    std::vector<double> weights(logits.size() + 1, 1.0);
    double total = 1.0;
    for (std::size_t i = 0; i < logits.size(); ++i) {
        weights[i] = std::exp(std::max(-700.0, std::min(700.0, logits[i])));
        total += weights[i];
    }
    for (double& value : weights) {
        value /= total;
    }
    return weights;
}

double logit(double value) {
    value = std::max(1.0e-12, std::min(1.0 - 1.0e-12, value));
    return std::log(value / (1.0 - value));
}

}  // namespace

std::vector<double> pack_predictive_electrolyte_variables(double beta_formula, const std::vector<double>& org_formula) {
    beta_formula = std::max(1.0e-12, std::min(1.0 - 1.0e-12, beta_formula));
    std::vector<double> out;
    out.push_back(std::log(beta_formula / (1.0 - beta_formula)));
    std::vector<double> logits = composition_to_logits(org_formula);
    out.insert(out.end(), logits.begin(), logits.end());
    return out;
}

void unpack_predictive_electrolyte_variables(
    const std::vector<double>& variables,
    std::size_t nformula,
    double& beta_formula,
    std::vector<double>& org_formula
) {
    if (variables.size() != nformula) {
        throw SolutionError("Unexpected predictive electrolyte variable vector size.");
    }
    beta_formula = 1.0 / (1.0 + std::exp(-std::max(-700.0, std::min(700.0, variables[0]))));
    std::vector<double> logits(variables.begin() + 1, variables.end());
    org_formula = logits_to_composition(logits);
}

ElectrolyteCandidateNative build_predictive_electrolyte_formula_candidate(
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const std::vector<double>& formula_variables,
    const EquilibriumOptionsNative& options
) {
    ElectrolyteCandidateNative out;
    const std::size_t nformula = basis.formula_feed.size();
    if (formula_variables.size() != nformula) {
        throw SolutionError("Unexpected predictive electrolyte variable vector size.");
    }
    unpack_predictive_electrolyte_variables(formula_variables, nformula, out.beta_formula, out.org_formula);
    std::vector<double> aq_formula_raw(basis.formula_feed.size(), 0.0);
    for (std::size_t i = 0; i < basis.formula_feed.size(); ++i) {
        aq_formula_raw[i] = (basis.formula_feed[i] - out.beta_formula * out.org_formula[i]) / (1.0 - out.beta_formula);
        if (!std::isfinite(aq_formula_raw[i]) || aq_formula_raw[i] <= options.min_composition) {
            throw SolutionError("Predictive electrolyte variables produced an infeasible dependent phase.");
        }
    }
    out.aq_formula = clip_normalize(aq_formula_raw, options.min_composition);
    auto aq_expanded = formula_to_explicit(out.aq_formula, basis, feed.size());
    auto org_expanded = formula_to_explicit(out.org_formula, basis, feed.size());
    out.aq_comp = aq_expanded.first;
    out.org_comp = org_expanded.first;
    double aq_scale = aq_expanded.second;
    double org_scale = org_expanded.second;
    out.beta_org = out.beta_formula * org_scale / ((1.0 - out.beta_formula) * aq_scale + out.beta_formula * org_scale);
    return out;
}

ElectrolyteCandidateNative evaluate_predictive_electrolyte_variables(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const std::vector<double>& variables,
    const EquilibriumOptionsNative& options,
    double gibbs_feed
) {
    (void)p;
    const std::size_t nformula = basis.formula_feed.size();
    if (variables.size() != nformula + 2) {
        throw SolutionError("Predictive electrolyte variables must include transformed formula variables plus aqueous/org log-density variables.");
    }
    std::vector<double> formula_variables(
        variables.begin(),
        variables.begin() + static_cast<std::ptrdiff_t>(nformula)
    );
    ElectrolyteCandidateNative out = build_predictive_electrolyte_formula_candidate(
        feed,
        basis,
        formula_variables,
        options
    );
    out.explicit_density = true;
    const double aq_density = std::exp(variables[nformula]);
    const double org_density = std::exp(variables[nformula + 1]);
    out.aq_state = phase_state_at_density(mixture, t, aq_density, out.aq_comp, "liq");
    out.org_state = phase_state_at_density(mixture, t, org_density, out.org_comp, "liq");
    out.residual = electrolyte_residual_vector(out.aq_comp, out.org_comp, out.aq_state, out.org_state, basis);
    out.material_residual.resize(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        out.material_residual[i] = (1.0 - out.beta_org) * out.aq_comp[i] + out.beta_org * out.org_comp[i] - feed[i];
    }
    const std::vector<double>& charges = mixture->args().z;
    double feed_charge = composition_charge(feed, charges);
    double aq_charge = composition_charge(out.aq_comp, charges);
    double org_charge = composition_charge(out.org_comp, charges);
    out.solver_residual_norm = max_abs(out.residual);
    out.material_error = max_abs(out.material_residual);
    out.charge_balance_error = std::max({std::abs(feed_charge), std::abs(aq_charge), std::abs(org_charge)});
    out.gibbs_feed = gibbs_feed;
    out.gibbs_split = (1.0 - out.beta_org) * electrolyte_gibbs_proxy(out.aq_comp, out.aq_state)
        + out.beta_org * electrolyte_gibbs_proxy(out.org_comp, out.org_state);
    out.gibbs_delta = out.gibbs_split - out.gibbs_feed;
    out.phase_distance_value = phase_distance(out.aq_comp, out.org_comp);
    out.objective = 0.0;
    for (double value : out.residual) {
        out.objective += 0.5 * value * value;
    }
    for (double value : out.material_residual) {
        out.objective += 0.5 * value * value;
    }
    return out;
}

std::vector<double> explicit_density_seed_from_electrolyte_formula_variables(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const std::vector<double>& formula_variables,
    const EquilibriumOptionsNative& options
);

ElectrolyteLLEResidualEvaluationNative evaluate_electrolyte_lle_residual_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& raw_feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& species,
    const std::vector<double>& variables,
    bool has_variables,
    const std::vector<double>& initial_aq,
    const std::vector<double>& initial_org,
    double initial_beta_org,
    bool has_initial_phases
) {
    if (!mixture->has_ionic()) {
        throw ValueError("electrolyte_lle residual evaluation requires an ion-containing mixture.");
    }
    std::vector<double> feed = normalize_feed(raw_feed, mixture->ncomp(), options.min_composition, "electrolyte_lle");
    const std::vector<double>& charges = mixture->args().z;
    if (std::abs(composition_charge(feed, charges)) > 1.0e-10) {
        throw ValueError("electrolyte_lle residual evaluation feed must be charge neutral.");
    }
    ElectrolyteBasisNative basis = build_electrolyte_basis_native(mixture, feed, species);
    std::vector<double> eval_variables = variables;
    if (has_variables && variables.size() != basis.formula_feed.size() + 2) {
        throw SolutionError(
            "electrolyte_lle residual variables must include transformed formula variables plus aqueous/org log-density variables."
        );
    }
    if (!has_variables) {
        std::vector<double> org_formula = basis.formula_feed;
        double beta_formula = 0.5;
        if (has_initial_phases) {
            if (initial_aq.size() != feed.size() || initial_org.size() != feed.size()) {
                throw ValueError("initial_phases aq/org length must match mixture component count.");
            }
            if (!(initial_beta_org > 0.0 && initial_beta_org < 1.0) || !std::isfinite(initial_beta_org)) {
                throw ValueError("initial_phases phase_fraction must be > 0 and < 1.");
            }
            std::vector<double> aq_comp_initial = clip_normalize(initial_aq, options.min_composition);
            std::vector<double> org_comp_initial = clip_normalize(initial_org, options.min_composition);
            double aq_charge_initial = composition_charge(aq_comp_initial, charges);
            double org_charge_initial = composition_charge(org_comp_initial, charges);
            if (std::max(std::abs(aq_charge_initial), std::abs(org_charge_initial)) > 1.0e-8) {
                throw ValueError("initial_phases aq and org must be charge neutral for electrolyte_lle residual evaluation.");
            }
            std::vector<double> material_residual_initial(feed.size(), 0.0);
            for (std::size_t i = 0; i < feed.size(); ++i) {
                material_residual_initial[i] = (1.0 - initial_beta_org) * aq_comp_initial[i]
                    + initial_beta_org * org_comp_initial[i] - feed[i];
            }
            if (max_abs(material_residual_initial) > 1.0e-7) {
                throw ValueError("initial_phases aq/org/phase_fraction must reconstruct the electrolyte_lle feed.");
            }
            std::vector<double> aq_formula = explicit_to_formula(aq_comp_initial, basis);
            org_formula = explicit_to_formula(org_comp_initial, basis);
            auto aq_expanded = formula_to_explicit(aq_formula, basis, feed.size());
            auto org_expanded = formula_to_explicit(org_formula, basis, feed.size());
            double numerator = initial_beta_org / org_expanded.second;
            double denominator = numerator + (1.0 - initial_beta_org) / aq_expanded.second;
            if (denominator > 0.0) {
                beta_formula = std::max(1.0e-12, std::min(1.0 - 1.0e-12, numerator / denominator));
            }
        }
        const std::vector<double> formula_variables = pack_predictive_electrolyte_variables(beta_formula, org_formula);
        eval_variables = explicit_density_seed_from_electrolyte_formula_variables(
            mixture,
            t,
            p,
            feed,
            basis,
            formula_variables,
            options
        );
    }

    PhaseStateNative feed_state = phase_state(mixture, t, p, feed, "liq", "feed");
    double gibbs_feed = electrolyte_gibbs_proxy(feed, feed_state);
    ElectrolyteCandidateNative candidate = evaluate_predictive_electrolyte_variables(
        mixture,
        t,
        p,
        feed,
        basis,
        eval_variables,
        options,
        gibbs_feed
    );
    std::vector<double> residual = candidate.residual;
    residual.insert(residual.end(), candidate.material_residual.begin(), candidate.material_residual.end());
    double objective = 0.0;
    for (double value : residual) {
        objective += 0.5 * value * value;
    }
    std::vector<double> jacobian = electrolyte_residual_jacobian_row_major(mixture, t, p, feed, basis, candidate);
    std::vector<double> gradient(eval_variables.size(), 0.0);
    for (std::size_t row = 0; row < residual.size(); ++row) {
        for (std::size_t col = 0; col < eval_variables.size(); ++col) {
            gradient[col] += jacobian[row * eval_variables.size() + col] * residual[row];
        }
    }

    ElectrolyteLLEResidualEvaluationNative out;
    out.variable_model = basis.variable_model + "_explicit_density";
    out.variables = eval_variables;
    out.lower_bounds.assign(eval_variables.size(), -100.0);
    out.upper_bounds.assign(eval_variables.size(), 100.0);
    out.residual = residual;
    out.jacobian_rows = static_cast<int>(residual.size());
    out.jacobian_cols = static_cast<int>(eval_variables.size());
    out.objective = objective;
    out.gradient = gradient;
    out.jacobian_row_major = jacobian;
    out.aq_composition = candidate.aq_comp;
    out.org_composition = candidate.org_comp;
    out.aq_ln_fugacity_coefficient = candidate.aq_state.ln_phi;
    out.org_ln_fugacity_coefficient = candidate.org_state.ln_phi;
    out.aq_density = candidate.aq_state.density;
    out.org_density = candidate.org_state.density;
    out.phase_fraction_org = candidate.beta_org;
    out.material_balance_error = candidate.material_error;
    out.charge_balance_error = candidate.charge_balance_error;
    out.phase_distance = candidate.phase_distance_value;
    out.gibbs_delta = candidate.gibbs_delta;
    out.diagnostics_string["residual_surface"] = "native_electrolyte_lle_transformed_variables";
    out.diagnostics_string["residual_blocks"] = "phase_equilibrium,material_balance";
    out.diagnostics_string["jacobian_backend"] = "cppad_explicit_density";
    out.diagnostics_string["derivative_backend"] = "cppad_explicit_density";
    out.diagnostics_string["density_backend"] = "explicit_log_density_pressure_constraint";
    out.diagnostics_string["jacobian_scope"] = "transformed_variables_explicit_density";
    out.diagnostics_bool["jacobian_available"] = true;
    out.diagnostics_bool["derivative_available"] = true;
    out.diagnostics_bool["phase_charge_enforced_by_basis"] = true;
    out.diagnostics_bool["material_balance_enforced_by_formula_transform"] = true;
    out.diagnostics_int["basis_rank"] = basis.basis_rank;
    out.diagnostics_int["phase_equilibrium_residual_size"] = static_cast<int>(candidate.residual.size());
    out.diagnostics_int["material_balance_residual_size"] = static_cast<int>(candidate.material_residual.size());
    out.diagnostics_int["residual_size"] = static_cast<int>(residual.size());
    out.diagnostics_int["variable_count"] = static_cast<int>(eval_variables.size());
    out.diagnostics_double["phase_equilibrium_residual_norm"] = candidate.solver_residual_norm;
    out.diagnostics_double["solver_residual_norm"] = max_abs(residual);
    out.diagnostics_double["fugacity_residual_norm"] = candidate.solver_residual_norm;
    out.diagnostics_double["material_balance_error"] = candidate.material_error;
    out.diagnostics_double["charge_balance_error"] = candidate.charge_balance_error;
    out.diagnostics_double["phase_distance"] = candidate.phase_distance_value;
    out.diagnostics_double["gibbs_feed"] = candidate.gibbs_feed;
    out.diagnostics_double["gibbs_split"] = candidate.gibbs_split;
    out.diagnostics_double["gibbs_delta"] = candidate.gibbs_delta;
    out.diagnostics_double["phase_charge_balance_feed"] = composition_charge(feed, basis.species_charges);
    out.diagnostics_double["phase_charge_balance_aq"] = composition_charge(candidate.aq_comp, basis.species_charges);
    out.diagnostics_double["phase_charge_balance_org"] = composition_charge(candidate.org_comp, basis.species_charges);
    out.diagnostics_vector["feed_composition"] = feed;
    out.diagnostics_vector["phase_equilibrium_residual"] = candidate.residual;
    out.diagnostics_vector["material_balance_residual"] = candidate.material_residual;
    out.diagnostics_vector["formula_feed"] = basis.formula_feed;
    out.diagnostics_vector["aq_formula"] = candidate.aq_formula;
    out.diagnostics_vector["org_formula"] = candidate.org_formula;
    out.diagnostics_vector["basis_vectors_row_major"] = electrolyte_basis_vectors_row_major(basis, feed.size());
    out.diagnostics_vector["species_charge_vector"] = basis.species_charges;
    return out;
}

std::vector<double> cap_formula_seed_to_material_feasible_region(
    std::vector<double> formula,
    const std::vector<double>& feed_formula,
    double beta_formula
) {
    const double safety = 0.995;
    std::vector<double> caps;
    caps.reserve(feed_formula.size());
    for (double value : feed_formula) {
        caps.push_back(std::min(1.0, safety * value / beta_formula));
    }
    formula = clip_normalize(formula, 1.0e-300);
    for (int iter = 0; iter < 24; ++iter) {
        double excess = 0.0;
        std::vector<bool> capped(formula.size(), false);
        for (std::size_t i = 0; i < formula.size(); ++i) {
            if (formula[i] > caps[i]) {
                excess += formula[i] - caps[i];
                formula[i] = caps[i];
                capped[i] = true;
            }
        }
        if (excess <= 1.0e-14) {
            break;
        }
        double capacity = 0.0;
        for (std::size_t i = 0; i < formula.size(); ++i) {
            if (!capped[i]) {
                capacity += std::max(0.0, caps[i] - formula[i]);
            }
        }
        if (capacity <= 0.0) {
            break;
        }
        for (std::size_t i = 0; i < formula.size(); ++i) {
            if (!capped[i]) {
                formula[i] += excess * std::max(0.0, caps[i] - formula[i]) / capacity;
            }
        }
    }
    return clip_normalize(formula, 1.0e-300);
}

std::vector<double> build_liquid_root_electrolyte_lle_initial_variables(
    const ElectrolyteBasisNative& basis,
    double shift_sign = 1.0
) {
    const double beta_formula = shift_sign > 0.0 ? 0.635 : 0.365;
    std::vector<double> org_formula = basis.formula_feed;
    if (basis.formula_feed.size() >= 3 && basis.neutral_indices.size() >= 2) {
        org_formula[0] *= shift_sign > 0.0 ? 0.667 : 1.548;
        org_formula[1] *= shift_sign > 0.0 ? 1.548 : 0.667;
        for (std::size_t pos = basis.neutral_indices.size(); pos < org_formula.size(); ++pos) {
            org_formula[pos] *= 0.106;
        }
    } else if (basis.formula_feed.size() >= 2) {
        org_formula[0] *= shift_sign > 0.0 ? 0.75 : 1.25;
        org_formula[1] *= shift_sign > 0.0 ? 1.25 : 0.75;
    }
    org_formula = cap_formula_seed_to_material_feasible_region(org_formula, basis.formula_feed, beta_formula);
    return pack_predictive_electrolyte_variables(beta_formula, org_formula);
}

std::vector<double> append_electrolyte_phase_density_variables(
    std::vector<double> variables,
    double aqueous_density,
    double organic_density
) {
    if (!std::isfinite(aqueous_density) || aqueous_density <= 0.0
        || !std::isfinite(organic_density) || organic_density <= 0.0) {
        throw SolutionError("Electrolyte LLE explicit-density initial point requires positive finite phase densities.");
    }
    variables.push_back(std::log(aqueous_density));
    variables.push_back(std::log(organic_density));
    return variables;
}

struct NamedInitialVariables {
    std::string seed_name;
    std::vector<double> variables;
};

std::vector<double> org_formula_jacobian_row_major(
    const std::vector<double>& org_formula,
    std::size_t nvar
) {
    const std::size_t nformula = org_formula.size();
    std::vector<double> jacobian(nformula * nvar, 0.0);
    for (std::size_t component = 0; component < nformula; ++component) {
        for (std::size_t var = 1; var < std::min(nvar, nformula); ++var) {
            const std::size_t logit_component = var - 1;
            const double indicator = component == logit_component ? 1.0 : 0.0;
            jacobian[component * nvar + var] =
                org_formula[component] * (indicator - org_formula[logit_component]);
        }
    }
    return jacobian;
}

std::vector<double> aqueous_formula_feasibility_values(
    const std::vector<double>& variables,
    const ElectrolyteBasisNative& basis
) {
    double beta_formula = 0.5;
    std::vector<double> org_formula;
    std::vector<double> formula_variables(
        variables.begin(),
        variables.begin() + static_cast<std::ptrdiff_t>(basis.formula_feed.size())
    );
    unpack_predictive_electrolyte_variables(formula_variables, basis.formula_feed.size(), beta_formula, org_formula);
    std::vector<double> out;
    out.reserve(basis.formula_feed.size());
    for (std::size_t i = 0; i < basis.formula_feed.size(); ++i) {
        out.push_back((basis.formula_feed[i] - beta_formula * org_formula[i]) / (1.0 - beta_formula));
    }
    return out;
}

std::vector<double> aqueous_formula_feasibility_jacobian_row_major(
    const std::vector<double>& variables,
    const ElectrolyteBasisNative& basis
) {
    double beta_formula = 0.5;
    std::vector<double> org_formula;
    std::vector<double> formula_variables(
        variables.begin(),
        variables.begin() + static_cast<std::ptrdiff_t>(basis.formula_feed.size())
    );
    unpack_predictive_electrolyte_variables(formula_variables, basis.formula_feed.size(), beta_formula, org_formula);
    const std::size_t nvar = variables.size();
    const std::size_t nformula = basis.formula_feed.size();
    std::vector<double> jacobian(nformula * nvar, 0.0);
    const double dbeta = beta_formula * (1.0 - beta_formula);
    const double denom = 1.0 - beta_formula;
    const std::vector<double> dorg = org_formula_jacobian_row_major(org_formula, nvar);
    for (std::size_t component = 0; component < nformula; ++component) {
        jacobian[component * nvar] =
            dbeta * (basis.formula_feed[component] - org_formula[component]) / (denom * denom);
        for (std::size_t var = 1; var < nvar; ++var) {
            jacobian[component * nvar + var] =
                -beta_formula * dorg[component * nvar + var] / denom;
        }
    }
    return jacobian;
}

std::vector<double> electrolyte_lle_residual_vector_from_candidate(
    const ElectrolyteCandidateNative& candidate
) {
    std::vector<double> residual = candidate.residual;
    residual.insert(residual.end(), candidate.material_residual.begin(), candidate.material_residual.end());
    return residual;
}

std::vector<double> electrolyte_lle_objective_gradient_from_candidate(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const ElectrolyteCandidateNative& candidate
) {
    const std::vector<double> residual = electrolyte_lle_residual_vector_from_candidate(candidate);
    const std::vector<double> jacobian = electrolyte_residual_jacobian_row_major(mixture, t, p, feed, basis, candidate);
    const std::size_t cols = basis.formula_feed.size() + (candidate.explicit_density ? 2 : 0);
    std::vector<double> gradient(cols, 0.0);
    for (std::size_t row = 0; row < residual.size(); ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
            gradient[col] += jacobian[row * cols + col] * residual[row];
        }
    }
    return gradient;
}

std::vector<double> electrolyte_pressure_constraint_jacobian(
    const PhaseStateCompositionSensitivityResult& sensitivity,
    const std::vector<double>& composition_dvar,
    double density,
    std::size_t density_var,
    double scale,
    std::size_t ncomp,
    std::size_t nvar
) {
    std::vector<double> jacobian(nvar, 0.0);
    for (std::size_t var = 0; var < nvar; ++var) {
        for (std::size_t species = 0; species < ncomp; ++species) {
            jacobian[var] += sensitivity.pressure_composition_fixed_density_derivative[species]
                * composition_dvar[species * nvar + var];
        }
        if (var == density_var) {
            jacobian[var] += density * sensitivity.pressure_density_derivative;
        }
        jacobian[var] /= scale;
    }
    return jacobian;
}

std::vector<double> electrolyte_pressure_constraint_hessian(
    const PhaseStateCompositionSensitivityResult& sensitivity,
    const std::vector<double>& composition_dvar,
    const std::vector<double>& composition_hessian,
    double density,
    std::size_t density_var,
    double scale,
    std::size_t ncomp,
    std::size_t nvar
) {
    std::vector<double> hessian(nvar * nvar, 0.0);
    for (std::size_t first = 0; first < nvar; ++first) {
        for (std::size_t second = 0; second < nvar; ++second) {
            double value = 0.0;
            for (std::size_t species = 0; species < ncomp; ++species) {
                value += sensitivity.pressure_composition_fixed_density_derivative[species]
                    * composition_hessian[species * nvar * nvar + first * nvar + second];
                for (std::size_t other = 0; other < ncomp; ++other) {
                    value += sensitivity.pressure_composition_fixed_density_hessian_row_major[species * ncomp + other]
                        * composition_dvar[species * nvar + first]
                        * composition_dvar[other * nvar + second];
                }
                if (first == density_var) {
                    value += density
                        * sensitivity.pressure_density_composition_cross_derivative[species]
                        * composition_dvar[species * nvar + second];
                }
                if (second == density_var) {
                    value += density
                        * sensitivity.pressure_density_composition_cross_derivative[species]
                        * composition_dvar[species * nvar + first];
                }
            }
            if (first == density_var && second == density_var) {
                value += density * density * sensitivity.pressure_density_second_derivative
                    + density * sensitivity.pressure_density_derivative;
            }
            hessian[first * nvar + second] = value / scale;
        }
    }
    return hessian;
}

class LiquidRootElectrolyteLleProblem final : public equilibrium_nlp::NlpProblem {
public:
    LiquidRootElectrolyteLleProblem(
        std::shared_ptr<ePCSAFTMixtureNative> mixture,
        double temperature,
        double target_pressure,
        std::vector<double> raw_feed,
        EquilibriumOptionsNative options,
        std::vector<std::string> species,
        double phase_distance_tolerance
    )
        : mixture_(std::move(mixture)),
          temperature_(temperature),
          target_pressure_(target_pressure),
          options_(std::move(options)),
          species_(std::move(species)),
          minimum_phase_distance_(std::max(phase_distance_tolerance, kLiquidLleMinimumPhaseDistance)) {
        if (!mixture_) {
            throw ValueError("Electrolyte LLE liquid-root NLP requires a native mixture.");
        }
        if (!mixture_->has_ionic()) {
            throw ValueError("Electrolyte LLE liquid-root NLP requires an ion-containing mixture.");
        }
        feed_ = normalize_feed(raw_feed, mixture_->ncomp(), options_.min_composition, "electrolyte_lle");
        const std::vector<double>& charges = mixture_->args().z;
        if (std::abs(composition_charge(feed_, charges)) > 1.0e-10) {
            throw ValueError("Electrolyte LLE liquid-root NLP feed must be charge neutral.");
        }
        basis_ = build_electrolyte_basis_native(mixture_, feed_, species_);
        PhaseStateNative feed_state = phase_state(mixture_, temperature_, target_pressure_, feed_, "liq", "feed");
        gibbs_feed_ = electrolyte_gibbs_proxy(feed_, feed_state);
        const std::vector<double> formula_initial = build_liquid_root_electrolyte_lle_initial_variables(basis_);
        initial_variables_ = explicit_density_seed_from_electrolyte_formula_variables(
            mixture_,
            temperature_,
            target_pressure_,
            feed_,
            basis_,
            formula_initial,
            options_
        );
        const ElectrolyteCandidateNative initial = candidate(initial_variables_);
        residual_constraint_count_ = initial.residual.size();
        initialize_pressure_scales(initial);
        select_separation_component(initial);
    }

    std::string name() const override {
        return "electrolyte_lle_eos";
    }

    int variable_count() const override {
        return static_cast<int>(initial_variables_.size());
    }

    int constraint_count() const override {
        return static_cast<int>(
            residual_constraint_count_ + 2 + 1 + basis_.formula_feed.size()
        );
    }

    int jacobian_nonzero_count() const override {
        return variable_count() * constraint_count();
    }

    equilibrium_nlp::NlpBounds bounds() const override {
        equilibrium_nlp::NlpBounds out;
        out.variable_lower.assign(initial_variables_.size(), -100.0);
        out.variable_upper.assign(initial_variables_.size(), 100.0);
        out.variable_lower[0] = logit(kLiquidLleMinimumRetainedPhaseFraction);
        out.variable_upper[0] = logit(1.0 - kLiquidLleMinimumRetainedPhaseFraction);
        out.variable_lower[aqueous_density_variable_index()] = std::log(aqueous_density_lower_bound_);
        out.variable_upper[aqueous_density_variable_index()] = std::log(aqueous_density_upper_bound_);
        out.variable_lower[organic_density_variable_index()] = std::log(organic_density_lower_bound_);
        out.variable_upper[organic_density_variable_index()] = std::log(organic_density_upper_bound_);
        out.constraint_lower.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        const std::size_t phase_distance_row = phase_distance_constraint_index();
        out.constraint_lower[phase_distance_row] = minimum_phase_distance_;
        out.constraint_upper[phase_distance_row] = 1.0e12;
        for (std::size_t row = phase_distance_row + 1; row < out.constraint_lower.size(); ++row) {
            out.constraint_lower[row] = options_.min_composition;
            out.constraint_upper[row] = 1.0e12;
        }
        return out;
    }

    std::vector<double> initial_point() const override {
        return initial_variables_;
    }

    double objective(const std::vector<double>& variables) const override {
        try {
            return candidate(variables).objective;
        } catch (const std::exception&) {
            return infeasible_objective_penalty(variables);
        }
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        try {
            return electrolyte_lle_objective_gradient_from_candidate(
                mixture_,
                temperature_,
                target_pressure_,
                feed_,
                basis_,
                candidate(variables)
            );
        } catch (const std::exception&) {
            return infeasible_objective_penalty_gradient(variables);
        }
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(constraint_count()));
        try {
            const ElectrolyteCandidateNative current = candidate(variables);
            out.insert(out.end(), current.residual.begin(), current.residual.end());
            out.push_back(pressure_constraint_value(current, false));
            out.push_back(pressure_constraint_value(current, true));
        } catch (const std::exception&) {
            out.insert(out.end(), residual_constraint_count_, 1.0e3);
            out.push_back(1.0e3);
            out.push_back(1.0e3);
        }
        out.push_back(phase_separation_from_variables(variables));
        const std::vector<double> aq_formula = aqueous_formula_feasibility_values(variables, basis_);
        out.insert(out.end(), aq_formula.begin(), aq_formula.end());
        return out;
    }

    equilibrium_nlp::NlpJacobianStructure jacobian_structure() const override {
        equilibrium_nlp::NlpJacobianStructure out;
        out.rows.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        out.cols.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < constraint_count(); ++row) {
            for (int col = 0; col < variable_count(); ++col) {
                out.rows.push_back(row);
                out.cols.push_back(col);
            }
        }
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        const std::size_t nvar = static_cast<std::size_t>(variable_count());
        std::vector<double> out(static_cast<std::size_t>(constraint_count()) * nvar, 0.0);
        try {
            const ElectrolyteCandidateNative current = candidate(variables);
            const ElectrolyteTransformDerivativesNative transform =
                electrolyte_transform_derivatives(feed_, basis_, current, nvar);
            const std::vector<double> residual_jacobian = electrolyte_residual_jacobian_row_major(
                mixture_,
                temperature_,
                target_pressure_,
                feed_,
                basis_,
                current
            );
            for (std::size_t row = 0; row < residual_constraint_count_; ++row) {
                for (std::size_t col = 0; col < nvar; ++col) {
                    out[row * nvar + col] = residual_jacobian[row * nvar + col];
                }
            }
            const PhaseStateCompositionSensitivityResult aq_sensitivity = electrolyte_phase_sensitivity(
                mixture_,
                temperature_,
                target_pressure_,
                current.aq_comp,
                current.aq_state.density,
                true
            );
            const PhaseStateCompositionSensitivityResult org_sensitivity = electrolyte_phase_sensitivity(
                mixture_,
                temperature_,
                target_pressure_,
                current.org_comp,
                current.org_state.density,
                true
            );
            const std::vector<double> aq_pressure_jacobian = electrolyte_pressure_constraint_jacobian(
                aq_sensitivity,
                transform.aq_comp_dvar,
                current.aq_state.density,
                aqueous_density_variable_index(),
                aqueous_pressure_constraint_scale_,
                feed_.size(),
                nvar
            );
            const std::vector<double> org_pressure_jacobian = electrolyte_pressure_constraint_jacobian(
                org_sensitivity,
                transform.org_comp_dvar,
                current.org_state.density,
                organic_density_variable_index(),
                organic_pressure_constraint_scale_,
                feed_.size(),
                nvar
            );
            std::copy(
                aq_pressure_jacobian.begin(),
                aq_pressure_jacobian.end(),
                out.begin() + static_cast<std::ptrdiff_t>(aqueous_pressure_constraint_index() * nvar)
            );
            std::copy(
                org_pressure_jacobian.begin(),
                org_pressure_jacobian.end(),
                out.begin() + static_cast<std::ptrdiff_t>(organic_pressure_constraint_index() * nvar)
            );
        } catch (const std::exception&) {
        }
        double beta_formula = 0.5;
        std::vector<double> org_formula;
        const std::vector<double> formula_variables(
            variables.begin(),
            variables.begin() + static_cast<std::ptrdiff_t>(basis_.formula_feed.size())
        );
        unpack_predictive_electrolyte_variables(formula_variables, basis_.formula_feed.size(), beta_formula, org_formula);
        const std::vector<double> org_jacobian = org_formula_jacobian_row_major(org_formula, nvar);
        const std::vector<double> aq_jacobian = aqueous_formula_feasibility_jacobian_row_major(variables, basis_);
        const std::size_t formula = static_cast<std::size_t>(separation_formula_index_);
        const std::size_t phase_distance_row = phase_distance_constraint_index();
        for (std::size_t col = 0; col < nvar; ++col) {
            out[phase_distance_row * nvar + col] = separation_sign_ * (
                org_jacobian[formula * nvar + col]
                - aq_jacobian[formula * nvar + col]
            );
        }
        for (std::size_t row = 0; row < basis_.formula_feed.size(); ++row) {
            for (std::size_t col = 0; col < nvar; ++col) {
                out[(phase_distance_row + 1 + row) * nvar + col] = aq_jacobian[row * nvar + col];
            }
        }
        return out;
    }

    bool has_exact_hessian() const override {
        return true;
    }

    int hessian_nonzero_count() const override {
        return equilibrium_nlp::LagrangianHessianAssembler(variable_count()).nonzero_count();
    }

    equilibrium_nlp::NlpHessianStructure hessian_structure() const override {
        return equilibrium_nlp::LagrangianHessianAssembler(variable_count()).structure();
    }

    std::vector<double> hessian_values(
        const std::vector<double>& variables,
        double objective_factor,
        const std::vector<double>& constraint_multipliers
    ) const override {
        const std::size_t nvar = static_cast<std::size_t>(variable_count());
        if (constraint_multipliers.size() != static_cast<std::size_t>(constraint_count())) {
            throw ValueError("Electrolyte LLE liquid-root NLP Hessian received an unexpected multiplier vector size.");
        }
        equilibrium_nlp::ObjectiveSecondOrderData objective;
        objective.variable_count = variable_count();
        objective.hessian_row_major.assign(nvar * nvar, 0.0);
        objective.backend = "cppad_explicit_density_chain_rule";
        equilibrium_nlp::ConstraintSecondOrderData constraints;
        constraints.constraint_count = constraint_count();
        constraints.variable_count = variable_count();
        constraints.hessian_tensor_row_major.assign(
            static_cast<std::size_t>(constraints.constraint_count) * nvar * nvar,
            0.0
        );
        constraints.has_hessian.assign(static_cast<std::size_t>(constraints.constraint_count), false);
        constraints.backend = "cppad_explicit_density_chain_rule";
        try {
            const ElectrolyteCandidateNative current = candidate(variables);
            const ElectrolyteTransformDerivativesNative transform =
                electrolyte_transform_derivatives(feed_, basis_, current, nvar);
            const ElectrolyteResidualSecondOrderNative second_order = electrolyte_residual_second_order_native(
                mixture_,
                temperature_,
                target_pressure_,
                feed_,
                basis_,
                current,
                separation_formula_index_,
                separation_sign_
            );
            const std::vector<double> residual = electrolyte_lle_residual_vector_from_candidate(current);
            if (second_order.residual_rows != residual.size() || second_order.nvar != nvar) {
                throw ValueError("Electrolyte LLE liquid-root Hessian helper returned an unexpected shape.");
            }
            equilibrium_nlp::ResidualSecondOrderData residual_second_order;
            residual_second_order.residual_count = static_cast<int>(residual.size());
            residual_second_order.variable_count = variable_count();
            residual_second_order.residuals = residual;
            residual_second_order.jacobian_row_major = second_order.jacobian_row_major;
            residual_second_order.residual_hessian_tensor_row_major =
                second_order.residual_hessian_tensor_row_major;
            residual_second_order.backend = "cppad_explicit_density_chain_rule";
            objective = equilibrium_nlp::residual_quadratic_objective_second_order(residual_second_order);
            for (std::size_t row = 0; row < current.residual.size(); ++row) {
                constraints.has_hessian[row] = true;
                std::copy(
                    second_order.residual_hessian_tensor_row_major.begin()
                        + static_cast<std::ptrdiff_t>(row * nvar * nvar),
                    second_order.residual_hessian_tensor_row_major.begin()
                        + static_cast<std::ptrdiff_t>((row + 1) * nvar * nvar),
                    constraints.hessian_tensor_row_major.begin()
                        + static_cast<std::ptrdiff_t>(row * nvar * nvar)
                );
            }
            const PhaseStateCompositionSensitivityResult aq_sensitivity = electrolyte_phase_sensitivity(
                mixture_,
                temperature_,
                target_pressure_,
                current.aq_comp,
                current.aq_state.density,
                true
            );
            const PhaseStateCompositionSensitivityResult org_sensitivity = electrolyte_phase_sensitivity(
                mixture_,
                temperature_,
                target_pressure_,
                current.org_comp,
                current.org_state.density,
                true
            );
            const std::vector<double> aq_pressure_hessian = electrolyte_pressure_constraint_hessian(
                aq_sensitivity,
                transform.aq_comp_dvar,
                transform.aq_comp_hessian_row_major,
                current.aq_state.density,
                aqueous_density_variable_index(),
                aqueous_pressure_constraint_scale_,
                feed_.size(),
                nvar
            );
            const std::vector<double> org_pressure_hessian = electrolyte_pressure_constraint_hessian(
                org_sensitivity,
                transform.org_comp_dvar,
                transform.org_comp_hessian_row_major,
                current.org_state.density,
                organic_density_variable_index(),
                organic_pressure_constraint_scale_,
                feed_.size(),
                nvar
            );
            constraints.has_hessian[aqueous_pressure_constraint_index()] = true;
            std::copy(
                aq_pressure_hessian.begin(),
                aq_pressure_hessian.end(),
                constraints.hessian_tensor_row_major.begin()
                    + static_cast<std::ptrdiff_t>(aqueous_pressure_constraint_index() * nvar * nvar)
            );
            constraints.has_hessian[organic_pressure_constraint_index()] = true;
            std::copy(
                org_pressure_hessian.begin(),
                org_pressure_hessian.end(),
                constraints.hessian_tensor_row_major.begin()
                    + static_cast<std::ptrdiff_t>(organic_pressure_constraint_index() * nvar * nvar)
            );
            const std::size_t phase_distance_row = phase_distance_constraint_index();
            constraints.has_hessian[phase_distance_row] = true;
            std::copy(
                second_order.phase_separation_hessian_row_major.begin(),
                second_order.phase_separation_hessian_row_major.end(),
                constraints.hessian_tensor_row_major.begin()
                    + static_cast<std::ptrdiff_t>(phase_distance_row * nvar * nvar)
            );
            for (std::size_t formula_row = 0; formula_row < basis_.formula_feed.size(); ++formula_row) {
                const std::size_t constraint_row = phase_distance_row + 1 + formula_row;
                constraints.has_hessian[constraint_row] = true;
                std::copy(
                    second_order.aqueous_formula_hessian_tensor_row_major.begin()
                        + static_cast<std::ptrdiff_t>(formula_row * nvar * nvar),
                    second_order.aqueous_formula_hessian_tensor_row_major.begin()
                        + static_cast<std::ptrdiff_t>((formula_row + 1) * nvar * nvar),
                    constraints.hessian_tensor_row_major.begin()
                        + static_cast<std::ptrdiff_t>(constraint_row * nvar * nvar)
                );
            }
        } catch (const std::exception&) {
            const ElectrolyteTransformDerivativesNative transform = infeasible_formula_transform_derivatives(variables);
            const std::vector<double> aq_formula = aqueous_formula_feasibility_values(variables, basis_);
            const double scale = 1.0e8;
            objective.backend = "analytic_infeasible_transform";
            constraints.backend = "analytic_infeasible_transform";
            for (std::size_t formula_row = 0; formula_row < aq_formula.size(); ++formula_row) {
                const double violation = std::max(0.0, options_.min_composition - aq_formula[formula_row]);
                if (violation <= 0.0) {
                    continue;
                }
                for (std::size_t first = 0; first < nvar; ++first) {
                    for (std::size_t second = 0; second < nvar; ++second) {
                        objective.hessian_row_major[first * nvar + second] += scale * (
                            transform.aq_formula_dvar[formula_row * nvar + first]
                            * transform.aq_formula_dvar[formula_row * nvar + second]
                            - violation * transform.aq_formula_hessian_row_major[
                                formula_row * nvar * nvar + first * nvar + second
                            ]
                        );
                    }
                }
            }
            const std::size_t phase_distance_row = phase_distance_constraint_index();
            const std::size_t formula = static_cast<std::size_t>(separation_formula_index_);
            constraints.has_hessian[phase_distance_row] = true;
            for (std::size_t first = 0; first < nvar; ++first) {
                for (std::size_t second = 0; second < nvar; ++second) {
                    constraints.hessian_tensor_row_major[
                        phase_distance_row * nvar * nvar + first * nvar + second
                    ] = separation_sign_ * (
                        transform.org_formula_hessian_row_major[formula * nvar * nvar + first * nvar + second]
                        - transform.aq_formula_hessian_row_major[formula * nvar * nvar + first * nvar + second]
                    );
                }
            }
            for (std::size_t formula_row = 0; formula_row < basis_.formula_feed.size(); ++formula_row) {
                const std::size_t constraint_row = phase_distance_row + 1 + formula_row;
                constraints.has_hessian[constraint_row] = true;
                std::copy(
                    transform.aq_formula_hessian_row_major.begin()
                        + static_cast<std::ptrdiff_t>(formula_row * nvar * nvar),
                    transform.aq_formula_hessian_row_major.begin()
                        + static_cast<std::ptrdiff_t>((formula_row + 1) * nvar * nvar),
                    constraints.hessian_tensor_row_major.begin()
                        + static_cast<std::ptrdiff_t>(constraint_row * nvar * nvar)
                );
            }
        }

        return equilibrium_nlp::LagrangianHessianAssembler(variable_count()).values(
            objective_factor,
            objective,
            constraints,
            constraint_multipliers
        );
    }

    std::string hessian_backend() const override {
        return "cppad_explicit_density_chain_rule";
    }

    equilibrium_nlp::NlpScaling scaling() const override {
        equilibrium_nlp::NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(initial_variables_.size(), 1.0);
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0);
        return out;
    }

    std::map<std::string, std::string> diagnostics() const override {
        return {
            {"derivative_backend", "cppad_explicit_density"},
            {"density_backend", "explicit_log_density_pressure_constraint"},
            {"variable_model", basis_.variable_model + "_explicit_density"},
        };
    }

    ElectrolyteCandidateNative candidate(const std::vector<double>& variables) const {
        return evaluate_predictive_electrolyte_variables(
            mixture_,
            temperature_,
            target_pressure_,
            feed_,
            basis_,
            variables,
            options_,
            gibbs_feed_
        );
    }

    std::vector<double> explicit_density_seed_from_formula_variables(
        const std::vector<double>& formula_variables
    ) const {
        return explicit_density_seed_from_electrolyte_formula_variables(
            mixture_,
            temperature_,
            target_pressure_,
            feed_,
            basis_,
            formula_variables,
            options_
        );
    }

    const ElectrolyteBasisNative& basis() const {
        return basis_;
    }

    const std::vector<double>& feed() const {
        return feed_;
    }

    double minimum_phase_distance() const {
        return minimum_phase_distance_;
    }

    double minimum_phase_fraction() const {
        return kLiquidLleMinimumRetainedPhaseFraction;
    }

    double phase_separation(const ElectrolyteCandidateNative& current) const {
        const std::size_t species = static_cast<std::size_t>(separation_species_index_);
        return separation_sign_ * (current.org_comp[species] - current.aq_comp[species]);
    }

    double phase_separation_from_variables(const std::vector<double>& variables) const {
        double beta_formula = 0.5;
        std::vector<double> org_formula;
        const std::vector<double> formula_variables(
            variables.begin(),
            variables.begin() + static_cast<std::ptrdiff_t>(basis_.formula_feed.size())
        );
        unpack_predictive_electrolyte_variables(formula_variables, basis_.formula_feed.size(), beta_formula, org_formula);
        const std::vector<double> aq_formula = aqueous_formula_feasibility_values(variables, basis_);
        const std::size_t formula = static_cast<std::size_t>(separation_formula_index_);
        return separation_sign_ * (org_formula[formula] - aq_formula[formula]);
    }

    double pressure_consistency_norm(const ElectrolyteCandidateNative& current) const {
        return std::max(
            std::abs(pressure_constraint_value(current, false) * aqueous_pressure_constraint_scale_),
            std::abs(pressure_constraint_value(current, true) * organic_pressure_constraint_scale_)
        );
    }

    double pressure_consistency_tolerance() const {
        return std::max(1.0e-2, 1.0e-8 * std::abs(target_pressure_));
    }

private:
    std::size_t aqueous_density_variable_index() const {
        return basis_.formula_feed.size();
    }

    std::size_t organic_density_variable_index() const {
        return basis_.formula_feed.size() + 1;
    }

    std::size_t aqueous_pressure_constraint_index() const {
        return residual_constraint_count_;
    }

    std::size_t organic_pressure_constraint_index() const {
        return residual_constraint_count_ + 1;
    }

    std::size_t phase_distance_constraint_index() const {
        return residual_constraint_count_ + 2;
    }

    void initialize_pressure_scales(const ElectrolyteCandidateNative& initial) {
        const PhaseStateCompositionSensitivityResult aq_sensitivity = electrolyte_phase_sensitivity(
            mixture_,
            temperature_,
            target_pressure_,
            initial.aq_comp,
            initial.aq_state.density,
            true
        );
        const PhaseStateCompositionSensitivityResult org_sensitivity = electrolyte_phase_sensitivity(
            mixture_,
            temperature_,
            target_pressure_,
            initial.org_comp,
            initial.org_state.density,
            true
        );
        aqueous_pressure_constraint_scale_ = electrolyte_pressure_constraint_scale(
            target_pressure_,
            initial.aq_state.density,
            aq_sensitivity
        );
        organic_pressure_constraint_scale_ = electrolyte_pressure_constraint_scale(
            target_pressure_,
            initial.org_state.density,
            org_sensitivity
        );
        aqueous_density_lower_bound_ = std::max(minimum_density_, initial.aq_state.density / density_bound_factor_);
        aqueous_density_upper_bound_ = std::min(maximum_density_, initial.aq_state.density * density_bound_factor_);
        organic_density_lower_bound_ = std::max(minimum_density_, initial.org_state.density / density_bound_factor_);
        organic_density_upper_bound_ = std::min(maximum_density_, initial.org_state.density * density_bound_factor_);
    }

    double pressure_constraint_value(const ElectrolyteCandidateNative& current, bool organic_phase) const {
        const std::vector<double>& composition = organic_phase ? current.org_comp : current.aq_comp;
        const double density = organic_phase ? current.org_state.density : current.aq_state.density;
        const double scale = organic_phase ? organic_pressure_constraint_scale_ : aqueous_pressure_constraint_scale_;
        const PhaseStateCompositionSensitivityResult sensitivity = electrolyte_phase_sensitivity(
            mixture_,
            temperature_,
            target_pressure_,
            composition,
            density,
            true
        );
        return (sensitivity.pressure - target_pressure_) / scale;
    }

    double infeasible_objective_penalty(const std::vector<double>& variables) const {
        const double scale = 1.0e8;
        double penalty = 1.0e6;
        const std::vector<double> aq_formula = aqueous_formula_feasibility_values(variables, basis_);
        for (double value : aq_formula) {
            const double violation = std::max(0.0, options_.min_composition - value);
            penalty += 0.5 * scale * violation * violation;
        }
        return std::isfinite(penalty) ? penalty : 1.0e100;
    }

    std::vector<double> infeasible_objective_penalty_gradient(const std::vector<double>& variables) const {
        const double scale = 1.0e8;
        const std::size_t nvar = static_cast<std::size_t>(variable_count());
        std::vector<double> gradient(nvar, 0.0);
        const std::vector<double> aq_formula = aqueous_formula_feasibility_values(variables, basis_);
        const std::vector<double> aq_jacobian = aqueous_formula_feasibility_jacobian_row_major(variables, basis_);
        for (std::size_t row = 0; row < aq_formula.size(); ++row) {
            const double violation = std::max(0.0, options_.min_composition - aq_formula[row]);
            if (violation <= 0.0) {
                continue;
            }
            for (std::size_t col = 0; col < nvar; ++col) {
                gradient[col] -= scale * violation * aq_jacobian[row * nvar + col];
            }
        }
        return gradient;
    }

    ElectrolyteTransformDerivativesNative infeasible_formula_transform_derivatives(
        const std::vector<double>& variables
    ) const {
        double beta_formula = 0.5;
        std::vector<double> org_formula;
        const std::vector<double> formula_variables(
            variables.begin(),
            variables.begin() + static_cast<std::ptrdiff_t>(basis_.formula_feed.size())
        );
        unpack_predictive_electrolyte_variables(formula_variables, basis_.formula_feed.size(), beta_formula, org_formula);
        const std::size_t ncomp = feed_.size();
        const std::size_t nformula = basis_.formula_feed.size();
        const std::size_t nvar = variables.size();
        ElectrolyteTransformDerivativesNative out;
        out.ncomp = static_cast<int>(ncomp);
        out.nformula = static_cast<int>(nformula);
        out.nvar = static_cast<int>(nvar);
        out.aq_formula_dvar.assign(nformula * nvar, 0.0);
        out.org_formula_dvar.assign(nformula * nvar, 0.0);
        out.aq_formula_hessian_row_major.assign(nformula * nvar * nvar, 0.0);
        out.org_formula_hessian_row_major.assign(nformula * nvar * nvar, 0.0);

        ScalarRouteDerivativeState beta_formula_state = make_constant_scalar_state(beta_formula, nvar);
        beta_formula_state.gradient[0] = beta_formula * (1.0 - beta_formula);
        beta_formula_state.hessian_row_major[0] =
            beta_formula_state.gradient[0] * (1.0 - 2.0 * beta_formula);

        std::vector<ScalarRouteDerivativeState> org_formula_states;
        org_formula_states.reserve(nformula);
        for (std::size_t component = 0; component < nformula; ++component) {
            ScalarRouteDerivativeState state = make_constant_scalar_state(org_formula[component], nvar);
            for (std::size_t var = 1; var < nformula; ++var) {
                const std::size_t alpha = var - 1;
                const double delta_component_alpha = component == alpha ? 1.0 : 0.0;
                state.gradient[var] =
                    org_formula[component] * (delta_component_alpha - org_formula[alpha]);
            }
            for (std::size_t var_i = 1; var_i < nformula; ++var_i) {
                const std::size_t alpha = var_i - 1;
                const double delta_component_alpha = component == alpha ? 1.0 : 0.0;
                for (std::size_t var_j = 1; var_j < nformula; ++var_j) {
                    const std::size_t beta = var_j - 1;
                    const double delta_component_beta = component == beta ? 1.0 : 0.0;
                    const double delta_alpha_beta = alpha == beta ? 1.0 : 0.0;
                    state.hessian_row_major[var_i * nvar + var_j] =
                        org_formula[component] * (
                            (delta_component_beta - org_formula[beta])
                            * (delta_component_alpha - org_formula[alpha])
                            - org_formula[alpha] * (delta_alpha_beta - org_formula[beta])
                        );
                }
            }
            org_formula_states.push_back(std::move(state));
        }

        ScalarRouteDerivativeState one_minus_beta = subtract_scalar_states(
            make_constant_scalar_state(1.0, nvar),
            beta_formula_state
        );
        std::vector<ScalarRouteDerivativeState> aq_formula_states;
        aq_formula_states.reserve(nformula);
        for (std::size_t component = 0; component < nformula; ++component) {
            ScalarRouteDerivativeState numerator = subtract_scalar_states(
                make_constant_scalar_state(basis_.formula_feed[component], nvar),
                multiply_scalar_states(beta_formula_state, org_formula_states[component])
            );
            aq_formula_states.push_back(divide_scalar_states(numerator, one_minus_beta));
        }

        for (std::size_t component = 0; component < nformula; ++component) {
            for (std::size_t first = 0; first < nvar; ++first) {
                out.org_formula_dvar[component * nvar + first] = org_formula_states[component].gradient[first];
                out.aq_formula_dvar[component * nvar + first] = aq_formula_states[component].gradient[first];
                for (std::size_t second = 0; second < nvar; ++second) {
                    out.org_formula_hessian_row_major[component * nvar * nvar + first * nvar + second] =
                        org_formula_states[component].hessian_row_major[first * nvar + second];
                    out.aq_formula_hessian_row_major[component * nvar * nvar + first * nvar + second] =
                        aq_formula_states[component].hessian_row_major[first * nvar + second];
                }
            }
        }
        return out;
    }

    void select_separation_component(const ElectrolyteCandidateNative& initial) {
        double best = 0.0;
        separation_species_index_ = basis_.neutral_indices.size() >= 2 ? basis_.neutral_indices[1] : 0;
        separation_formula_index_ = basis_.neutral_indices.size() >= 2 ? 1 : 0;
        separation_sign_ = 1.0;
        for (std::size_t pos = 0; pos < basis_.neutral_indices.size(); ++pos) {
            const std::size_t species = static_cast<std::size_t>(basis_.neutral_indices[pos]);
            const double diff = initial.org_comp[species] - initial.aq_comp[species];
            if (std::abs(diff) > best) {
                best = std::abs(diff);
                separation_species_index_ = static_cast<int>(species);
                separation_formula_index_ = static_cast<int>(pos);
                separation_sign_ = diff >= 0.0 ? 1.0 : -1.0;
            }
        }
    }

    std::shared_ptr<ePCSAFTMixtureNative> mixture_;
    double temperature_ = 0.0;
    double target_pressure_ = 0.0;
    EquilibriumOptionsNative options_;
    std::vector<std::string> species_;
    std::vector<double> feed_;
    ElectrolyteBasisNative basis_;
    std::vector<double> initial_variables_;
    std::size_t residual_constraint_count_ = 0;
    double gibbs_feed_ = 0.0;
    double minimum_phase_distance_ = kLiquidLleMinimumPhaseDistance;
    double minimum_density_ = 1.0e-12;
    double maximum_density_ = 1.0e8;
    double density_bound_factor_ = 20.0;
    double aqueous_density_lower_bound_ = 1.0e-12;
    double aqueous_density_upper_bound_ = 1.0e8;
    double organic_density_lower_bound_ = 1.0e-12;
    double organic_density_upper_bound_ = 1.0e8;
    double aqueous_pressure_constraint_scale_ = 1.0;
    double organic_pressure_constraint_scale_ = 1.0;
    int separation_species_index_ = 0;
    int separation_formula_index_ = 0;
    double separation_sign_ = 1.0;
};

equilibrium_nlp::NeutralTwoPhaseEosNlpContract liquid_root_contract_from_problem(
    const LiquidRootElectrolyteLleProblem& problem
) {
    equilibrium_nlp::validate_nlp_problem_shape(problem);
    const std::vector<double> initial = problem.initial_point();
    const equilibrium_nlp::NlpBounds bounds = problem.bounds();
    const equilibrium_nlp::NlpJacobianStructure structure = problem.jacobian_structure();

    equilibrium_nlp::NeutralTwoPhaseEosNlpContract out;
    out.problem_name = problem.name();
    out.derivative_backend = "cppad_explicit_density";
    out.variable_model = problem.basis().variable_model + "_explicit_density";
    out.density_backend = "explicit_log_density_pressure_constraint";
    out.phase_count = 2;
    out.species_count = static_cast<int>(problem.feed().size());
    out.variable_count = problem.variable_count();
    out.constraint_count = problem.constraint_count();
    out.jacobian_nonzero_count = problem.jacobian_nonzero_count();
    out.initial_point = initial;
    out.variable_lower_bounds = bounds.variable_lower;
    out.variable_upper_bounds = bounds.variable_upper;
    out.constraint_lower_bounds = bounds.constraint_lower;
    out.constraint_upper_bounds = bounds.constraint_upper;
    out.objective_at_initial = problem.objective(initial);
    out.gradient_at_initial = problem.objective_gradient(initial);
    out.constraints_at_initial = problem.constraints(initial);
    out.jacobian_rows = structure.rows;
    out.jacobian_cols = structure.cols;
    out.jacobian_values_at_initial = problem.jacobian_values(initial);
    return out;
}

equilibrium_nlp::NeutralTwoPhaseEosPostsolve liquid_root_postsolve_from_candidate(
    const LiquidRootElectrolyteLleProblem& problem,
    const ElectrolyteCandidateNative& candidate,
    double material_tolerance,
    double charge_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    equilibrium_nlp::NeutralTwoPhaseEosPostsolve out;
    out.derivative_backend = "cppad_explicit_density";
    out.phase_count = 2;
    out.species_count = static_cast<int>(candidate.aq_comp.size());
    out.material_balance_norm = candidate.material_error;
    out.pressure_consistency_norm = problem.pressure_consistency_norm(candidate);
    out.chemical_potential_consistency_norm = candidate.solver_residual_norm;
    out.ln_fugacity_consistency_norm = candidate.solver_residual_norm;
    out.charge_balance_norm = candidate.charge_balance_error;
    out.phase_distance = candidate.phase_distance_value;
    out.objective = candidate.objective;
    out.gibbs_feed = candidate.gibbs_feed;
    out.gibbs_split = candidate.gibbs_split;
    out.gibbs_delta = candidate.gibbs_delta;
    out.minimum_phase_fraction = problem.minimum_phase_fraction();
    out.density_backend = "explicit_log_density_pressure_constraint";
    out.constraints = electrolyte_lle_residual_vector_from_candidate(candidate);
    out.constraints.push_back(problem.pressure_consistency_norm(candidate));
    out.phase_amount_totals = {1.0 - candidate.beta_org, candidate.beta_org};
    out.phase_compositions = {candidate.aq_comp, candidate.org_comp};
    out.phase_densities = {candidate.aq_state.density, candidate.org_state.density};
    out.phase_ln_fugacity_coefficients = {candidate.aq_state.ln_phi, candidate.org_state.ln_phi};
    out.phase_volumes = {
        out.phase_amount_totals[0] / candidate.aq_state.density,
        out.phase_amount_totals[1] / candidate.org_state.density,
    };
    const double effective_chemical_tolerance = std::max(
        chemical_potential_tolerance,
        2.0 * std::sqrt(chemical_potential_tolerance)
    );
    const double effective_phase_distance = std::max(
        phase_distance_tolerance,
        problem.minimum_phase_distance()
    );
    const bool finite_liquid_densities = std::isfinite(candidate.aq_state.density)
        && std::isfinite(candidate.org_state.density)
        && candidate.aq_state.density > 0.0
        && candidate.org_state.density > 0.0;
    const bool retained_phase_split = out.phase_amount_totals[0] >= problem.minimum_phase_fraction()
        && out.phase_amount_totals[1] >= problem.minimum_phase_fraction();
    out.accepted = out.material_balance_norm <= material_tolerance
        && out.pressure_consistency_norm <= problem.pressure_consistency_tolerance()
        && out.charge_balance_norm <= charge_tolerance
        && out.ln_fugacity_consistency_norm <= effective_chemical_tolerance
        && out.chemical_potential_consistency_norm <= effective_chemical_tolerance
        && out.phase_distance >= effective_phase_distance
        && finite_liquid_densities
        && retained_phase_split;
    if (out.accepted) {
        out.rejection_reason = "accepted";
    } else if (out.material_balance_norm > material_tolerance) {
        out.rejection_reason = "material_balance";
    } else if (out.pressure_consistency_norm > problem.pressure_consistency_tolerance()) {
        out.rejection_reason = "pressure_consistency";
    } else if (out.charge_balance_norm > charge_tolerance) {
        out.rejection_reason = "charge_balance";
    } else if (out.ln_fugacity_consistency_norm > effective_chemical_tolerance) {
        out.rejection_reason = "ln_fugacity_consistency";
    } else if (out.chemical_potential_consistency_norm > effective_chemical_tolerance) {
        out.rejection_reason = "chemical_potential_consistency";
    } else if (!finite_liquid_densities) {
        out.rejection_reason = "liquid_density";
    } else if (!retained_phase_split) {
        out.rejection_reason = "phase_fraction";
    } else {
        out.rejection_reason = "phase_distance";
    }
    return out;
}

std::vector<std::vector<double>> phase_amounts_from_liquid_root_candidate(
    const ElectrolyteCandidateNative& candidate
) {
    std::vector<std::vector<double>> out(2);
    out[0].reserve(candidate.aq_comp.size());
    out[1].reserve(candidate.org_comp.size());
    for (double value : candidate.aq_comp) {
        out[0].push_back((1.0 - candidate.beta_org) * value);
    }
    for (double value : candidate.org_comp) {
        out[1].push_back(candidate.beta_org * value);
    }
    return out;
}

std::vector<double> explicit_density_seed_from_electrolyte_formula_variables(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const std::vector<double>& formula_variables,
    const EquilibriumOptionsNative& options
) {
    ElectrolyteCandidateNative rooted = build_predictive_electrolyte_formula_candidate(
        feed,
        basis,
        formula_variables,
        options
    );
    rooted.aq_state = phase_state(mixture, t, p, rooted.aq_comp, "liq", "aq_density_seed");
    rooted.org_state = phase_state(mixture, t, p, rooted.org_comp, "liq", "org_density_seed");
    return append_electrolyte_phase_density_variables(
        formula_variables,
        rooted.aq_state.density,
        rooted.org_state.density
    );
}

int neutral_route_quality(const equilibrium_nlp::NeutralTwoPhaseEosRouteResult& result) {
    if (result.accepted) {
        return 3;
    }
    if (result.solver_accepted) {
        return 2;
    }
    if (result.ran) {
        return 1;
    }
    return 0;
}

equilibrium_nlp::RouteSeedAttempt neutral_seed_attempt_from_result(
    const equilibrium_nlp::NeutralTwoPhaseEosRouteResult& result
) {
    equilibrium_nlp::RouteSeedAttempt out;
    out.seed_name = result.seed_name;
    out.status = result.status;
    out.solver_status = result.solver_status;
    out.application_status = result.application_status;
    out.solver_accepted = result.solver_accepted;
    out.accepted = result.accepted;
    out.iteration_count = result.iteration_count;
    out.objective = result.objective;
    out.phase_distance = result.postsolve.phase_distance;
    out.material_balance_norm = result.postsolve.material_balance_norm;
    out.charge_balance_norm = result.postsolve.charge_balance_norm;
    out.pressure_consistency_norm = result.postsolve.pressure_consistency_norm;
    out.chemical_potential_consistency_norm = result.postsolve.chemical_potential_consistency_norm;
    return out;
}

equilibrium_nlp::NeutralTwoPhaseEosNlpContract evaluate_electrolyte_lle_liquid_root_nlp_contract_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& species,
    double phase_distance_tolerance
) {
    const LiquidRootElectrolyteLleProblem problem(
        mixture,
        t,
        p,
        feed,
        options,
        species,
        phase_distance_tolerance
    );
    return liquid_root_contract_from_problem(problem);
}

equilibrium_nlp::NeutralTwoPhaseEosRouteResult solve_electrolyte_lle_liquid_root_route_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& equilibrium_options,
    const std::vector<std::string>& species,
    const equilibrium_nlp::IpoptSolveOptions& solve_options,
    double material_tolerance,
    double charge_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
) {
    const equilibrium_nlp::IpoptAdapterInfo adapter = equilibrium_nlp::native_ipopt_adapter_info();
    equilibrium_nlp::NeutralTwoPhaseEosRouteResult best;
    best.compiled = adapter.compiled;
    best.adapter_available = adapter.adapter_available;
    best.adapter_kind = adapter.adapter_kind;
    best.exact_gradient_required = adapter.exact_gradient_required;
    best.exact_jacobian_required = adapter.exact_jacobian_required;
    best.problem_name = "electrolyte_lle_eos";
    best.derivative_backend = "cppad_explicit_density";
    best.postsolve.derivative_backend = "cppad_explicit_density";
    if (!adapter.compiled) {
        best.status = "ipopt_dependency_required";
        return best;
    }

    const LiquidRootElectrolyteLleProblem problem(
        mixture,
        t,
        p,
        feed,
        equilibrium_options,
        species,
        phase_distance_tolerance
    );
    const std::vector<NamedInitialVariables> seeds = {
        {
            "canonical_formula_shift",
            problem.explicit_density_seed_from_formula_variables(
                build_liquid_root_electrolyte_lle_initial_variables(problem.basis(), 1.0)
            )
        },
        {
            "mirrored_formula_shift",
            problem.explicit_density_seed_from_formula_variables(
                build_liquid_root_electrolyte_lle_initial_variables(problem.basis(), -1.0)
            )
        },
    };
    bool have_best = false;
    std::vector<equilibrium_nlp::RouteSeedAttempt> attempts;
    attempts.reserve(seeds.size() + (solve_options.initial_variables.empty() ? 0 : 1));

    auto run_attempt = [&](
        const std::string& seed_name,
        const equilibrium_nlp::IpoptSolveOptions& attempt_options
    ) {
        const equilibrium_nlp::IpoptSolveResult solve = equilibrium_nlp::solve_ipopt_nlp(problem, attempt_options);
        equilibrium_nlp::NeutralTwoPhaseEosRouteResult result;
        result.compiled = adapter.compiled;
        result.adapter_available = adapter.adapter_available;
        result.adapter_kind = adapter.adapter_kind;
        result.exact_gradient_required = adapter.exact_gradient_required;
        result.exact_jacobian_required = adapter.exact_jacobian_required;
        result.problem_name = "electrolyte_lle_eos";
        result.derivative_backend = "cppad_explicit_density";
        result.postsolve.derivative_backend = "cppad_explicit_density";
        result.initial_point_strategy = "deterministic_seed_sweep";
        result.seed_name = seed_name;
        result.ran = solve.solver_ran;
        result.solver_accepted = solve.accepted;
        result.accepted = solve.accepted;
        result.solver_status = solve.solver_status;
        result.application_status = solve.application_status;
        equilibrium_nlp::apply_ipopt_solve_metadata(result, solve);
        const auto last_exception = solve.diagnostics_string.find("last_callback_exception");
        if (last_exception != solve.diagnostics_string.end()) {
            result.last_callback_exception = last_exception->second;
        }
        result.objective = solve.objective;
        result.variables = solve.variables;
        result.constraints = solve.constraints;
        if (!solve.accepted) {
            result.status = "solver_rejected";
            attempts.push_back(neutral_seed_attempt_from_result(result));
            if (!have_best || neutral_route_quality(result) > neutral_route_quality(best)) {
                best = result;
                have_best = true;
            }
            return result;
        }

        const ElectrolyteCandidateNative candidate = problem.candidate(solve.variables);
        result.phase_amounts = phase_amounts_from_liquid_root_candidate(candidate);
        result.postsolve = liquid_root_postsolve_from_candidate(
            problem,
            candidate,
            material_tolerance,
            charge_tolerance,
            chemical_potential_tolerance,
            phase_distance_tolerance
        );
        result.phase_volumes = result.postsolve.phase_volumes;
        result.accepted = result.postsolve.accepted;
        result.status = result.accepted ? "accepted" : "postsolve_rejected";
        attempts.push_back(neutral_seed_attempt_from_result(result));
        if (!have_best || neutral_route_quality(result) > neutral_route_quality(best)) {
            best = result;
            have_best = true;
        }
        return result;
    };

    if (!solve_options.initial_variables.empty()) {
        equilibrium_nlp::IpoptSolveOptions continuation_options = solve_options;
        if (continuation_options.initial_variables.size() == problem.basis().formula_feed.size()) {
            continuation_options.initial_variables =
                problem.explicit_density_seed_from_formula_variables(continuation_options.initial_variables);
        }
        const equilibrium_nlp::NeutralTwoPhaseEosRouteResult continuation =
            run_attempt("continuation_state", continuation_options);
        if (continuation.accepted) {
            best.seed_attempts = attempts;
            return best;
        }
    }

    for (const auto& seed : seeds) {
        equilibrium_nlp::IpoptSolveOptions attempt_options = solve_options;
        attempt_options.initial_variables = seed.variables;
        attempt_options.initial_bound_lower_multipliers.clear();
        attempt_options.initial_bound_upper_multipliers.clear();
        attempt_options.initial_constraint_multipliers.clear();
        const equilibrium_nlp::NeutralTwoPhaseEosRouteResult attempt = run_attempt(seed.seed_name, attempt_options);
        if (attempt.accepted) {
            break;
        }
    }

    best.initial_point_strategy = "deterministic_seed_sweep";
    best.seed_attempts = attempts;
    return best;
}
