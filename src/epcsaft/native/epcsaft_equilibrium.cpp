#include "epcsaft_equilibrium.h"
#include "equilibrium/equilibrium_helpers.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

PhaseStateCompositionSensitivityResult phase_state_ln_fugacity_composition_sensitivity_cpp(
    double t,
    double p,
    std::vector<double> x,
    int phase,
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

struct ElectrolyteTransformJacobianNative {
    std::vector<double> aq_comp_dvar;
    std::vector<double> org_comp_dvar;
    std::vector<double> beta_org_dvar;
    int ncomp = 0;
    int nvar = 0;
};

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

ElectrolyteTransformJacobianNative electrolyte_transform_jacobian(
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const ElectrolyteCandidateNative& candidate
) {
    const std::size_t ncomp = feed.size();
    const std::size_t nformula = basis.formula_feed.size();
    const std::size_t nvar = nformula;
    ElectrolyteTransformJacobianNative out;
    out.ncomp = static_cast<int>(ncomp);
    out.nvar = static_cast<int>(nvar);
    out.aq_comp_dvar.assign(ncomp * nvar, 0.0);
    out.org_comp_dvar.assign(ncomp * nvar, 0.0);
    out.beta_org_dvar.assign(nvar, 0.0);

    std::vector<double> dbeta_formula_dvar(nvar, 0.0);
    dbeta_formula_dvar[0] = candidate.beta_formula * (1.0 - candidate.beta_formula);

    std::vector<double> dorg_formula_dvar(nformula * nvar, 0.0);
    for (std::size_t var = 1; var < nvar; ++var) {
        const std::size_t logit_index = var - 1;
        for (std::size_t k = 0; k < nformula; ++k) {
            double delta = (k == logit_index) ? 1.0 : 0.0;
            dorg_formula_dvar[k * nvar + var] = candidate.org_formula[k] * (delta - candidate.org_formula[logit_index]);
        }
    }

    std::vector<double> daq_formula_dvar(nformula * nvar, 0.0);
    const double denom = 1.0 - candidate.beta_formula;
    for (std::size_t k = 0; k < nformula; ++k) {
        daq_formula_dvar[k * nvar] =
            dbeta_formula_dvar[0] * (basis.formula_feed[k] - candidate.org_formula[k]) / (denom * denom);
        for (std::size_t var = 1; var < nvar; ++var) {
            daq_formula_dvar[k * nvar + var] =
                -candidate.beta_formula / denom * dorg_formula_dvar[k * nvar + var];
        }
    }

    std::vector<double> aq_formula_to_explicit = formula_to_explicit_jacobian_row_major(candidate.aq_formula, basis, ncomp);
    std::vector<double> org_formula_to_explicit = formula_to_explicit_jacobian_row_major(candidate.org_formula, basis, ncomp);
    for (std::size_t i = 0; i < ncomp; ++i) {
        for (std::size_t var = 0; var < nvar; ++var) {
            for (std::size_t k = 0; k < nformula; ++k) {
                out.aq_comp_dvar[i * nvar + var] += aq_formula_to_explicit[i * nformula + k] * daq_formula_dvar[k * nvar + var];
                out.org_comp_dvar[i * nvar + var] += org_formula_to_explicit[i * nformula + k] * dorg_formula_dvar[k * nvar + var];
            }
        }
    }

    std::vector<double> column_sums = formula_expansion_column_sums(basis, ncomp);
    double aq_scale = 0.0;
    double org_scale = 0.0;
    for (std::size_t k = 0; k < nformula; ++k) {
        aq_scale += column_sums[k] * candidate.aq_formula[k];
        org_scale += column_sums[k] * candidate.org_formula[k];
    }
    double beta_denominator = (1.0 - candidate.beta_formula) * aq_scale + candidate.beta_formula * org_scale;
    if (!(beta_denominator > 0.0) || !std::isfinite(beta_denominator)) {
        throw SolutionError("Electrolyte transformed-variable Jacobian produced an invalid phase-fraction denominator.");
    }
    double beta_numerator = candidate.beta_formula * org_scale;
    for (std::size_t var = 0; var < nvar; ++var) {
        double d_aq_scale = 0.0;
        double d_org_scale = 0.0;
        for (std::size_t k = 0; k < nformula; ++k) {
            d_aq_scale += column_sums[k] * daq_formula_dvar[k * nvar + var];
            d_org_scale += column_sums[k] * dorg_formula_dvar[k * nvar + var];
        }
        double d_beta_formula = dbeta_formula_dvar[var];
        double d_num = d_beta_formula * org_scale + candidate.beta_formula * d_org_scale;
        double d_den = d_beta_formula * (org_scale - aq_scale)
            + (1.0 - candidate.beta_formula) * d_aq_scale
            + candidate.beta_formula * d_org_scale;
        out.beta_org_dvar[var] = (d_num * beta_denominator - beta_numerator * d_den)
            / (beta_denominator * beta_denominator);
    }
    return out;
}

double phase_log_fugacity_derivative_for_species(
    const std::vector<double>& composition,
    const std::vector<double>& dcomp_dvar,
    const std::vector<double>& lnphi_jacobian,
    std::size_t species_index,
    std::size_t var,
    std::size_t ncomp,
    std::size_t nvar
) {
    double value = dcomp_dvar[species_index * nvar + var] / composition[species_index];
    for (std::size_t j = 0; j < ncomp; ++j) {
        value += lnphi_jacobian[species_index * ncomp + j] * dcomp_dvar[j * nvar + var];
    }
    return value;
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
    const std::size_t nvar = basis.formula_feed.size();
    ElectrolyteTransformJacobianNative transform = electrolyte_transform_jacobian(feed, basis, candidate);
    PhaseStateCompositionSensitivityResult aq_sensitivity =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, candidate.aq_comp, phase_token_to_int("liq"), mixture->args());
    PhaseStateCompositionSensitivityResult org_sensitivity =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, candidate.org_comp, phase_token_to_int("liq"), mixture->args());
    if (!aq_sensitivity.supported || !org_sensitivity.supported) {
        throw SolutionError("Electrolyte residual Jacobian requires supported phase-state fugacity composition sensitivities.");
    }
    const std::size_t phase_rows = basis.neutral_indices.size() + basis.salt_pairs.size();
    const std::size_t rows = phase_rows + ncomp;
    std::vector<double> jacobian(rows * nvar, 0.0);
    std::size_t row = 0;
    auto add_species_contribution = [&](std::size_t output_row, std::size_t species_index, double coefficient) {
        for (std::size_t var = 0; var < nvar; ++var) {
            double org_value = phase_log_fugacity_derivative_for_species(
                candidate.org_comp,
                transform.org_comp_dvar,
                org_sensitivity.jacobian_row_major,
                species_index,
                var,
                ncomp,
                nvar
            );
            double aq_value = phase_log_fugacity_derivative_for_species(
                candidate.aq_comp,
                transform.aq_comp_dvar,
                aq_sensitivity.jacobian_row_major,
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
    ElectrolyteCandidateNative out;
    unpack_predictive_electrolyte_variables(variables, basis.formula_feed.size(), out.beta_formula, out.org_formula);
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
    out.aq_state = phase_state(mixture, t, p, out.aq_comp, "liq", "aq");
    out.org_state = phase_state(mixture, t, p, out.org_comp, "liq", "org");
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
        eval_variables = pack_predictive_electrolyte_variables(beta_formula, org_formula);
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
    out.variable_model = basis.variable_model;
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
    out.diagnostics_string["jacobian_backend"] = "cppad_implicit";
    out.diagnostics_string["derivative_backend"] = "cppad_implicit";
    out.diagnostics_string["jacobian_scope"] = "transformed_variables_phase_state_implicit_density";
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
    const ElectrolyteBasisNative& basis
) {
    const double beta_formula = 0.635;
    std::vector<double> org_formula = basis.formula_feed;
    if (basis.formula_feed.size() >= 3 && basis.neutral_indices.size() >= 2) {
        org_formula[0] *= 0.667;
        org_formula[1] *= 1.548;
        for (std::size_t pos = basis.neutral_indices.size(); pos < org_formula.size(); ++pos) {
            org_formula[pos] *= 0.106;
        }
    } else if (basis.formula_feed.size() >= 2) {
        org_formula[0] *= 0.75;
        org_formula[1] *= 1.25;
    }
    org_formula = cap_formula_seed_to_material_feasible_region(org_formula, basis.formula_feed, beta_formula);
    return pack_predictive_electrolyte_variables(beta_formula, org_formula);
}

std::vector<double> org_formula_jacobian_row_major(
    const std::vector<double>& org_formula,
    std::size_t nvar
) {
    const std::size_t nformula = org_formula.size();
    std::vector<double> jacobian(nformula * nvar, 0.0);
    for (std::size_t component = 0; component < nformula; ++component) {
        for (std::size_t var = 1; var < nvar; ++var) {
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
    unpack_predictive_electrolyte_variables(variables, basis.formula_feed.size(), beta_formula, org_formula);
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
    unpack_predictive_electrolyte_variables(variables, basis.formula_feed.size(), beta_formula, org_formula);
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
    const std::size_t cols = basis.formula_feed.size();
    std::vector<double> gradient(cols, 0.0);
    for (std::size_t row = 0; row < residual.size(); ++row) {
        for (std::size_t col = 0; col < cols; ++col) {
            gradient[col] += jacobian[row * cols + col] * residual[row];
        }
    }
    return gradient;
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
        initial_variables_ = build_liquid_root_electrolyte_lle_initial_variables(basis_);
        const ElectrolyteCandidateNative initial = candidate(initial_variables_);
        residual_constraint_count_ = initial.residual.size();
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
            residual_constraint_count_ + 1 + basis_.formula_feed.size()
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
        out.constraint_lower.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(constraint_count()), 0.0);
        const std::size_t phase_distance_row = residual_constraint_count_;
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
        } catch (const std::exception&) {
            out.insert(out.end(), residual_constraint_count_, 1.0e3);
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
        } catch (const std::exception&) {
        }
        double beta_formula = 0.5;
        std::vector<double> org_formula;
        unpack_predictive_electrolyte_variables(variables, basis_.formula_feed.size(), beta_formula, org_formula);
        const std::vector<double> org_jacobian = org_formula_jacobian_row_major(org_formula, nvar);
        const std::vector<double> aq_jacobian = aqueous_formula_feasibility_jacobian_row_major(variables, basis_);
        const std::size_t formula = static_cast<std::size_t>(separation_formula_index_);
        const std::size_t phase_distance_row = residual_constraint_count_;
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

    equilibrium_nlp::NlpScaling scaling() const override {
        equilibrium_nlp::NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(initial_variables_.size(), 1.0);
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0);
        return out;
    }

    std::map<std::string, std::string> diagnostics() const override {
        return {
            {"derivative_backend", "cppad_implicit"},
            {"density_backend", "liquid_pressure_root"},
            {"variable_model", basis_.variable_model},
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
        unpack_predictive_electrolyte_variables(variables, basis_.formula_feed.size(), beta_formula, org_formula);
        const std::vector<double> aq_formula = aqueous_formula_feasibility_values(variables, basis_);
        const std::size_t formula = static_cast<std::size_t>(separation_formula_index_);
        return separation_sign_ * (org_formula[formula] - aq_formula[formula]);
    }

private:
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
    out.derivative_backend = "cppad_implicit";
    out.variable_model = problem.basis().variable_model;
    out.density_backend = "liquid_pressure_root";
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
    out.derivative_backend = "cppad_implicit";
    out.phase_count = 2;
    out.species_count = static_cast<int>(candidate.aq_comp.size());
    out.material_balance_norm = candidate.material_error;
    out.pressure_consistency_norm = 0.0;
    out.chemical_potential_consistency_norm = candidate.solver_residual_norm;
    out.ln_fugacity_consistency_norm = candidate.solver_residual_norm;
    out.charge_balance_norm = candidate.charge_balance_error;
    out.phase_distance = candidate.phase_distance_value;
    out.objective = candidate.objective;
    out.gibbs_feed = candidate.gibbs_feed;
    out.gibbs_split = candidate.gibbs_split;
    out.gibbs_delta = candidate.gibbs_delta;
    out.minimum_phase_fraction = problem.minimum_phase_fraction();
    out.density_backend = "liquid_pressure_root";
    out.constraints = electrolyte_lle_residual_vector_from_candidate(candidate);
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
    equilibrium_nlp::NeutralTwoPhaseEosRouteResult out;
    out.compiled = adapter.compiled;
    out.adapter_available = adapter.adapter_available;
    out.adapter_kind = adapter.adapter_kind;
    out.exact_gradient_required = adapter.exact_gradient_required;
    out.exact_jacobian_required = adapter.exact_jacobian_required;
    out.problem_name = "electrolyte_lle_eos";
    out.derivative_backend = "cppad_implicit";
    out.postsolve.derivative_backend = "cppad_implicit";
    if (!adapter.compiled) {
        out.status = "ipopt_dependency_required";
        return out;
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
    const equilibrium_nlp::IpoptSolveResult solve = equilibrium_nlp::solve_ipopt_nlp(problem, solve_options);
    out.ran = solve.solver_ran;
    out.solver_accepted = solve.accepted;
    out.accepted = solve.accepted;
    out.solver_status = solve.solver_status;
    out.application_status = solve.application_status;
    const auto last_exception = solve.diagnostics_string.find("last_callback_exception");
    if (last_exception != solve.diagnostics_string.end()) {
        out.last_callback_exception = last_exception->second;
    }
    out.objective = solve.objective;
    out.variables = solve.variables;
    out.constraints = solve.constraints;
    if (!solve.accepted) {
        out.status = "solver_rejected";
        return out;
    }

    const ElectrolyteCandidateNative candidate = problem.candidate(solve.variables);
    out.phase_amounts = phase_amounts_from_liquid_root_candidate(candidate);
    out.postsolve = liquid_root_postsolve_from_candidate(
        problem,
        candidate,
        material_tolerance,
        charge_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
    out.phase_volumes = out.postsolve.phase_volumes;
    out.accepted = out.postsolve.accepted;
    out.status = out.accepted ? "accepted" : "postsolve_rejected";
    return out;
}
