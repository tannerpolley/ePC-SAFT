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
    out.objective = l2_norm(out.residual);
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
