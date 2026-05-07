#include "epcsaft_equilibrium.h"
#include "equilibrium/equilibrium_helpers.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <sstream>

namespace {

struct PhaseStateNative {
    std::shared_ptr<ePCSAFTStateNative> state;
    std::vector<double> ln_phi;
    double density = 0.0;
};

struct LLESeedNative {
    std::string seed_name;
    double beta = 0.5;
    std::vector<double> comp1;
    std::vector<double> comp2;
};

struct LLECandidateNative {
    double beta = 0.5;
    std::vector<double> comp1;
    std::vector<double> comp2;
    PhaseStateNative state1;
    PhaseStateNative state2;
    std::vector<double> fugacity_residual;
    std::vector<double> material_residual;
    std::vector<double> residual;
    double fugacity_residual_norm = 0.0;
    double material_error = 0.0;
    int iteration = 0;
    double objective = std::numeric_limits<double>::infinity();
};

struct LLEAttemptNative {
    std::string status;
    std::string seed_name;
    int attempt_count = 0;
    std::string message;
    LLECandidateNative candidate;
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
    int basis_rank = 0;
    std::string variable_model = "ascani_transformed_salt_pairs";
};

struct ElectrolyteTpdSearchNative {
    std::vector<StabilityTrialNative> trials;
    int multistart_count = 0;
    int polish_iterations = 0;
};

bool rachford_rice_beta(const std::vector<double>& feed, const std::vector<double>& k_values, double& beta, std::string& message);
std::pair<std::vector<double>, std::vector<double>> phase_compositions(
    const std::vector<double>& feed,
    const std::vector<double>& k_values,
    double beta,
    double min_composition
);

using namespace epcsaft::native::equilibrium;

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

void append_last_density_failure(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    std::vector<DensitySolveDiagnostics>* failures
) {
    if (failures == nullptr) {
        return;
    }
    DensitySolveDiagnostics diagnostics = mixture->last_density_diagnostics();
    if (diagnostics.validity_gate != "failed") {
        return;
    }
    failures->push_back(diagnostics);
}

EquilibriumPhaseNative phase_from_state(
    const std::string& label,
    const std::vector<double>& composition,
    double phase_fraction,
    const PhaseStateNative& state,
    bool include_phase_diagnostics
) {
    EquilibriumPhaseNative phase;
    phase.label = label;
    phase.composition = composition;
    phase.density = state.density;
    phase.temperature = state.state->temperature();
    phase.pressure = state.state->pressure();
    phase.phase_fraction = phase_fraction;
    phase.ln_fugacity_coefficient = state.ln_phi;
    if (include_phase_diagnostics) {
        phase.diagnostics_string["phase"] = label;
        phase.diagnostics_double["density"] = state.density;
        phase.diagnostics_string["fugacity_coefficient_terms"] = "native";
    }
    return phase;
}

double tpd_value(
    const std::vector<double>& composition,
    const std::vector<double>& trial_ln_phi,
    const std::vector<double>& feed,
    const std::vector<double>& parent_ln_phi
) {
    double value = 0.0;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        value += composition[i] * (std::log(composition[i]) + trial_ln_phi[i] - std::log(feed[i]) - parent_ln_phi[i]);
    }
    return value;
}

StabilityTrialNative solve_tpd_trial(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const std::vector<double>& parent_ln_phi,
    const std::string& parent_phase,
    const std::string& trial_phase,
    const std::string& seed_name,
    const std::vector<double>& seed_composition,
    const EquilibriumOptionsNative& options,
    double threshold
) {
    std::vector<double> composition = clip_normalize(seed_composition, options.min_composition);
    double best_tpd = std::numeric_limits<double>::infinity();
    std::vector<double> best_composition = composition;
    int best_iteration = 0;
    bool converged = false;
    double max_delta = std::numeric_limits<double>::infinity();
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        PhaseStateNative trial = phase_state(mixture, t, p, composition, trial_phase);
        double candidate_tpd = tpd_value(composition, trial.ln_phi, feed, parent_ln_phi);
        if (candidate_tpd < best_tpd) {
            best_tpd = candidate_tpd;
            best_composition = composition;
            best_iteration = iteration;
        }
        std::vector<double> target_weights(feed.size(), 0.0);
        for (std::size_t i = 0; i < feed.size(); ++i) {
            target_weights[i] = std::log(feed[i]) + parent_ln_phi[i] - trial.ln_phi[i];
        }
        std::vector<double> target = composition_from_log_weights(target_weights, options.min_composition);
        std::vector<double> next(composition.size(), 0.0);
        for (std::size_t i = 0; i < composition.size(); ++i) {
            next[i] = (1.0 - options.damping) * composition[i] + options.damping * target[i];
        }
        next = clip_normalize(next, options.min_composition);
        max_delta = phase_distance(next, composition);
        if (max_delta <= std::max(options.tolerance, 1.0e-10)) {
            converged = true;
            break;
        }
        composition = next;
    }
    StabilityTrialNative out;
    out.parent_phase = parent_phase;
    out.trial_phase = trial_phase;
    out.seed_name = seed_name;
    out.composition = best_composition;
    out.tpd = best_tpd;
    out.iterations = best_iteration;
    out.converged = converged;
    out.unstable = best_tpd < -threshold;
    out.diagnostics_double["tpd_threshold"] = threshold;
    out.diagnostics_double["final_max_composition_delta"] = max_delta;
    out.diagnostics_string["message"] = out.unstable ? "negative TPD trial" : "non-negative TPD trial";
    return out;
}

std::vector<std::pair<std::string, std::vector<double>>> tpd_seeds(const std::vector<double>& feed, const EquilibriumOptionsNative& options) {
    std::vector<std::pair<std::string, std::vector<double>>> seeds;
    seeds.push_back({"feed", feed});
    for (std::size_t i = 0; i < feed.size(); ++i) {
        seeds.push_back({"component_" + std::to_string(i) + "_rich", component_rich_composition(feed.size(), i, options.min_composition)});
    }
    return seeds;
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

std::vector<int> species_indices_from_labels(
    const std::vector<std::string>& species,
    const std::vector<std::string>& labels,
    const std::string& field_name
) {
    if (species.empty()) {
        throw ValueError(field_name + " requires species labels from the Python mixture.");
    }
    std::vector<int> out;
    for (const std::string& label : labels) {
        auto it = std::find(species.begin(), species.end(), label);
        if (it == species.end()) {
            throw ValueError(field_name + " contains unknown species label: " + label);
        }
        int index = static_cast<int>(std::distance(species.begin(), it));
        if (std::find(out.begin(), out.end(), index) == out.end()) {
            out.push_back(index);
        }
    }
    return out;
}

std::vector<double> normalize_nonnegative_composition(
    const std::vector<double>& values,
    std::size_t ncomp,
    double min_composition,
    const std::string& field_name
) {
    if (values.size() != ncomp) {
        throw ValueError(field_name + " length must match mixture component count.");
    }
    double total = 0.0;
    std::vector<double> out(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i])) {
            throw ValueError(field_name + " must contain only finite values.");
        }
        if (values[i] < 0.0) {
            throw ValueError(field_name + " must be non-negative.");
        }
        out[i] = std::max(values[i], min_composition);
        total += out[i];
    }
    if (!std::isfinite(total) || total <= 0.0) {
        throw ValueError(field_name + " must have a positive finite sum.");
    }
    for (double& value : out) {
        value /= total;
    }
    return out;
}

std::vector<double> vapor_full_composition(
    const std::vector<int>& vapor_indices,
    const std::vector<double>& y_vap,
    std::size_t ncomp,
    double min_composition
) {
    std::vector<double> out(ncomp, min_composition);
    for (std::size_t pos = 0; pos < vapor_indices.size(); ++pos) {
        out[static_cast<std::size_t>(vapor_indices[pos])] = std::max(y_vap[pos], min_composition);
    }
    double total = std::accumulate(out.begin(), out.end(), 0.0);
    for (double& value : out) {
        value /= total;
    }
    return out;
}

std::vector<double> subset_double_vector(const std::vector<double>& values, const std::vector<int>& indices, std::size_t ncomp) {
    if (values.empty() || values.size() != ncomp) {
        return values;
    }
    std::vector<double> out;
    out.reserve(indices.size());
    for (int index : indices) {
        out.push_back(values[static_cast<std::size_t>(index)]);
    }
    return out;
}

std::vector<int> subset_int_vector(const std::vector<int>& values, const std::vector<int>& indices, std::size_t ncomp) {
    if (values.empty() || values.size() != ncomp) {
        return values;
    }
    std::vector<int> out;
    out.reserve(indices.size());
    for (int index : indices) {
        out.push_back(values[static_cast<std::size_t>(index)]);
    }
    return out;
}

std::vector<double> subset_double_matrix(const std::vector<double>& values, const std::vector<int>& indices, std::size_t ncomp) {
    if (values.empty() || values.size() != ncomp * ncomp) {
        return values;
    }
    std::vector<double> out(indices.size() * indices.size(), 0.0);
    for (std::size_t row = 0; row < indices.size(); ++row) {
        for (std::size_t col = 0; col < indices.size(); ++col) {
            out[row * indices.size() + col] = values[static_cast<std::size_t>(indices[row]) * ncomp + static_cast<std::size_t>(indices[col])];
        }
    }
    return out;
}

std::vector<int> subset_int_matrix(const std::vector<int>& values, const std::vector<int>& indices, std::size_t ncomp) {
    if (values.empty() || values.size() != ncomp * ncomp) {
        return values;
    }
    std::vector<int> out(indices.size() * indices.size(), 0);
    for (std::size_t row = 0; row < indices.size(); ++row) {
        for (std::size_t col = 0; col < indices.size(); ++col) {
            out[row * indices.size() + col] = values[static_cast<std::size_t>(indices[row]) * ncomp + static_cast<std::size_t>(indices[col])];
        }
    }
    return out;
}

std::shared_ptr<ePCSAFTMixtureNative> vapor_submixture(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<int>& vapor_indices
) {
    const add_args& args = mixture->args();
    std::size_t ncomp = mixture->ncomp();
    add_args sub = args;
    sub.m = subset_double_vector(args.m, vapor_indices, ncomp);
    sub.s = subset_double_vector(args.s, vapor_indices, ncomp);
    sub.e = subset_double_vector(args.e, vapor_indices, ncomp);
    sub.e_assoc = subset_double_vector(args.e_assoc, vapor_indices, ncomp);
    sub.vol_a = subset_double_vector(args.vol_a, vapor_indices, ncomp);
    sub.z = std::vector<double>(vapor_indices.size(), 0.0);
    sub.dielc = subset_double_vector(args.dielc, vapor_indices, ncomp);
    sub.mw = subset_double_vector(args.mw, vapor_indices, ncomp);
    sub.mixed_rel_perm_a = subset_double_vector(args.mixed_rel_perm_a, vapor_indices, ncomp);
    sub.mixed_rel_perm_b = subset_double_vector(args.mixed_rel_perm_b, vapor_indices, ncomp);
    sub.mixed_rel_perm_c = subset_double_vector(args.mixed_rel_perm_c, vapor_indices, ncomp);
    sub.mixed_rel_perm_mask = subset_int_vector(args.mixed_rel_perm_mask, vapor_indices, ncomp);
    sub.d_born = subset_double_vector(args.d_born, vapor_indices, ncomp);
    sub.f_solv = subset_double_vector(args.f_solv, vapor_indices, ncomp);
    sub.assoc_num = subset_int_vector(args.assoc_num, vapor_indices, ncomp);
    sub.k_ij = subset_double_matrix(args.k_ij, vapor_indices, ncomp);
    sub.l_ij = subset_double_matrix(args.l_ij, vapor_indices, ncomp);
    sub.assoc_matrix = subset_int_matrix(args.assoc_matrix, vapor_indices, ncomp);
    sub.k_hb = subset_double_matrix(args.k_hb, vapor_indices, ncomp);
    sub.mixed_rel_perm_water_index = -1;
    for (std::size_t pos = 0; pos < vapor_indices.size(); ++pos) {
        if (vapor_indices[pos] == args.mixed_rel_perm_water_index) {
            sub.mixed_rel_perm_water_index = static_cast<int>(pos);
            break;
        }
    }
    return std::make_shared<ePCSAFTMixtureNative>(sub);
}

struct ElectrolyteBubbleEvaluationNative {
    bool finite = false;
    double p = 0.0;
    double objective = std::numeric_limits<double>::infinity();
    double residual_norm = std::numeric_limits<double>::infinity();
    double charge_residual = 0.0;
    int vapor_iterations = 0;
    std::string message;
    std::vector<double> y_vap;
    std::vector<double> y_full;
    std::vector<double> partial_pressures;
    std::vector<double> fugacity_residual;
    PhaseStateNative liquid;
    PhaseStateNative vapor;
};

ElectrolyteBubbleEvaluationNative evaluate_electrolyte_bubble_pressure(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::shared_ptr<ePCSAFTMixtureNative>& vapor_mixture,
    double t,
    double p,
    const std::vector<double>& x_liq,
    const std::vector<int>& vapor_indices,
    const ElectrolyteBubbleOptionsNative& options,
    const std::vector<double>& y_seed
) {
    ElectrolyteBubbleEvaluationNative out;
    out.p = p;
    out.charge_residual = composition_charge(x_liq, mixture->args().z);
    if (!(p > 0.0) || !std::isfinite(p)) {
        out.message = "pressure is not positive finite";
        return out;
    }
    try {
        mixture->clear_runtime_caches();
        vapor_mixture->clear_runtime_caches();
        std::vector<double> y(vapor_indices.size(), 1.0 / static_cast<double>(vapor_indices.size()));
        if (y_seed.size() == vapor_indices.size()) {
            y = clip_normalize(y_seed, options.min_composition);
        }
        PhaseStateNative liquid = phase_state(mixture, t, p, x_liq, "liq", "electrolyte_bubble_liq");
        PhaseStateNative vapor;
        std::vector<double> raw(y.size(), 0.0);
        double sum_raw = 0.0;
        int iteration = 0;
        for (iteration = 1; iteration <= options.max_vapor_iterations; ++iteration) {
            vapor = phase_state(vapor_mixture, t, p, y, "vap", "electrolyte_bubble_vap");
            if (vapor.density > 0.5 * liquid.density) {
                out.message = "vapor density root is liquid-like";
                return out;
            }
            sum_raw = 0.0;
            for (std::size_t pos = 0; pos < vapor_indices.size(); ++pos) {
                std::size_t i = static_cast<std::size_t>(vapor_indices[pos]);
                double ln_k = liquid.ln_phi[i] - vapor.ln_phi[pos];
                raw[pos] = std::max(x_liq[i], options.min_composition) * std::exp(std::clamp(ln_k, -700.0, 700.0));
                sum_raw += raw[pos];
            }
            if (!std::isfinite(sum_raw) || sum_raw <= 0.0) {
                out.message = "bubble objective produced invalid K-values";
                return out;
            }
            double max_delta = 0.0;
            std::vector<double> y_next(y.size(), 0.0);
            for (std::size_t pos = 0; pos < y.size(); ++pos) {
                y_next[pos] = std::max(raw[pos] / sum_raw, options.min_composition);
                max_delta = std::max(max_delta, std::abs(y_next[pos] - y[pos]));
            }
            y_next = clip_normalize(y_next, options.min_composition);
            y = y_next;
            if (max_delta <= options.vapor_tolerance) {
                break;
            }
        }
        vapor = phase_state(vapor_mixture, t, p, y, "vap", "electrolyte_bubble_vap");
        if (vapor.density > 0.5 * liquid.density) {
            out.message = "vapor density root is liquid-like";
            return out;
        }
        sum_raw = 0.0;
        for (std::size_t pos = 0; pos < vapor_indices.size(); ++pos) {
            std::size_t i = static_cast<std::size_t>(vapor_indices[pos]);
            double ln_k = liquid.ln_phi[i] - vapor.ln_phi[pos];
            raw[pos] = std::max(x_liq[i], options.min_composition) * std::exp(std::clamp(ln_k, -700.0, 700.0));
            sum_raw += raw[pos];
        }
        out.y_vap = clip_normalize(raw, options.min_composition);
        out.y_full = out.y_vap;
        out.partial_pressures.resize(out.y_vap.size(), 0.0);
        out.fugacity_residual.resize(out.y_vap.size(), 0.0);
        for (std::size_t pos = 0; pos < vapor_indices.size(); ++pos) {
            std::size_t i = static_cast<std::size_t>(vapor_indices[pos]);
            out.partial_pressures[pos] = out.y_vap[pos] * p;
            out.fugacity_residual[pos] = std::log(std::max(out.y_vap[pos], options.min_composition))
                + vapor.ln_phi[pos]
                - std::log(std::max(x_liq[i], options.min_composition))
                - liquid.ln_phi[i];
        }
        out.finite = std::isfinite(sum_raw) && sum_raw > 0.0;
        out.objective = sum_raw - 1.0;
        out.residual_norm = max_abs(out.fugacity_residual);
        out.vapor_iterations = iteration;
        out.liquid = liquid;
        out.vapor = vapor;
        out.message = "finite";
        return out;
    } catch (const std::exception& exc) {
        out.message = exc.what();
        return out;
    }
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

StabilityTrialNative electrolyte_trial_from_formula(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const PhaseStateNative& parent,
    const ElectrolyteBasisNative& basis,
    const std::vector<double>& formula,
    const std::string& seed_name,
    int iterations,
    double threshold
) {
    auto expanded = formula_to_explicit(formula, basis, feed.size());
    std::vector<double> composition = expanded.first;
    PhaseStateNative trial_state = phase_state(mixture, t, p, composition, "liq");
    double tpd = tpd_value(composition, trial_state.ln_phi, feed, parent.ln_phi);
    StabilityTrialNative trial;
    trial.parent_phase = "liq";
    trial.trial_phase = "liq";
    trial.seed_name = seed_name;
    trial.composition = composition;
    trial.tpd = tpd;
    trial.iterations = iterations;
    trial.converged = true;
    trial.unstable = tpd < -threshold;
    trial.diagnostics_double["tpd_threshold"] = threshold;
    trial.diagnostics_double["tpd_objective_value"] = tpd;
    trial.diagnostics_string["stability_analysis"] = "electrolyte_tpd";
    trial.diagnostics_string["tpd_method"] = "native_tpd_global_search";
    trial.diagnostics_string["message"] = trial.unstable ? "negative transformed electrolyte TPD trial" : "non-negative transformed electrolyte TPD trial";
    return trial;
}

std::vector<std::pair<std::string, std::vector<double>>> electrolyte_formula_seeds(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const EquilibriumOptionsNative& options
) {
    std::vector<std::pair<std::string, std::vector<double>>> seeds;
    seeds.push_back({"formula_feed", basis.formula_feed});
    for (std::size_t i = 0; i < basis.formula_feed.size(); ++i) {
        seeds.push_back({"formula_component_" + std::to_string(i) + "_rich", component_rich_composition(basis.formula_feed.size(), i, options.min_composition)});
    }
    std::vector<double> k(feed.size(), 1.0);
    const auto& charges = mixture->args().z;
    int aq_index = basis.neutral_indices[0];
    int org_index = basis.neutral_indices.size() > 1 ? basis.neutral_indices[1] : basis.neutral_indices[0];
    k[static_cast<std::size_t>(aq_index)] = 0.08;
    k[static_cast<std::size_t>(org_index)] = 8.0;
    for (std::size_t i = 0; i < feed.size(); ++i) {
        if (std::abs(charges[i]) > 1.0e-12) {
            k[i] = 0.02;
        }
    }
    double beta = 0.0;
    std::string message;
    if (rachford_rice_beta(feed, k, beta, message)) {
        auto comps = phase_compositions(feed, k, beta, options.min_composition);
        seeds.push_back({"partition_seed", explicit_to_formula(comps.second, basis)});
    }
    return seeds;
}

std::vector<std::pair<std::string, std::vector<double>>> electrolyte_endpoint_phase_seeds(
    const ElectrolyteBasisNative& basis,
    const EquilibriumOptionsNative& options
) {
    std::vector<std::pair<std::string, std::vector<double>>> seeds;
    if (basis.neutral_indices.size() < 2 || basis.salt_pairs.empty()) {
        return seeds;
    }
    const std::size_t n_formula = basis.formula_feed.size();
    const std::size_t aq_index = 0;
    const std::size_t org_index = 1;
    const std::size_t salt_index = basis.neutral_indices.size();
    for (double salt_share : {0.001, 0.01, 0.05, 0.15}) {
        std::vector<double> org_formula(n_formula, options.min_composition);
        org_formula[org_index] = 0.80;
        org_formula[aq_index] = 0.20;
        org_formula[salt_index] = salt_share;
        seeds.push_back({"org_endpoint_salt_" + std::to_string(static_cast<int>(salt_share * 1000)), clip_normalize(org_formula, options.min_composition)});

        std::vector<double> aq_formula(n_formula, options.min_composition);
        aq_formula[aq_index] = 0.95;
        aq_formula[org_index] = 0.05;
        aq_formula[salt_index] = salt_share;
        seeds.push_back({"aq_endpoint_salt_" + std::to_string(static_cast<int>(salt_share * 1000)), clip_normalize(aq_formula, options.min_composition)});
    }
    return seeds;
}

EquilibriumOptionsNative precheck_options(const EquilibriumOptionsNative& options) {
    EquilibriumOptionsNative out = options;
    out.max_iterations = std::min(options.max_iterations, 40);
    out.tolerance = std::max(options.tolerance, 1.0e-8);
    out.include_phase_diagnostics = false;
    out.stability_precheck = true;
    return out;
}

void merge_stability_diagnostics(EquilibriumResultNative& result, const StabilityResultNative& stability, const EquilibriumOptionsNative& pre_opts) {
    int unstable_trial_count = 0;
    for (const auto& trial : stability.trials) {
        if (trial.unstable) {
            unstable_trial_count += 1;
        }
    }
    result.diagnostics_string["stability_analysis"] = "neutral_tpd";
    result.diagnostics_bool["stability_checked"] = true;
    result.diagnostics_bool["stability_stable"] = stability.stable;
    result.diagnostics_double["min_tpd"] = stability.min_tpd;
    result.diagnostics_string["parent_phase"] = stability.parent_phase;
    result.diagnostics_string["trial_phase"] = stability.trial_phase;
    result.diagnostics_vector["trial_composition"] = stability.trial_composition;
    result.diagnostics_int["unstable_trial_count"] = unstable_trial_count;
    result.diagnostics_int["trial_count"] = static_cast<int>(stability.trials.size());
    result.diagnostics_int["stability_max_iterations"] = pre_opts.max_iterations;
    result.diagnostics_double["stability_tolerance"] = pre_opts.tolerance;
}

void skipped_stability_diagnostics(EquilibriumResultNative& result) {
    result.diagnostics_string["stability_analysis"] = "not_run";
    result.diagnostics_bool["stability_checked"] = false;
    result.diagnostics_string["stability_message"] = "stability precheck skipped by options.stability_precheck=False";
}

bool no_split_stable_from_precheck(const EquilibriumResultNative& result) {
    auto checked = result.diagnostics_bool.find("stability_checked");
    if (checked != result.diagnostics_bool.end() && !checked->second) {
        return false;
    }
    auto stable = result.diagnostics_bool.find("stability_stable");
    return stable == result.diagnostics_bool.end() ? true : stable->second;
}

bool rachford_rice_beta(const std::vector<double>& feed, const std::vector<double>& k_values, double& beta, std::string& message) {
    auto rr = [&](double b) {
        double value = 0.0;
        for (std::size_t i = 0; i < feed.size(); ++i) {
            value += feed[i] * (k_values[i] - 1.0) / (1.0 + b * (k_values[i] - 1.0));
        }
        return value;
    };
    double f0 = rr(0.0);
    double f1 = rr(1.0);
    if (f0 <= 0.0) {
        beta = 0.0;
        message = "no two-phase Rachford-Rice bracket; liquid-like single phase";
        return false;
    }
    if (f1 >= 0.0) {
        beta = 1.0;
        message = "no two-phase Rachford-Rice bracket; vapor-like single phase";
        return false;
    }
    double lo = 0.0;
    double hi = 1.0;
    for (int i = 0; i < 100; ++i) {
        double mid = 0.5 * (lo + hi);
        if (rr(mid) > 0.0) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    beta = 0.5 * (lo + hi);
    message.clear();
    return true;
}

std::pair<std::vector<double>, std::vector<double>> phase_compositions(
    const std::vector<double>& feed,
    const std::vector<double>& k_values,
    double beta,
    double min_composition
) {
    std::vector<double> x_liq(feed.size(), 0.0);
    std::vector<double> y_vap(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        x_liq[i] = feed[i] / (1.0 + beta * (k_values[i] - 1.0));
        y_vap[i] = k_values[i] * x_liq[i];
    }
    return {clip_normalize(x_liq, min_composition), clip_normalize(y_vap, min_composition)};
}

std::vector<double> blend_lle_seed(const std::vector<double>& feed, const std::vector<double>& target, double strength, double min_composition) {
    std::vector<double> out(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        out[i] = (1.0 - strength) * feed[i] + strength * target[i];
    }
    return clip_normalize(out, min_composition);
}

std::pair<std::vector<double>, std::vector<double>> feed_perturb_guess(const std::vector<double>& feed, const EquilibriumOptionsNative& options) {
    std::vector<double> comp1 = feed;
    std::vector<double> comp2 = feed;
    double min_feed = *std::min_element(feed.begin(), feed.end());
    double delta = std::min(0.2, 0.5 * min_feed);
    if (feed.size() > 1 && delta > options.min_composition) {
        comp1[0] = std::max(options.min_composition, comp1[0] - delta);
        comp1[1] += delta;
        comp2[0] += delta;
        comp2[1] = std::max(options.min_composition, comp2[1] - delta);
        comp1 = clip_normalize(comp1, options.min_composition);
        comp2 = clip_normalize(comp2, options.min_composition);
    }
    return {comp1, comp2};
}

std::vector<LLESeedNative> default_lle_guesses(const std::vector<double>& feed, const EquilibriumOptionsNative& options) {
    std::vector<LLESeedNative> guesses;
    if (feed.size() > 1) {
        for (double strength : {0.9, 0.7}) {
            for (std::size_t component1 = 0; component1 < feed.size(); ++component1) {
                for (std::size_t component2 = 0; component2 < feed.size(); ++component2) {
                    if (component1 == component2) {
                        continue;
                    }
                    LLESeedNative seed;
                    seed.seed_name = "auto_pair_" + std::to_string(component2) + "_" + std::to_string(component1) + "_s" + std::to_string(static_cast<int>(100 * strength));
                    seed.beta = 0.5;
                    seed.comp1 = blend_lle_seed(feed, component_rich_composition(feed.size(), component2, options.min_composition), strength, options.min_composition);
                    seed.comp2 = blend_lle_seed(feed, component_rich_composition(feed.size(), component1, options.min_composition), strength, options.min_composition);
                    guesses.push_back(seed);
                }
            }
        }
    }
    auto perturb = feed_perturb_guess(feed, options);
    guesses.push_back({"auto_feed_perturb", 0.5, perturb.first, perturb.second});
    return guesses;
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

std::vector<double> pack_lle_variables(double beta, const std::vector<double>& comp1, const std::vector<double>& comp2) {
    beta = std::max(1.0e-12, std::min(1.0 - 1.0e-12, beta));
    std::vector<double> out;
    out.push_back(std::log(beta / (1.0 - beta)));
    auto logits1 = composition_to_logits(comp1);
    auto logits2 = composition_to_logits(comp2);
    out.insert(out.end(), logits1.begin(), logits1.end());
    out.insert(out.end(), logits2.begin(), logits2.end());
    return out;
}

void unpack_lle_variables(const std::vector<double>& variables, std::size_t ncomp, double& beta, std::vector<double>& comp1, std::vector<double>& comp2) {
    if (variables.size() != 1 + 2 * (ncomp - 1)) {
        throw SolutionError("Unexpected LLE variable vector size.");
    }
    beta = 1.0 / (1.0 + std::exp(-std::max(-700.0, std::min(700.0, variables[0]))));
    std::vector<double> logits1(variables.begin() + 1, variables.begin() + static_cast<std::ptrdiff_t>(ncomp));
    std::vector<double> logits2(variables.begin() + static_cast<std::ptrdiff_t>(ncomp), variables.end());
    comp1 = logits_to_composition(logits1);
    comp2 = logits_to_composition(logits2);
}

LLECandidateNative evaluate_lle_variables(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const std::vector<double>& variables,
    const EquilibriumOptionsNative& options
) {
    LLECandidateNative out;
    unpack_lle_variables(variables, feed.size(), out.beta, out.comp1, out.comp2);
    out.state1 = phase_state(mixture, t, p, out.comp1, "liq");
    out.state2 = phase_state(mixture, t, p, out.comp2, "liq");
    out.fugacity_residual.resize(feed.size(), 0.0);
    out.material_residual.resize(feed.size(), 0.0);
    out.residual.resize(feed.size() * 2, 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        out.fugacity_residual[i] = std::log(out.comp2[i]) + out.state2.ln_phi[i] - std::log(out.comp1[i]) - out.state1.ln_phi[i];
        out.material_residual[i] = (1.0 - out.beta) * out.comp1[i] + out.beta * out.comp2[i] - feed[i];
        out.residual[i] = out.fugacity_residual[i];
        out.residual[i + feed.size()] = out.material_residual[i];
    }
    out.fugacity_residual_norm = max_abs(out.fugacity_residual);
    out.material_error = max_abs(out.material_residual);
    (void)options;
    return out;
}

bool lle_degenerate(const LLECandidateNative& candidate, const EquilibriumOptionsNative& options) {
    return candidate.beta <= options.min_composition
        || candidate.beta >= 1.0 - options.min_composition
        || phase_distance(candidate.comp1, candidate.comp2) <= split_distance_tolerance(options);
}

bool lle_converged(const LLECandidateNative& candidate, const EquilibriumOptionsNative& options) {
    return candidate.fugacity_residual_norm <= options.tolerance
        && candidate.material_error <= std::max(options.tolerance, 1.0e-10)
        && !lle_degenerate(candidate, options);
}

std::vector<double> newton_step(
    const std::function<std::vector<double>(const std::vector<double>&)>& residual_fn,
    const std::vector<double>& variables,
    const std::vector<double>& residual
) {
    Eigen::MatrixXd jacobian(static_cast<Eigen::Index>(residual.size()), static_cast<Eigen::Index>(variables.size()));
    for (std::size_t column = 0; column < variables.size(); ++column) {
        double step = 1.0e-5 * std::max(1.0, std::abs(variables[column]));
        std::vector<double> forward = variables;
        std::vector<double> backward = variables;
        forward[column] += step;
        backward[column] -= step;
        std::vector<double> fwd = residual_fn(forward);
        std::vector<double> bwd = residual_fn(backward);
        for (std::size_t row = 0; row < residual.size(); ++row) {
            jacobian(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(column)) = (fwd[row] - bwd[row]) / (2.0 * step);
        }
    }
    Eigen::VectorXd rhs(static_cast<Eigen::Index>(residual.size()));
    for (std::size_t row = 0; row < residual.size(); ++row) {
        rhs(static_cast<Eigen::Index>(row)) = -residual[row];
    }
    Eigen::VectorXd delta = jacobian.colPivHouseholderQr().solve(rhs);
    std::vector<double> out(static_cast<std::size_t>(delta.size()), 0.0);
    for (Eigen::Index i = 0; i < delta.size(); ++i) {
        out[static_cast<std::size_t>(i)] = delta(i);
    }
    return out;
}

std::vector<std::pair<std::string, std::vector<double>>> deterministic_formula_multistart_variables(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const EquilibriumOptionsNative& options
) {
    std::vector<std::pair<std::string, std::vector<double>>> starts;
    for (const auto& seed : electrolyte_formula_seeds(mixture, feed, basis, options)) {
        starts.push_back({seed.first, composition_to_logits(clip_normalize(seed.second, options.min_composition))});
    }
    const std::size_t nvar = basis.formula_feed.empty() ? 0 : basis.formula_feed.size() - 1;
    starts.push_back({"logit_center", std::vector<double>(nvar, 0.0)});
    for (std::size_t i = 0; i < nvar; ++i) {
        for (double magnitude : {2.0, 4.0}) {
            std::vector<double> positive(nvar, 0.0);
            positive[i] = magnitude;
            starts.push_back({"logit_axis_" + std::to_string(i) + "_pos_" + std::to_string(static_cast<int>(magnitude)), positive});
            std::vector<double> negative(nvar, 0.0);
            negative[i] = -magnitude;
            starts.push_back({"logit_axis_" + std::to_string(i) + "_neg_" + std::to_string(static_cast<int>(magnitude)), negative});
        }
    }
    return starts;
}

std::pair<std::vector<double>, int> polish_formula_tpd_variables(
    const std::function<double(const std::vector<double>&)>& objective,
    const std::vector<double>& start,
    int max_iterations
) {
    if (start.empty() || max_iterations <= 0) {
        return {start, 0};
    }
    const std::size_t n = start.size();
    std::vector<std::vector<double>> simplex(n + 1, start);
    for (std::size_t i = 0; i < n; ++i) {
        simplex[i + 1][i] += 0.75;
    }
    std::vector<double> values(n + 1, 0.0);
    auto refresh = [&]() {
        for (std::size_t i = 0; i < simplex.size(); ++i) {
            values[i] = objective(simplex[i]);
        }
    };
    refresh();
    int iterations = 0;
    for (; iterations < max_iterations; ++iterations) {
        std::vector<std::size_t> order(simplex.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return values[a] < values[b]; });
        std::vector<std::vector<double>> sorted_simplex;
        std::vector<double> sorted_values;
        for (std::size_t idx : order) {
            sorted_simplex.push_back(simplex[idx]);
            sorted_values.push_back(values[idx]);
        }
        simplex = sorted_simplex;
        values = sorted_values;
        double spread = 0.0;
        for (double value : values) {
            spread = std::max(spread, std::abs(value - values.front()));
        }
        if (spread <= 1.0e-10) {
            break;
        }
        std::vector<double> centroid(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t row = 0; row < n; ++row) {
                centroid[i] += simplex[row][i];
            }
            centroid[i] /= static_cast<double>(n);
        }
        std::vector<double> reflected(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            reflected[i] = centroid[i] + (centroid[i] - simplex.back()[i]);
        }
        double reflected_value = objective(reflected);
        if (reflected_value < values.front()) {
            std::vector<double> expanded(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) {
                expanded[i] = centroid[i] + 2.0 * (reflected[i] - centroid[i]);
            }
            double expanded_value = objective(expanded);
            simplex.back() = expanded_value < reflected_value ? expanded : reflected;
            values.back() = std::min(expanded_value, reflected_value);
            continue;
        }
        if (reflected_value < values[n - 1]) {
            simplex.back() = reflected;
            values.back() = reflected_value;
            continue;
        }
        std::vector<double> contracted(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            contracted[i] = centroid[i] + 0.5 * (simplex.back()[i] - centroid[i]);
        }
        double contracted_value = objective(contracted);
        if (contracted_value < values.back()) {
            simplex.back() = contracted;
            values.back() = contracted_value;
            continue;
        }
        for (std::size_t row = 1; row < simplex.size(); ++row) {
            for (std::size_t i = 0; i < n; ++i) {
                simplex[row][i] = simplex.front()[i] + 0.5 * (simplex[row][i] - simplex.front()[i]);
            }
        }
        refresh();
    }
    std::size_t best = static_cast<std::size_t>(std::min_element(values.begin(), values.end()) - values.begin());
    return {simplex[best], iterations};
}

ElectrolyteTpdSearchNative electrolyte_tpd_global_search(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const PhaseStateNative& parent,
    const ElectrolyteBasisNative& basis,
    const EquilibriumOptionsNative& options,
    double threshold
) {
    ElectrolyteTpdSearchNative search;
    auto starts = deterministic_formula_multistart_variables(mixture, feed, basis, options);
    search.multistart_count = static_cast<int>(starts.size());
    int trial_index = 0;
    auto objective = [&](const std::vector<double>& variables) {
        try {
            return electrolyte_trial_from_formula(
                mixture,
                t,
                p,
                feed,
                parent,
                basis,
                logits_to_composition(variables),
                "tpd_objective",
                1,
                threshold
            ).tpd;
        } catch (const std::exception&) {
            return 1.0e6;
        }
    };
    for (const auto& start : starts) {
        trial_index += 1;
        search.trials.push_back(electrolyte_trial_from_formula(
            mixture,
            t,
            p,
            feed,
            parent,
            basis,
            logits_to_composition(start.second),
            start.first,
            trial_index,
            threshold
        ));
        auto polished = polish_formula_tpd_variables(objective, start.second, options.max_iterations);
        search.polish_iterations += polished.second;
        trial_index += std::max(1, polished.second);
        search.trials.push_back(electrolyte_trial_from_formula(
            mixture,
            t,
            p,
            feed,
            parent,
            basis,
            logits_to_composition(polished.first),
            "polished_" + start.first,
            trial_index,
            threshold
        ));
    }
    return search;
}

LLEAttemptNative solve_lle_attempt(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const LLESeedNative& seed,
    int attempt_count
) {
    LLEAttemptNative attempt;
    attempt.seed_name = seed.seed_name;
    attempt.attempt_count = attempt_count;
    if (phase_distance(seed.comp1, seed.comp2) <= split_distance_tolerance(options)) {
        attempt.status = "degenerate";
        attempt.message = "no V2 LLE split found; initial liquid phases are compositionally identical";
        attempt.candidate.beta = seed.beta;
        attempt.candidate.comp1 = seed.comp1;
        attempt.candidate.comp2 = seed.comp2;
        attempt.candidate.iteration = 0;
        attempt.candidate.objective = 0.0;
        return attempt;
    }
    std::vector<double> variables = pack_lle_variables(seed.beta, seed.comp1, seed.comp2);
    LLECandidateNative best;
    double best_objective = std::numeric_limits<double>::infinity();
    bool has_best = false;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        LLECandidateNative current = evaluate_lle_variables(mixture, t, p, feed, variables, options);
        double objective = l2_norm(current.residual);
        current.iteration = iteration;
        current.objective = objective;
        if (objective < best_objective) {
            best = current;
            best_objective = objective;
            has_best = true;
        }
        if (lle_converged(current, options)) {
            attempt.status = "converged";
            attempt.candidate = current;
            return attempt;
        }
        if (lle_degenerate(current, options)) {
            attempt.status = "degenerate";
            attempt.message = "no V2 LLE split found; phase split collapsed during iteration";
            attempt.candidate = current;
            return attempt;
        }
        auto residual_fn = [&](const std::vector<double>& candidate) {
            return evaluate_lle_variables(mixture, t, p, feed, candidate, options).residual;
        };
        std::vector<double> step = newton_step(residual_fn, variables, current.residual);
        bool accepted = false;
        for (double scale : damping_schedule(options.damping)) {
            std::vector<double> candidate_vars = variables;
            for (std::size_t i = 0; i < candidate_vars.size(); ++i) {
                candidate_vars[i] += scale * step[i];
            }
            LLECandidateNative candidate = evaluate_lle_variables(mixture, t, p, feed, candidate_vars, options);
            if (l2_norm(candidate.residual) < objective) {
                variables = candidate_vars;
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            attempt.status = "failed";
            attempt.message = "residual improvement stalled";
            attempt.candidate = has_best ? best : current;
            return attempt;
        }
    }
    attempt.candidate = best;
    if (lle_degenerate(best, options)) {
        attempt.status = "degenerate";
        attempt.message = "no V2 LLE split found; best candidate collapsed to one liquid phase";
    } else {
        attempt.status = "failed";
        attempt.message = "maximum iterations reached";
    }
    return attempt;
}

EquilibriumResultNative lle_two_phase_result(
    const LLECandidateNative& best,
    const EquilibriumOptionsNative& options,
    const std::string& seed_name,
    int attempt_count
) {
    EquilibriumResultNative result;
    result.backend = "neutral_lle";
    result.problem_kind = "lle_flash";
    result.stable = false;
    result.split_detected = true;
    result.phases.push_back(phase_from_state("liq1", best.comp1, 1.0 - best.beta, best.state1, options.include_phase_diagnostics));
    result.phases.push_back(phase_from_state("liq2", best.comp2, best.beta, best.state2, options.include_phase_diagnostics));
    result.diagnostics_int["iterations"] = best.iteration;
    result.diagnostics_double["fugacity_residual_norm"] = best.fugacity_residual_norm;
    result.diagnostics_vector["fugacity_residual"] = best.fugacity_residual;
    result.diagnostics_double["material_balance_error"] = best.material_error;
    result.diagnostics_double["liquid2_phase_fraction"] = best.beta;
    result.diagnostics_double["phase_distance"] = phase_distance(best.comp1, best.comp2);
    result.diagnostics_string["seed_name"] = seed_name;
    result.diagnostics_int["attempt_count"] = attempt_count;
    result.diagnostics_string["message"] = "converged";
    result.diagnostics_bool["point_solver_split_detected"] = true;
    result.diagnostics_string["point_solver_message"] = "converged";
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_string["nonlinear_solver"] = "native_finite_difference_newton";
    result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
    result.diagnostics_string["jacobian_backend"] = "finite_difference";
    result.diagnostics_bool["jacobian_available"] = true;
    if (options.jacobian_backend == "auto" || options.jacobian_backend == "autodiff") {
        result.diagnostics_bool["finite_difference_fallback_used"] = true;
        result.diagnostics_bool["jacobian_fallback_used"] = true;
        result.diagnostics_string["finite_difference_fallback_reason"] =
            "autodiff neutral LLE residual jacobian is not implemented for native state calls yet";
        result.diagnostics_string["jacobian_fallback_reason"] =
            "autodiff neutral LLE residual jacobian is not implemented for native state calls yet";
    } else {
        result.diagnostics_bool["finite_difference_fallback_used"] = false;
        result.diagnostics_bool["jacobian_fallback_used"] = false;
        result.diagnostics_string["jacobian_fallback_reason"] = "";
    }
    result.diagnostics_bool["hessian_available"] = false;
    result.diagnostics_string["hessian_backend"] = "not_implemented";
    result.diagnostics_bool["hessian_fallback_used"] = false;
    result.diagnostics_string["hessian_fallback_reason"] =
        "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
    return result;
}

EquilibriumResultNative lle_no_split_result(
    const std::vector<double>& feed,
    const PhaseStateNative& feed_state,
    const EquilibriumOptionsNative& options,
    const LLEAttemptNative& attempt
) {
    EquilibriumResultNative result;
    result.backend = "neutral_lle";
    result.problem_kind = "lle_flash";
    result.split_detected = false;
    result.phases.push_back(phase_from_state("liq", feed, 1.0, feed_state, options.include_phase_diagnostics));
    result.diagnostics_int["iterations"] = attempt.candidate.iteration;
    result.diagnostics_double["fugacity_residual_norm"] = attempt.candidate.fugacity_residual_norm;
    result.diagnostics_double["material_balance_error"] = attempt.candidate.material_error;
    result.diagnostics_double["liquid2_phase_fraction"] = attempt.candidate.beta;
    result.diagnostics_double["phase_distance"] = phase_distance(attempt.candidate.comp1, attempt.candidate.comp2);
    result.diagnostics_string["seed_name"] = attempt.seed_name;
    result.diagnostics_int["attempt_count"] = attempt.attempt_count;
    result.diagnostics_string["message"] = attempt.message;
    result.diagnostics_bool["point_solver_split_detected"] = false;
    result.diagnostics_string["point_solver_message"] = attempt.message;
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
    result.diagnostics_string["jacobian_backend"] = "finite_difference";
    result.diagnostics_bool["jacobian_available"] = true;
    if (options.jacobian_backend == "auto" || options.jacobian_backend == "autodiff") {
        result.diagnostics_bool["finite_difference_fallback_used"] = true;
        result.diagnostics_bool["jacobian_fallback_used"] = true;
        result.diagnostics_string["finite_difference_fallback_reason"] =
            "autodiff neutral LLE residual jacobian is not implemented for native state calls yet";
        result.diagnostics_string["jacobian_fallback_reason"] =
            "autodiff neutral LLE residual jacobian is not implemented for native state calls yet";
    } else {
        result.diagnostics_bool["finite_difference_fallback_used"] = false;
        result.diagnostics_bool["jacobian_fallback_used"] = false;
        result.diagnostics_string["jacobian_fallback_reason"] = "";
    }
    result.diagnostics_bool["hessian_available"] = false;
    result.diagnostics_string["hessian_backend"] = "not_implemented";
    result.diagnostics_bool["hessian_fallback_used"] = false;
    result.diagnostics_string["hessian_fallback_reason"] =
        "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
    return result;
}

}  // namespace

StabilityResultNative neutral_stability_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& raw_feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& parent_phases,
    const std::vector<std::string>& trial_phases
) {
    if (mixture->has_ionic()) {
        throw ValueError("Neutral stability does not support ion-containing mixtures.");
    }
    std::vector<double> feed = normalize_feed(raw_feed, mixture->ncomp(), options.min_composition, "stability");
    std::vector<StabilityTrialNative> trials;
    double threshold = std::max(options.tolerance, 1.0e-8);
    for (const std::string& parent_phase : parent_phases) {
        PhaseStateNative parent = phase_state(mixture, t, p, feed, parent_phase);
        for (const std::string& trial_phase : trial_phases) {
            for (const auto& seed : tpd_seeds(feed, options)) {
                trials.push_back(solve_tpd_trial(mixture, t, p, feed, parent.ln_phi, parent_phase, trial_phase, seed.first, seed.second, options, threshold));
            }
        }
    }
    if (trials.empty()) {
        throw ValueError("stability analysis requires at least one parent and trial phase.");
    }
    auto best_iter = std::min_element(trials.begin(), trials.end(), [](const auto& a, const auto& b) {
        return a.tpd < b.tpd;
    });
    StabilityResultNative result;
    result.backend = "neutral_tpd";
    result.problem_kind = "stability";
    result.trials = trials;
    result.stable = best_iter->tpd >= -threshold;
    result.min_tpd = best_iter->tpd;
    result.parent_phase = best_iter->parent_phase;
    result.trial_phase = best_iter->trial_phase;
    result.trial_composition = best_iter->composition;
    result.diagnostics_string["stability_analysis"] = "neutral_tpd";
    result.diagnostics_bool["stable"] = result.stable;
    result.diagnostics_double["tpd_threshold"] = threshold;
    result.diagnostics_double["min_tpd"] = result.min_tpd;
    result.diagnostics_int["trial_count"] = static_cast<int>(trials.size());
    result.diagnostics_string["min_seed_name"] = best_iter->seed_name;
    result.diagnostics_string["message"] = result.stable ? "no negative TPD trial found" : "unstable trial phase detected";
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_string["tpd_method"] = "native_tpd_fixed_point";
    return result;
}

StabilityResultNative electrolyte_stability_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& raw_feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& species
) {
    if (!mixture->has_ionic()) {
        throw ValueError("electrolyte_stability requires an ion-containing mixture.");
    }
    std::vector<double> feed = normalize_feed(raw_feed, mixture->ncomp(), options.min_composition, "electrolyte_stability");
    const std::vector<double>& charges = mixture->args().z;
    double feed_charge = composition_charge(feed, charges);
    if (std::abs(feed_charge) > 1.0e-10) {
        throw ValueError("electrolyte_stability feed must be charge neutral.");
    }
    ElectrolyteBasisNative basis = build_electrolyte_basis_native(mixture, feed, species);
    PhaseStateNative parent = phase_state(mixture, t, p, feed, "liq");
    double threshold = std::max(options.tolerance, 1.0e-8);
    ElectrolyteTpdSearchNative search = electrolyte_tpd_global_search(mixture, t, p, feed, parent, basis, options, threshold);
    if (search.trials.empty()) {
        throw SolutionError("electrolyte TPD could not generate a trial phase.");
    }
    auto best_iter = std::min_element(search.trials.begin(), search.trials.end(), [](const auto& a, const auto& b) {
        return a.tpd < b.tpd;
    });
    StabilityResultNative result;
    result.backend = "electrolyte_tpd";
    result.problem_kind = "electrolyte_stability";
    result.trials = search.trials;
    result.stable = best_iter->tpd >= -threshold;
    result.min_tpd = best_iter->tpd;
    result.parent_phase = "liq";
    result.trial_phase = best_iter->trial_phase;
    result.trial_composition = best_iter->composition;
    result.diagnostics_string["stability_analysis"] = "electrolyte_tpd";
    result.diagnostics_bool["stability_checked"] = true;
    result.diagnostics_bool["stability_stable"] = result.stable;
    result.diagnostics_string["tpd_method"] = "native_tpd_global_search";
    result.diagnostics_string["variable_model"] = basis.variable_model;
    result.diagnostics_string["seed_name"] = best_iter->seed_name;
    result.diagnostics_string["tpd_best_seed_name"] = best_iter->seed_name;
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_double["tpd_threshold"] = threshold;
    result.diagnostics_double["tpd_objective_value"] = best_iter->tpd;
    result.diagnostics_double["min_tpd"] = best_iter->tpd;
    result.diagnostics_double["phase_charge_balance_feed"] = feed_charge;
    result.diagnostics_double["phase_charge_balance_trial"] = composition_charge(best_iter->composition, charges);
    result.diagnostics_int["basis_rank"] = basis.basis_rank;
    result.diagnostics_int["trial_count"] = static_cast<int>(search.trials.size());
    result.diagnostics_int["tpd_trial_count"] = static_cast<int>(search.trials.size());
    result.diagnostics_int["tpd_multistart_count"] = search.multistart_count;
    result.diagnostics_int["tpd_polish_iterations"] = search.polish_iterations;
    result.diagnostics_int["requested_max_iterations"] = options.max_iterations;
    result.diagnostics_int["effective_max_iterations"] = options.max_iterations;
    return result;
}

EquilibriumResultNative tp_flash_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& raw_feed,
    const EquilibriumOptionsNative& options
) {
    if (mixture->has_ionic()) {
        throw ValueError("Neutral equilibrium does not support ion-containing mixtures.");
    }
    std::vector<double> feed = normalize_feed(raw_feed, mixture->ncomp(), options.min_composition, "tp_flash");
    EquilibriumResultNative result;
    result.backend = "neutral_vle";
    result.problem_kind = "tp_flash";

    if (options.stability_precheck) {
        EquilibriumOptionsNative stab_opts = precheck_options(options);
        StabilityResultNative stability = neutral_stability_native(mixture, t, p, feed, stab_opts, {"liq", "vap"}, {"liq", "vap"});
        merge_stability_diagnostics(result, stability, stab_opts);
        mixture->clear_runtime_caches();
    } else {
        skipped_stability_diagnostics(result);
    }

    PhaseStateNative liquid_seed = phase_state(mixture, t, p, feed, "liq");
    PhaseStateNative vapor_seed = phase_state(mixture, t, p, feed, "vap");
    std::vector<double> ln_k(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        ln_k[i] = liquid_seed.ln_phi[i] - vapor_seed.ln_phi[i];
    }

    double beta = 0.0;
    std::string no_split_message;
    std::vector<double> k_values(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        k_values[i] = std::exp(ln_k[i]);
    }
    if (!rachford_rice_beta(feed, k_values, beta, no_split_message)) {
        std::string phase_label = beta <= 0.0 ? "liq" : "vap";
        const PhaseStateNative& phase_state_payload = beta <= 0.0 ? liquid_seed : vapor_seed;
        result.phases.push_back(phase_from_state(phase_label, feed, 1.0, phase_state_payload, options.include_phase_diagnostics));
        result.stable = no_split_stable_from_precheck(result);
        result.split_detected = false;
        result.diagnostics_int["iterations"] = 0;
        result.diagnostics_double["fugacity_residual_norm"] = 0.0;
        result.diagnostics_double["material_balance_error"] = 0.0;
        result.diagnostics_double["vapor_fraction"] = beta;
        result.diagnostics_string["message"] = no_split_message;
        result.diagnostics_bool["point_solver_split_detected"] = false;
        result.diagnostics_string["point_solver_message"] = no_split_message;
        result.diagnostics_string["solver_language"] = "c++";
        result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
        return result;
    }

    std::vector<double> best_x;
    std::vector<double> best_y;
    PhaseStateNative best_liquid;
    PhaseStateNative best_vapor;
    std::vector<double> best_residual;
    double best_residual_norm = std::numeric_limits<double>::infinity();
    double best_material_error = std::numeric_limits<double>::infinity();
    int best_iteration = 0;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        for (std::size_t i = 0; i < feed.size(); ++i) {
            k_values[i] = std::exp(ln_k[i]);
        }
        if (!rachford_rice_beta(feed, k_values, beta, no_split_message)) {
            throw SolutionError("TP flash lost its two-phase Rachford-Rice bracket after iteration.");
        }
        auto comps = phase_compositions(feed, k_values, beta, options.min_composition);
        PhaseStateNative liquid = phase_state(mixture, t, p, comps.first, "liq");
        PhaseStateNative vapor = phase_state(mixture, t, p, comps.second, "vap");
        std::vector<double> fugacity_residual(feed.size(), 0.0);
        std::vector<double> material_residual(feed.size(), 0.0);
        for (std::size_t i = 0; i < feed.size(); ++i) {
            fugacity_residual[i] = std::log(comps.second[i]) + vapor.ln_phi[i] - std::log(comps.first[i]) - liquid.ln_phi[i];
            material_residual[i] = (1.0 - beta) * comps.first[i] + beta * comps.second[i] - feed[i];
        }
        best_x = comps.first;
        best_y = comps.second;
        best_liquid = liquid;
        best_vapor = vapor;
        best_residual = fugacity_residual;
        best_residual_norm = max_abs(fugacity_residual);
        best_material_error = max_abs(material_residual);
        best_iteration = iteration;
        if (best_residual_norm <= options.tolerance && best_material_error <= std::max(options.tolerance, 1.0e-10)) {
            result.phases.push_back(phase_from_state("liq", best_x, 1.0 - beta, best_liquid, options.include_phase_diagnostics));
            result.phases.push_back(phase_from_state("vap", best_y, beta, best_vapor, options.include_phase_diagnostics));
            result.stable = false;
            result.split_detected = true;
            result.diagnostics_int["iterations"] = best_iteration;
            result.diagnostics_double["fugacity_residual_norm"] = best_residual_norm;
            result.diagnostics_vector["fugacity_residual"] = best_residual;
            result.diagnostics_double["material_balance_error"] = best_material_error;
            result.diagnostics_double["vapor_fraction"] = beta;
            result.diagnostics_string["message"] = "converged";
            result.diagnostics_bool["point_solver_split_detected"] = true;
            result.diagnostics_string["point_solver_message"] = "converged";
            result.diagnostics_string["solver_language"] = "c++";
            result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
            result.diagnostics_string["nonlinear_solver"] = "native_successive_substitution";
            return result;
        }
        for (std::size_t i = 0; i < feed.size(); ++i) {
            double ln_k_target = liquid.ln_phi[i] - vapor.ln_phi[i];
            ln_k[i] = (1.0 - options.damping) * ln_k[i] + options.damping * ln_k_target;
        }
    }
    throw SolutionError("neutral TP flash did not converge after native iterations.");
}

EquilibriumResultNative electrolyte_bubble_pressure_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    const std::vector<double>& raw_x_liq,
    const ElectrolyteBubbleOptionsNative& options,
    const std::vector<std::string>& species,
    const std::vector<std::string>& vapor_species
) {
    if (!mixture->has_ionic()) {
        throw ValueError("electrolyte_bubble_pressure requires an ion-containing mixture.");
    }
    if (vapor_species.empty()) {
        throw ValueError("electrolyte_bubble_pressure requires at least one neutral vapor species.");
    }
    if (options.min_pressure <= 0.0 || options.max_pressure <= options.min_pressure || options.initial_pressure <= 0.0) {
        throw ValueError("electrolyte_bubble_pressure pressure bounds and initial pressure must be positive and ordered.");
    }
    if (options.pressure_factor <= 1.0) {
        throw ValueError("electrolyte_bubble_pressure pressure_factor must be greater than 1.");
    }
    const std::vector<double>& charges = mixture->args().z;
    if (charges.size() != mixture->ncomp()) {
        throw ValueError("mixture parameters must include one charge value per species in params['z'].");
    }
    std::vector<double> x_liq = normalize_nonnegative_composition(
        raw_x_liq,
        mixture->ncomp(),
        options.min_composition,
        "x_liq"
    );
    double charge_residual = composition_charge(x_liq, charges);
    if (std::abs(charge_residual) > options.charge_tolerance) {
        throw ValueError("electrolyte_bubble_pressure liquid composition must be charge neutral.");
    }
    std::vector<int> vapor_indices = species_indices_from_labels(species, vapor_species, "vapor_species");
    for (int index : vapor_indices) {
        if (std::abs(charges[static_cast<std::size_t>(index)]) > 1.0e-12) {
            throw ValueError("electrolyte_bubble_pressure vapor_species must be neutral; ions are liquid-only.");
        }
    }
    std::shared_ptr<ePCSAFTMixtureNative> vapor_mix = vapor_submixture(mixture, vapor_indices);
    std::vector<double> y_seed(vapor_indices.size(), 1.0 / static_cast<double>(vapor_indices.size()));
    if (!options.initial_y_vap.empty()) {
        if (options.initial_y_vap.size() != vapor_indices.size()) {
            throw ValueError("initial_y_vap length must match vapor_species length.");
        }
        y_seed = clip_normalize(options.initial_y_vap, options.min_composition);
    }

    auto evaluate = [&](double pressure, const std::vector<double>& seed) {
        return evaluate_electrolyte_bubble_pressure(mixture, vapor_mix, t, pressure, x_liq, vapor_indices, options, seed);
    };

    int state_failure_count = 0;
    std::vector<double> history_p;
    std::vector<double> history_objective;
    auto append_history = [&](const ElectrolyteBubbleEvaluationNative& eval) {
        if (history_p.size() < 40) {
            history_p.push_back(eval.p);
            history_objective.push_back(eval.finite ? eval.objective : 1.0e300);
        }
    };
    auto better = [](const ElectrolyteBubbleEvaluationNative& a, const ElectrolyteBubbleEvaluationNative& b) {
        if (!a.finite) {
            return false;
        }
        if (!b.finite) {
            return true;
        }
        return std::abs(a.objective) + 1.0e-5 < std::abs(b.objective);
    };

    ElectrolyteBubbleEvaluationNative best;
    double initial_p = std::clamp(options.initial_pressure, options.min_pressure, options.max_pressure);
    ElectrolyteBubbleEvaluationNative initial = evaluate(initial_p, y_seed);
    if (!initial.finite) {
        ++state_failure_count;
    }
    append_history(initial);
    if (better(initial, best)) {
        best = initial;
        y_seed = initial.y_vap;
    }

    ElectrolyteBubbleEvaluationNative low = initial;
    ElectrolyteBubbleEvaluationNative high = initial;
    bool bracketed = false;
    double log_step = std::log(options.pressure_factor);
    for (int expansion = 1; expansion <= options.max_bracket_expansions; ++expansion) {
        double low_p = std::max(options.min_pressure, initial_p * std::exp(-log_step * expansion));
        double high_p = std::min(options.max_pressure, initial_p * std::exp(log_step * expansion));
        ElectrolyteBubbleEvaluationNative low_eval = evaluate(low_p, y_seed);
        ElectrolyteBubbleEvaluationNative high_eval = evaluate(high_p, y_seed);
        if (!low_eval.finite) {
            ++state_failure_count;
        }
        if (!high_eval.finite) {
            ++state_failure_count;
        }
        append_history(low_eval);
        append_history(high_eval);
        if (better(low_eval, best)) {
            best = low_eval;
        }
        if (better(high_eval, best)) {
            best = high_eval;
        }
        if (low_eval.finite && high_eval.finite && low_eval.objective * high_eval.objective <= 0.0) {
            low = low_eval;
            high = high_eval;
            bracketed = true;
            break;
        }
        if (low_p <= options.min_pressure && high_p >= options.max_pressure) {
            break;
        }
    }

    bool converged = false;
    std::string message = "electrolyte bubble pressure did not bracket a pressure root";
    int iterations = 0;
    if (bracketed) {
        message = "electrolyte bubble pressure did not converge";
        double log_low = std::log(low.p);
        double log_high = std::log(high.p);
        for (iterations = 1; iterations <= options.max_iterations; ++iterations) {
            double log_mid = 0.5 * (log_low + log_high);
            ElectrolyteBubbleEvaluationNative mid = evaluate(std::exp(log_mid), best.finite ? best.y_vap : y_seed);
            if (!mid.finite) {
                ++state_failure_count;
                append_history(mid);
                log_high = log_mid;
                continue;
            }
            append_history(mid);
            if (better(mid, best)) {
                best = mid;
            }
            if (std::abs(mid.objective) <= options.tolerance || mid.residual_norm <= options.tolerance) {
                best = mid;
                converged = true;
                message = "converged";
                break;
            }
            if (low.objective * mid.objective <= 0.0) {
                high = mid;
                log_high = log_mid;
            } else {
                low = mid;
                log_low = log_mid;
            }
        }
    }
    double acceptance_tolerance = std::max(options.tolerance, 5.0e-2);
    bool accepted_by_diagnostic_envelope = false;
    if (!converged && best.finite && best.residual_norm <= acceptance_tolerance) {
        converged = true;
        accepted_by_diagnostic_envelope = true;
        message = "accepted by electrolyte bubble diagnostic envelope";
    }

    EquilibriumResultNative result;
    result.backend = "electrolyte_vle";
    result.problem_kind = "electrolyte_bubble_pressure";
    result.stable = false;
    result.split_detected = converged;
    result.diagnostics_bool["success"] = converged;
    result.diagnostics_bool["bracketed"] = bracketed;
    result.diagnostics_bool["accepted_by_diagnostic_envelope"] = accepted_by_diagnostic_envelope;
    result.diagnostics_int["iterations"] = iterations;
    result.diagnostics_int["state_failure_count"] = state_failure_count;
    result.diagnostics_int["vapor_species_count"] = static_cast<int>(vapor_indices.size());
    result.diagnostics_double["charge_residual"] = charge_residual;
    result.diagnostics_double["best_P"] = best.finite ? best.p : 1.0e300;
    result.diagnostics_double["best_objective"] = best.finite ? best.objective : 1.0e300;
    result.diagnostics_double["best_fugacity_residual_norm"] = best.finite ? best.residual_norm : 1.0e300;
    result.diagnostics_double["fugacity_residual_norm"] = best.finite ? best.residual_norm : 1.0e300;
    result.diagnostics_double["acceptance_tolerance"] = acceptance_tolerance;
    result.diagnostics_string["message"] = message;
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_electrolyte_bubble_native";
    result.diagnostics_string["nonlinear_solver"] = "native_log_pressure_bisection";
    result.diagnostics_vector["vapor_species_indices"] = std::vector<double>(vapor_indices.begin(), vapor_indices.end());
    result.diagnostics_vector["vapor_history_pressure"] = history_p;
    result.diagnostics_vector["vapor_history_objective"] = history_objective;
    if (best.finite) {
        result.phases.push_back(phase_from_state("liq", x_liq, 1.0, best.liquid, false));
        result.phases.push_back(phase_from_state("vap", best.y_full, 0.0, best.vapor, false));
        result.diagnostics_vector["best_y_vap"] = best.y_vap;
        result.diagnostics_vector["best_partial_pressures"] = best.partial_pressures;
        result.diagnostics_vector["fugacity_residual"] = best.fugacity_residual;
    }
    return result;
}

EquilibriumResultNative lle_flash_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& raw_feed,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& initial_liq1,
    const std::vector<double>& initial_liq2,
    double initial_beta,
    bool has_initial_phases
) {
    if (mixture->has_ionic()) {
        throw ValueError("Neutral equilibrium does not support ion-containing mixtures.");
    }
    std::vector<double> feed = normalize_feed(raw_feed, mixture->ncomp(), options.min_composition, "lle_flash");
    EquilibriumResultNative result_prefix;
    StabilityResultNative stability;
    EquilibriumOptionsNative stab_opts;
    std::vector<double> unstable_trial_composition;
    if (options.stability_precheck) {
        stab_opts = precheck_options(options);
        stability = neutral_stability_native(mixture, t, p, feed, stab_opts, {"liq"}, {"liq"});
        merge_stability_diagnostics(result_prefix, stability, stab_opts);
        mixture->clear_runtime_caches();
        for (const auto& trial : stability.trials) {
            if (trial.unstable && (unstable_trial_composition.empty() || trial.tpd < stability.min_tpd + 1.0e-16)) {
                unstable_trial_composition = trial.composition;
            }
        }
    } else {
        skipped_stability_diagnostics(result_prefix);
    }

    std::vector<LLESeedNative> seeds;
    if (has_initial_phases) {
        if (initial_liq1.size() != feed.size() || initial_liq2.size() != feed.size()) {
            throw ValueError("initial_phases liq1/liq2 length must match mixture component count.");
        }
        if (!(initial_beta > 0.0 && initial_beta < 1.0) || !std::isfinite(initial_beta)) {
            throw ValueError("initial_phases phase_fraction must be > 0 and < 1.");
        }
        seeds.push_back({"user", initial_beta, clip_normalize(initial_liq1, options.min_composition), clip_normalize(initial_liq2, options.min_composition)});
    } else {
        if (!unstable_trial_composition.empty() && phase_distance(feed, unstable_trial_composition) > split_distance_tolerance(options)) {
            for (double beta : {0.5, 0.25, 0.75, 0.1, 0.9}) {
                std::vector<double> comp1(feed.size(), 0.0);
                bool feasible = true;
                for (std::size_t i = 0; i < feed.size(); ++i) {
                    comp1[i] = (feed[i] - beta * unstable_trial_composition[i]) / (1.0 - beta);
                    feasible = feasible && std::isfinite(comp1[i]) && comp1[i] >= options.min_composition;
                }
                if (feasible) {
                    seeds.push_back({"tpd_liq_trial", beta, clip_normalize(comp1, options.min_composition), clip_normalize(unstable_trial_composition, options.min_composition)});
                    break;
                }
            }
        }
        auto defaults = default_lle_guesses(feed, options);
        seeds.insert(seeds.end(), defaults.begin(), defaults.end());
    }
    PhaseStateNative feed_state = phase_state(mixture, t, p, feed, "liq");
    std::vector<LLEAttemptNative> degenerate_attempts;
    std::vector<LLEAttemptNative> failed_attempts;
    int attempt_count = 0;
    for (const auto& seed : seeds) {
        attempt_count += 1;
        LLEAttemptNative attempt = solve_lle_attempt(mixture, t, p, feed, options, seed, attempt_count);
        if (attempt.status == "converged") {
            EquilibriumResultNative result = lle_two_phase_result(attempt.candidate, options, attempt.seed_name, attempt.attempt_count);
            result.diagnostics_double.insert(result_prefix.diagnostics_double.begin(), result_prefix.diagnostics_double.end());
            result.diagnostics_int.insert(result_prefix.diagnostics_int.begin(), result_prefix.diagnostics_int.end());
            result.diagnostics_bool.insert(result_prefix.diagnostics_bool.begin(), result_prefix.diagnostics_bool.end());
            result.diagnostics_string.insert(result_prefix.diagnostics_string.begin(), result_prefix.diagnostics_string.end());
            result.diagnostics_vector.insert(result_prefix.diagnostics_vector.begin(), result_prefix.diagnostics_vector.end());
            return result;
        }
        if (attempt.status == "degenerate") {
            degenerate_attempts.push_back(attempt);
        } else {
            failed_attempts.push_back(attempt);
        }
    }
    if (!failed_attempts.empty()) {
        auto best_failed = std::min_element(failed_attempts.begin(), failed_attempts.end(), [](const auto& a, const auto& b) {
            return a.candidate.objective < b.candidate.objective;
        });
        throw SolutionError("neutral LLE flash did not converge after native attempt(s); best_seed=" + best_failed->seed_name + ", reason=" + best_failed->message);
    }
    auto best_degenerate = std::min_element(degenerate_attempts.begin(), degenerate_attempts.end(), [](const auto& a, const auto& b) {
        return a.candidate.objective < b.candidate.objective;
    });
    EquilibriumResultNative result = lle_no_split_result(feed, feed_state, options, *best_degenerate);
    result.stable = no_split_stable_from_precheck(result_prefix);
    result.diagnostics_double.insert(result_prefix.diagnostics_double.begin(), result_prefix.diagnostics_double.end());
    result.diagnostics_int.insert(result_prefix.diagnostics_int.begin(), result_prefix.diagnostics_int.end());
    result.diagnostics_bool.insert(result_prefix.diagnostics_bool.begin(), result_prefix.diagnostics_bool.end());
    result.diagnostics_string.insert(result_prefix.diagnostics_string.begin(), result_prefix.diagnostics_string.end());
    result.diagnostics_vector.insert(result_prefix.diagnostics_vector.begin(), result_prefix.diagnostics_vector.end());
    return result;
}

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

std::vector<double> predictive_infeasible_residual(const std::vector<double>& variables, const ElectrolyteBasisNative& basis, const EquilibriumOptionsNative& options) {
    std::vector<double> residual(basis.formula_feed.size(), 1.0e3);
    try {
        double beta_formula = 0.5;
        std::vector<double> org_formula;
        unpack_predictive_electrolyte_variables(variables, basis.formula_feed.size(), beta_formula, org_formula);
        for (std::size_t i = 0; i < basis.formula_feed.size(); ++i) {
            double aq = (basis.formula_feed[i] - beta_formula * org_formula[i]) / (1.0 - beta_formula);
            residual[i] = 10.0 + 1.0e4 * std::max(options.min_composition - aq, 0.0);
        }
    } catch (const std::exception&) {
    }
    return residual;
}

bool predictive_electrolyte_accepted(const ElectrolyteCandidateNative& candidate, const EquilibriumOptionsNative& options) {
    return candidate.solver_residual_norm <= options.tolerance
        && candidate.material_error <= std::max(options.tolerance, 1.0e-10)
        && candidate.charge_balance_error <= 1.0e-8
        && candidate.gibbs_delta < 0.0
        && candidate.phase_distance_value > std::max(1.0e-4, split_distance_tolerance(options));
}

std::string electrolyte_failure_reason(const ElectrolyteCandidateNative& candidate, const EquilibriumOptionsNative& options, bool has_best) {
    if (!has_best) {
        return "no electrolyte LLE candidate was evaluated";
    }
    if (candidate.solver_residual_norm > options.tolerance) {
        return "nonlinear residual did not converge";
    }
    if (candidate.material_error > std::max(options.tolerance, 1.0e-10)) {
        return "material balance residual did not converge";
    }
    if (candidate.charge_balance_error > 1.0e-8) {
        return "charge balance residual did not converge";
    }
    if (candidate.gibbs_delta >= 0.0) {
        return "Gibbs split was not favored by current thermodynamic surface";
    }
    return "candidate collapsed to one phase";
}

EquilibriumAttemptDiagnosticsNative electrolyte_attempt_diagnostics(
    const std::string& seed_name,
    const ElectrolyteCandidateNative& candidate,
    const EquilibriumOptionsNative& options,
    bool accepted
) {
    EquilibriumAttemptDiagnosticsNative out;
    out.seed_name = seed_name;
    out.rejection_reason = accepted ? "accepted" : electrolyte_failure_reason(candidate, options, true);
    out.beta_org = candidate.beta_org;
    out.phase_distance = candidate.phase_distance_value;
    out.solver_residual_norm = candidate.solver_residual_norm;
    out.material_balance_error = candidate.material_error;
    out.charge_balance_error = candidate.charge_balance_error;
    out.gibbs_delta = candidate.gibbs_delta;
    out.iterations = candidate.iteration;
    return out;
}

EquilibriumResultNative electrolyte_lle_failure_result(
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const StabilityResultNative& stability,
    const ElectrolyteCandidateNative& best,
    const std::string& best_seed,
    bool has_best,
    const EquilibriumOptionsNative& options,
    const std::vector<EquilibriumAttemptDiagnosticsNative>& attempts
) {
    EquilibriumResultNative result;
    result.backend = "electrolyte_lle";
    result.problem_kind = "electrolyte_lle_flash";
    result.stable = true;
    result.split_detected = false;
    std::string reason = electrolyte_failure_reason(best, options, has_best);
    result.diagnostics_string["phase_equilibrium_model"] = "electrolyte_lle_v5_native_charge_constrained_solve";
    result.diagnostics_string["equilibrium_route"] = "electrolyte_lle";
    result.diagnostics_string["route_reason"] = "ion-containing mixture";
    result.diagnostics_string["stability_analysis"] = "electrolyte_tpd";
    result.diagnostics_string["variable_model"] = basis.variable_model;
    result.diagnostics_string["acceptance_gate"] = "predictive_solve_failed";
    result.diagnostics_string["solver_seed_name"] = best_seed;
    result.diagnostics_string["solver_method"] = "native_transformed_newton";
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_string["tpd_method"] = "native_tpd_global_search";
    result.diagnostics_string["gibbs_seed_method"] = "native_golden_section";
    result.diagnostics_string["best_failure_reason"] = reason;
    result.diagnostics_string["message"] = "electrolyte LLE flash did not converge";
    result.diagnostics_string["density_diagnostics_mode"] = options.density_diagnostics;
    result.diagnostics_string["density_validity_gate"] = "not_evaluated";
    result.diagnostics_string["density_warm_start_source"] = "";
    result.diagnostics_bool["experimental_coupled_density_lle"] = options.experimental_coupled_density_lle;
    result.diagnostics_bool["coupled_density_lle_attempted"] = options.experimental_coupled_density_lle;
    result.diagnostics_bool["density_fallback_used"] = false;
    result.diagnostics_int["density_failure_count"] = 0;
    result.diagnostics_bool["stability_checked"] = true;
    result.diagnostics_bool["stability_stable"] = stability.stable;
    result.diagnostics_bool["unstable_feed_collapsed_all_candidates"] = stability.min_tpd < -std::max(options.tolerance, 1.0e-8)
        && reason == "candidate collapsed to one phase";
    result.diagnostics_int["basis_rank"] = basis.basis_rank;
    result.diagnostics_int["seed_attempt_count"] = static_cast<int>(attempts.size());
    result.diagnostics_int["iterations"] = has_best ? best.iteration : 0;
    result.diagnostics_int["repeated_stability_iterations"] = 1;
    result.diagnostics_int["requested_max_iterations"] = options.max_iterations;
    result.diagnostics_int["effective_max_iterations"] = options.max_iterations;
    result.diagnostics_double["solver_residual_norm"] = has_best ? best.solver_residual_norm : 1.0e300;
    result.diagnostics_double["fugacity_residual_norm"] = has_best ? best.solver_residual_norm : 1.0e300;
    result.diagnostics_double["material_balance_error"] = has_best ? best.material_error : 1.0e300;
    result.diagnostics_double["charge_balance_error"] = has_best ? best.charge_balance_error : 1.0e300;
    result.diagnostics_double["gibbs_feed"] = has_best ? best.gibbs_feed : 0.0;
    result.diagnostics_double["gibbs_split"] = has_best ? best.gibbs_split : 0.0;
    result.diagnostics_double["gibbs_delta"] = has_best ? best.gibbs_delta : 0.0;
    result.diagnostics_double["phase_distance"] = has_best ? best.phase_distance_value : 0.0;
    result.diagnostics_double["stability_min_tpd"] = stability.min_tpd;
    result.diagnostics_vector["feed_composition"] = feed;
    result.diagnostics_vector["fugacity_residual"] = best.residual;
    if (has_best && best.phase_distance_value > std::max(1.0e-4, split_distance_tolerance(options))) {
        result.diagnostics_vector["best_noncollapsed_candidate_aq"] = best.aq_comp;
        result.diagnostics_vector["best_noncollapsed_candidate_org"] = best.org_comp;
        result.diagnostics_string["best_noncollapsed_candidate"] = "available";
    } else {
        result.diagnostics_string["best_noncollapsed_candidate"] = "none";
    }
    result.attempt_diagnostics = attempts;
    return result;
}

std::vector<double> nelder_mead_variables(
    const std::function<double(const std::vector<double>&)>& objective,
    const std::vector<double>& start,
    int max_iterations
) {
    const std::size_t n = start.size();
    std::vector<std::vector<double>> simplex(n + 1, start);
    for (std::size_t i = 0; i < n; ++i) {
        simplex[i + 1][i] += 0.25;
    }
    std::vector<double> values(n + 1, 0.0);
    auto refresh = [&]() {
        for (std::size_t i = 0; i < simplex.size(); ++i) {
            values[i] = objective(simplex[i]);
        }
    };
    refresh();
    for (int iter = 0; iter < max_iterations; ++iter) {
        std::vector<std::size_t> order(simplex.size());
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return values[a] < values[b]; });
        std::vector<std::vector<double>> sorted_simplex;
        std::vector<double> sorted_values;
        for (std::size_t idx : order) {
            sorted_simplex.push_back(simplex[idx]);
            sorted_values.push_back(values[idx]);
        }
        simplex = sorted_simplex;
        values = sorted_values;
        std::vector<double> centroid(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                centroid[j] += simplex[i][j];
            }
        }
        for (double& value : centroid) {
            value /= static_cast<double>(n);
        }
        auto transform = [&](double scale) {
            std::vector<double> out(n, 0.0);
            for (std::size_t j = 0; j < n; ++j) {
                out[j] = centroid[j] + scale * (centroid[j] - simplex[n][j]);
            }
            return out;
        };
        std::vector<double> reflected = transform(1.0);
        double reflected_value = objective(reflected);
        if (reflected_value < values[0]) {
            std::vector<double> expanded = transform(2.0);
            double expanded_value = objective(expanded);
            if (expanded_value < reflected_value) {
                simplex[n] = expanded;
                values[n] = expanded_value;
            } else {
                simplex[n] = reflected;
                values[n] = reflected_value;
            }
        } else if (reflected_value < values[n - 1]) {
            simplex[n] = reflected;
            values[n] = reflected_value;
        } else {
            std::vector<double> contracted(n, 0.0);
            for (std::size_t j = 0; j < n; ++j) {
                contracted[j] = centroid[j] + 0.5 * (simplex[n][j] - centroid[j]);
            }
            double contracted_value = objective(contracted);
            if (contracted_value < values[n]) {
                simplex[n] = contracted;
                values[n] = contracted_value;
            } else {
                for (std::size_t i = 1; i < simplex.size(); ++i) {
                    for (std::size_t j = 0; j < n; ++j) {
                        simplex[i][j] = simplex[0][j] + 0.5 * (simplex[i][j] - simplex[0][j]);
                    }
                }
                refresh();
            }
        }
    }
    std::size_t best = static_cast<std::size_t>(std::min_element(values.begin(), values.end()) - values.begin());
    return simplex[best];
}

std::vector<double> gibbs_seed_variables_from_trial(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& trial_composition
) {
    std::vector<double> org_formula = explicit_to_formula(trial_composition, basis);
    double upper = 1.0 - 1.0e-8;
    for (std::size_t i = 0; i < basis.formula_feed.size(); ++i) {
        if (org_formula[i] > basis.formula_feed[i]) {
            upper = std::min(upper, 0.999 * basis.formula_feed[i] / org_formula[i]);
        }
    }
    if (upper <= 1.0e-8) {
        throw SolutionError("No feasible phase-fraction interval for Gibbs seed.");
    }
    auto objective = [&](double beta_formula) {
        std::vector<double> aq_formula_raw(basis.formula_feed.size(), 0.0);
        for (std::size_t i = 0; i < basis.formula_feed.size(); ++i) {
            aq_formula_raw[i] = (basis.formula_feed[i] - beta_formula * org_formula[i]) / (1.0 - beta_formula);
            if (!std::isfinite(aq_formula_raw[i]) || aq_formula_raw[i] <= options.min_composition) {
                return 1.0e6;
            }
        }
        std::vector<double> aq_formula = clip_normalize(aq_formula_raw, options.min_composition);
        auto aq_expanded = formula_to_explicit(aq_formula, basis, feed.size());
        auto org_expanded = formula_to_explicit(org_formula, basis, feed.size());
        double beta_org = beta_formula * org_expanded.second / ((1.0 - beta_formula) * aq_expanded.second + beta_formula * org_expanded.second);
        PhaseStateNative aq_state = phase_state(mixture, t, p, aq_expanded.first, "liq", "aq");
        PhaseStateNative org_state = phase_state(mixture, t, p, org_expanded.first, "liq", "org");
        return (1.0 - beta_org) * electrolyte_gibbs_proxy(aq_expanded.first, aq_state)
            + beta_org * electrolyte_gibbs_proxy(org_expanded.first, org_state);
    };
    double lo = 1.0e-8;
    double hi = upper;
    const double phi = 0.5 * (3.0 - std::sqrt(5.0));
    double x1 = lo + phi * (hi - lo);
    double x2 = hi - phi * (hi - lo);
    double f1 = objective(x1);
    double f2 = objective(x2);
    for (int iter = 0; iter < 80; ++iter) {
        if (f1 < f2) {
            hi = x2;
            x2 = x1;
            f2 = f1;
            x1 = lo + phi * (hi - lo);
            f1 = objective(x1);
        } else {
            lo = x1;
            x1 = x2;
            f1 = f2;
            x2 = hi - phi * (hi - lo);
            f2 = objective(x2);
        }
    }
    double beta_formula = 0.5 * (lo + hi);
    return pack_predictive_electrolyte_variables(beta_formula, org_formula);
}

ElectrolyteCandidateNative solve_predictive_electrolyte_attempt(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& seed_variables,
    double gibbs_feed,
    std::vector<DensitySolveDiagnostics>* density_failures
) {
    std::vector<double> variables = seed_variables;
    auto objective_fn = [&](const std::vector<double>& candidate) {
        try {
            return l2_norm(evaluate_predictive_electrolyte_variables(mixture, t, p, feed, basis, candidate, options, gibbs_feed).residual);
        } catch (const SolutionError&) {
            append_last_density_failure(mixture, density_failures);
            return l2_norm(predictive_infeasible_residual(candidate, basis, options));
        }
    };
    if (basis.formula_feed.size() > 3) {
        variables = nelder_mead_variables(objective_fn, variables, std::max(40, options.max_iterations * 4));
    }
    ElectrolyteCandidateNative best;
    bool has_best = false;
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        ElectrolyteCandidateNative current;
        try {
            current = evaluate_predictive_electrolyte_variables(mixture, t, p, feed, basis, variables, options, gibbs_feed);
        } catch (const SolutionError&) {
            append_last_density_failure(mixture, density_failures);
            std::vector<double> residual = predictive_infeasible_residual(variables, basis, options);
            current.objective = l2_norm(residual);
            current.residual = residual;
        }
        current.iteration = iteration;
        if (!has_best || current.objective < best.objective) {
            best = current;
            has_best = true;
        }
        if (has_best && predictive_electrolyte_accepted(best, options)) {
            return best;
        }
        auto residual_fn = [&](const std::vector<double>& candidate) {
            try {
                return evaluate_predictive_electrolyte_variables(mixture, t, p, feed, basis, candidate, options, gibbs_feed).residual;
            } catch (const SolutionError&) {
                append_last_density_failure(mixture, density_failures);
                return predictive_infeasible_residual(candidate, basis, options);
            }
        };
        std::vector<double> base_residual = residual_fn(variables);
        std::vector<double> step = newton_step(residual_fn, variables, base_residual);
        bool accepted = false;
        double current_objective = l2_norm(base_residual);
        for (double scale : damping_schedule(options.damping)) {
            std::vector<double> candidate_vars = variables;
            for (std::size_t i = 0; i < candidate_vars.size(); ++i) {
                candidate_vars[i] += scale * step[i];
            }
            std::vector<double> candidate_residual = residual_fn(candidate_vars);
            if (l2_norm(candidate_residual) < current_objective) {
                variables = candidate_vars;
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            return best;
        }
    }
    return best;
}

EquilibriumResultNative electrolyte_lle_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& raw_feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& species,
    const std::vector<double>& initial_aq,
    const std::vector<double>& initial_org,
    double initial_beta_org,
    bool has_initial_phases
) {
    if (!mixture->has_ionic()) {
        throw ValueError("electrolyte_lle requires an ion-containing mixture.");
    }
    std::vector<double> feed = normalize_feed(raw_feed, mixture->ncomp(), options.min_composition, "electrolyte_lle");
    const std::vector<double>& charges = mixture->args().z;
    double feed_charge = composition_charge(feed, charges);
    if (std::abs(feed_charge) > 1.0e-10) {
        throw ValueError("electrolyte_lle feed must be charge neutral.");
    }
    ElectrolyteBasisNative basis = build_electrolyte_basis_native(mixture, feed, species);
    EquilibriumOptionsNative solver_options = options;
    PhaseStateNative feed_state = phase_state(mixture, t, p, feed, "liq", "feed");
    double gibbs_feed = electrolyte_gibbs_proxy(feed, feed_state);
    std::vector<std::pair<std::string, std::vector<double>>> seed_variables;
    if (has_initial_phases) {
        if (initial_aq.size() != feed.size() || initial_org.size() != feed.size()) {
            throw ValueError("initial_phases aq/org length must match mixture component count.");
        }
        if (!(initial_beta_org > 0.0 && initial_beta_org < 1.0) || !std::isfinite(initial_beta_org)) {
            throw ValueError("initial_phases phase_fraction must be > 0 and < 1.");
        }
        std::vector<double> aq_comp_initial = clip_normalize(initial_aq, solver_options.min_composition);
        std::vector<double> org_comp_initial = clip_normalize(initial_org, solver_options.min_composition);
        double aq_charge_initial = composition_charge(aq_comp_initial, charges);
        double org_charge_initial = composition_charge(org_comp_initial, charges);
        if (std::max(std::abs(aq_charge_initial), std::abs(org_charge_initial)) > 1.0e-8) {
            throw ValueError("initial_phases aq and org must be charge neutral for electrolyte_lle.");
        }
        std::vector<double> material_residual_initial(feed.size(), 0.0);
        for (std::size_t i = 0; i < feed.size(); ++i) {
            material_residual_initial[i] = (1.0 - initial_beta_org) * aq_comp_initial[i] + initial_beta_org * org_comp_initial[i] - feed[i];
        }
        if (max_abs(material_residual_initial) > 1.0e-7) {
            throw ValueError("initial_phases aq/org/phase_fraction must reconstruct the electrolyte_lle feed.");
        }
        std::vector<double> aq_formula = explicit_to_formula(aq_comp_initial, basis);
        std::vector<double> org_formula = explicit_to_formula(org_comp_initial, basis);
        double beta_formula = initial_beta_org;
        auto aq_expanded = formula_to_explicit(aq_formula, basis, feed.size());
        auto org_expanded = formula_to_explicit(org_formula, basis, feed.size());
        double numerator = initial_beta_org / org_expanded.second;
        double denominator = numerator + (1.0 - initial_beta_org) / aq_expanded.second;
        if (denominator > 0.0) {
            beta_formula = std::max(1.0e-12, std::min(1.0 - 1.0e-12, numerator / denominator));
        }
        seed_variables.emplace_back("initial_phases", pack_predictive_electrolyte_variables(beta_formula, org_formula));
    }
    StabilityResultNative stability = electrolyte_stability_native(mixture, t, p, feed, precheck_options(solver_options), species);
    if (!stability.trial_composition.empty() && phase_distance(stability.trial_composition, feed) > split_distance_tolerance(solver_options)) {
        seed_variables.emplace_back("native_tpd_trial", pack_predictive_electrolyte_variables(0.5, explicit_to_formula(stability.trial_composition, basis)));
        try {
            seed_variables.emplace_back("native_gibbs_tpd_trial", gibbs_seed_variables_from_trial(mixture, t, p, feed, basis, solver_options, stability.trial_composition));
        } catch (const SolutionError&) {
        }
    }
    for (const auto& seed : electrolyte_formula_seeds(mixture, feed, basis, solver_options)) {
        for (double beta : {0.02, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90}) {
            seed_variables.emplace_back(seed.first + "_b" + std::to_string(static_cast<int>(100 * beta)), pack_predictive_electrolyte_variables(beta, seed.second));
        }
    }
    for (const auto& seed : electrolyte_endpoint_phase_seeds(basis, solver_options)) {
        for (double beta : {0.02, 0.05, 0.10, 0.25, 0.50}) {
            seed_variables.emplace_back(seed.first + "_b" + std::to_string(static_cast<int>(100 * beta)), pack_predictive_electrolyte_variables(beta, seed.second));
        }
    }
    ElectrolyteCandidateNative best;
    std::string best_seed;
    bool has_best = false;
    bool accepted = false;
    std::vector<EquilibriumAttemptDiagnosticsNative> attempts;
    std::vector<DensitySolveDiagnostics> density_failures;
    for (const auto& seed : seed_variables) {
        ElectrolyteCandidateNative candidate = solve_predictive_electrolyte_attempt(
            mixture,
            t,
            p,
            feed,
            basis,
            solver_options,
            seed.second,
            gibbs_feed,
            &density_failures
        );
        bool candidate_accepted = predictive_electrolyte_accepted(candidate, solver_options);
        attempts.push_back(electrolyte_attempt_diagnostics(seed.first, candidate, solver_options, candidate_accepted));
        if (!has_best || candidate.objective < best.objective) {
            best = candidate;
            best_seed = seed.first;
            has_best = true;
        }
        if (candidate_accepted) {
            best = candidate;
            best_seed = seed.first;
            accepted = true;
            break;
        }
    }
    if (!has_best || !accepted) {
        EquilibriumResultNative result = electrolyte_lle_failure_result(feed, basis, stability, best, best_seed, has_best, solver_options, attempts);
        result.density_diagnostics = density_failures;
        return result;
    }
    bool labels_swapped = false;
    std::string phase_label_basis = "composition_order";
    if (basis.neutral_indices.size() >= 2) {
        std::size_t aq_index = static_cast<std::size_t>(basis.neutral_indices[0]);
        std::string aq_name = species.size() > aq_index ? species[aq_index] : std::string("component_") + std::to_string(aq_index);
        std::string aq_name_lower = aq_name;
        std::transform(aq_name_lower.begin(), aq_name_lower.end(), aq_name_lower.begin(), [](unsigned char ch) {
            return static_cast<char>(std::tolower(ch));
        });
        if (best_seed == "initial_phases") {
            phase_label_basis = "initial_aq_org";
        } else if (aq_name_lower == "h2o" || aq_name_lower == "water") {
            phase_label_basis = "water_rich:" + aq_name;
        } else {
            phase_label_basis = "high_dielectric_rich:" + aq_name;
        }
        if (best.aq_comp[aq_index] < best.org_comp[aq_index]) {
            std::swap(best.aq_comp, best.org_comp);
            std::swap(best.aq_formula, best.org_formula);
            std::swap(best.aq_state, best.org_state);
            best.beta_org = 1.0 - best.beta_org;
            for (double& value : best.residual) {
                value = -value;
            }
            labels_swapped = true;
        }
    }

    EquilibriumResultNative result;
    result.backend = "electrolyte_lle";
    result.problem_kind = "electrolyte_lle_flash";
    result.stable = false;
    result.split_detected = true;
    result.phases.push_back(phase_from_state("aq", best.aq_comp, 1.0 - best.beta_org, best.aq_state, options.include_phase_diagnostics));
    result.phases.push_back(phase_from_state("org", best.org_comp, best.beta_org, best.org_state, options.include_phase_diagnostics));
    result.diagnostics_string["phase_equilibrium_model"] = "electrolyte_lle_v5_native_charge_constrained_solve";
    result.diagnostics_string["equilibrium_route"] = "electrolyte_lle";
    result.diagnostics_string["route_reason"] = "ion-containing mixture";
    result.diagnostics_string["stability_analysis"] = "electrolyte_tpd";
    result.diagnostics_string["variable_model"] = basis.variable_model;
    result.diagnostics_string["acceptance_gate"] = "predictive_nonlinear_solve";
    result.diagnostics_string["solver_seed_name"] = best_seed;
    result.diagnostics_string["solver_method"] = "native_transformed_newton";
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_string["tpd_method"] = "native_tpd_global_search";
    result.diagnostics_string["gibbs_seed_method"] = "native_golden_section";
    result.diagnostics_string["phase_label_basis"] = phase_label_basis;
    result.diagnostics_string["density_diagnostics_mode"] = options.density_diagnostics;
    result.diagnostics_string["density_validity_gate"] = "passed";
    result.diagnostics_string["density_warm_start_source"] = mixture->last_density_diagnostics().warm_start_source;
    result.diagnostics_bool["phase_labels_swapped"] = labels_swapped;
    result.diagnostics_bool["experimental_coupled_density_lle"] = options.experimental_coupled_density_lle;
    result.diagnostics_bool["coupled_density_lle_attempted"] = options.experimental_coupled_density_lle;
    result.diagnostics_bool["density_fallback_used"] = false;
    result.diagnostics_bool["stability_checked"] = true;
    result.diagnostics_bool["stability_stable"] = stability.stable;
    result.diagnostics_int["basis_rank"] = basis.basis_rank;
    result.diagnostics_int["seed_attempt_count"] = static_cast<int>(attempts.size());
    result.diagnostics_int["density_failure_count"] = 0;
    result.diagnostics_int["iterations"] = best.iteration;
    result.diagnostics_int["repeated_stability_iterations"] = 1;
    result.diagnostics_int["requested_max_iterations"] = options.max_iterations;
    result.diagnostics_int["effective_max_iterations"] = solver_options.max_iterations;
    result.diagnostics_double["solver_residual_norm"] = best.solver_residual_norm;
    result.diagnostics_double["fugacity_residual_norm"] = best.solver_residual_norm;
    result.diagnostics_double["material_balance_error"] = best.material_error;
    result.diagnostics_double["charge_balance_error"] = best.charge_balance_error;
    result.diagnostics_double["gibbs_feed"] = best.gibbs_feed;
    result.diagnostics_double["gibbs_split"] = best.gibbs_split;
    result.diagnostics_double["gibbs_delta"] = best.gibbs_delta;
    result.diagnostics_double["phase_distance"] = best.phase_distance_value;
    result.diagnostics_double["stability_min_tpd"] = stability.min_tpd;
    result.diagnostics_vector["feed_composition"] = feed;
    result.diagnostics_vector["fugacity_residual"] = best.residual;
    result.diagnostics_string["best_noncollapsed_candidate"] = "accepted";
    result.attempt_diagnostics = attempts;
    result.density_diagnostics = density_failures;
    return result;
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
    if (options.jacobian_backend != "finite_difference") {
        throw ValueError(
            "electrolyte LLE residual jacobian does not yet have analytic/autodiff coverage; "
            "finite difference is only available when jacobian_backend='finite_difference' is requested explicitly."
        );
    }
    std::vector<double> feed = normalize_feed(raw_feed, mixture->ncomp(), options.min_composition, "electrolyte_lle");
    const std::vector<double>& charges = mixture->args().z;
    double feed_charge = composition_charge(feed, charges);
    if (std::abs(feed_charge) > 1.0e-10) {
        throw ValueError("electrolyte_lle feed must be charge neutral.");
    }
    if (options.min_composition <= 0.0) {
        throw ValueError("electrolyte_lle residual evaluation options contain invalid numerical controls.");
    }
    ElectrolyteBasisNative basis = build_electrolyte_basis_native(mixture, feed, species);
    PhaseStateNative feed_state = phase_state(mixture, t, p, feed, "liq", "feed");
    double gibbs_feed = electrolyte_gibbs_proxy(feed, feed_state);

    std::vector<double> active_variables = variables;
    if (has_variables) {
        if (active_variables.size() != basis.formula_feed.size()) {
            throw ValueError("electrolyte_lle residual variables length must match transformed basis dimension.");
        }
    } else if (has_initial_phases) {
        if (initial_aq.size() != feed.size() || initial_org.size() != feed.size()) {
            throw ValueError("initial_phases aq/org length must match mixture component count.");
        }
        if (!(initial_beta_org > 0.0 && initial_beta_org < 1.0) || !std::isfinite(initial_beta_org)) {
            throw ValueError("initial_phases phase_fraction must be > 0 and < 1.");
        }
        std::vector<double> aq_comp_initial = clip_normalize(initial_aq, options.min_composition);
        std::vector<double> org_comp_initial = clip_normalize(initial_org, options.min_composition);
        std::vector<double> aq_formula = explicit_to_formula(aq_comp_initial, basis);
        std::vector<double> org_formula = explicit_to_formula(org_comp_initial, basis);
        auto aq_expanded = formula_to_explicit(aq_formula, basis, feed.size());
        auto org_expanded = formula_to_explicit(org_formula, basis, feed.size());
        double beta_formula = initial_beta_org;
        double numerator = initial_beta_org / org_expanded.second;
        double denominator = numerator + (1.0 - initial_beta_org) / aq_expanded.second;
        if (denominator > 0.0) {
            beta_formula = std::max(1.0e-12, std::min(1.0 - 1.0e-12, numerator / denominator));
        }
        active_variables = pack_predictive_electrolyte_variables(beta_formula, org_formula);
    } else {
        auto seeds = electrolyte_formula_seeds(mixture, feed, basis, options);
        std::vector<double> org_formula = seeds.empty() ? basis.formula_feed : seeds.front().second;
        active_variables = pack_predictive_electrolyte_variables(0.5, org_formula);
    }
    for (double value : active_variables) {
        if (!std::isfinite(value)) {
            throw ValueError("electrolyte_lle residual variables must be finite.");
        }
    }

    auto residual_fn = [&](const std::vector<double>& candidate) {
        try {
            return evaluate_predictive_electrolyte_variables(
                mixture,
                t,
                p,
                feed,
                basis,
                candidate,
                options,
                gibbs_feed
            ).residual;
        } catch (const SolutionError&) {
            return predictive_infeasible_residual(candidate, basis, options);
        }
    };
    ElectrolyteCandidateNative current = evaluate_predictive_electrolyte_variables(
        mixture,
        t,
        p,
        feed,
        basis,
        active_variables,
        options,
        gibbs_feed
    );
    std::vector<double> base_residual = current.residual;
    const std::size_t nvar = active_variables.size();
    const std::size_t nres = base_residual.size();
    std::vector<double> jacobian(nres * nvar, 0.0);
    const double step = std::max(1.0e-7, 10.0 * options.min_composition);
    for (std::size_t col = 0; col < nvar; ++col) {
        std::vector<double> shifted = active_variables;
        shifted[col] += step;
        std::vector<double> value = residual_fn(shifted);
        for (std::size_t row = 0; row < nres; ++row) {
            jacobian[row * nvar + col] = (value[row] - base_residual[row]) / step;
        }
    }
    std::vector<double> gradient(nvar, 0.0);
    for (std::size_t col = 0; col < nvar; ++col) {
        for (std::size_t row = 0; row < nres; ++row) {
            gradient[col] += jacobian[row * nvar + col] * base_residual[row];
        }
    }
    double objective = 0.0;
    for (double value : base_residual) {
        objective += 0.5 * value * value;
    }

    ElectrolyteLLEResidualEvaluationNative out;
    out.variables = active_variables;
    out.lower_bounds.assign(nvar, -30.0);
    out.upper_bounds.assign(nvar, 30.0);
    out.residual = base_residual;
    out.jacobian_row_major = jacobian;
    out.jacobian_rows = static_cast<int>(nres);
    out.jacobian_cols = static_cast<int>(nvar);
    out.gradient = gradient;
    out.objective = objective;
    out.aq_composition = current.aq_comp;
    out.org_composition = current.org_comp;
    out.aq_ln_fugacity_coefficient = current.aq_state.ln_phi;
    out.org_ln_fugacity_coefficient = current.org_state.ln_phi;
    out.aq_density = current.aq_state.density;
    out.org_density = current.org_state.density;
    out.phase_fraction_org = current.beta_org;
    out.material_balance_error = current.material_error;
    out.charge_balance_error = current.charge_balance_error;
    out.phase_distance = current.phase_distance_value;
    out.gibbs_delta = current.gibbs_delta;
    out.diagnostics_string["solver_language"] = "c++";
    out.diagnostics_string["native_entrypoint"] = "_evaluate_electrolyte_lle_residual_native";
    out.diagnostics_string["phase_equilibrium_model"] = "electrolyte_lle_v5_native_charge_constrained_solve";
    out.diagnostics_string["variable_model"] = basis.variable_model;
    out.diagnostics_string["jacobian_backend"] = "finite_difference";
    out.diagnostics_string["derivative_backend_selected"] = "finite_difference";
    out.diagnostics_string["derivative_capability_path"] = "electrolyte_lle:explicit_finite_difference:transformed_formula_variables";
    out.diagnostics_string["unsupported_derivative_reason"] = "";
    out.diagnostics_string["hessian_backend"] = "gauss_newton";
    out.diagnostics_string["finite_difference_scheme"] = "forward";
    out.diagnostics_string["finite_difference_variable_space"] = "transformed_formula_variables";
    out.diagnostics_string["finite_difference_step_rule"] = "absolute_transformed_variable_step";
    out.diagnostics_string["hessian_kind"] = "approximate_least_squares_gauss_newton";
    out.diagnostics_string["hessian_structure"] = "dense_lower_triangular";
    out.diagnostics_bool["jacobian_available"] = true;
    out.diagnostics_bool["finite_difference_allowed"] = true;
    out.diagnostics_bool["explicit_finite_difference"] = true;
    out.diagnostics_bool["finite_difference_fallback_used"] = false;
    out.diagnostics_bool["hessian_available"] = true;
    out.diagnostics_bool["exact_hessian_available"] = false;
    out.diagnostics_bool["hessian_callback_available"] = true;
    out.diagnostics_bool["hessian_includes_second_residual_derivatives"] = false;
    out.diagnostics_bool["sparse_hessian_available"] = false;
    out.diagnostics_int["basis_rank"] = basis.basis_rank;
    out.diagnostics_double["solver_residual_norm"] = current.solver_residual_norm;
    out.diagnostics_double["fugacity_residual_norm"] = current.solver_residual_norm;
    out.diagnostics_double["material_balance_error"] = current.material_error;
    out.diagnostics_double["charge_balance_error"] = current.charge_balance_error;
    out.diagnostics_double["phase_distance"] = current.phase_distance_value;
    out.diagnostics_double["gibbs_delta"] = current.gibbs_delta;
    out.diagnostics_double["objective"] = objective;
    out.diagnostics_double["finite_difference_base_step"] = options.min_composition;
    out.diagnostics_double["finite_difference_effective_step"] = step;
    out.diagnostics_vector["feed_composition"] = feed;
    out.diagnostics_vector["fugacity_residual"] = base_residual;
    return out;
}
