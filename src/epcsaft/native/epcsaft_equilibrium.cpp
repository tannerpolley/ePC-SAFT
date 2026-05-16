#include "epcsaft_equilibrium.h"
#include "equilibrium/equilibrium_helpers.h"
#include "equilibrium/residual_solver.h"

#ifdef EPCSAFT_HAS_CERES
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
#endif

#include <Eigen/Dense>

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <utility>

PhaseStateCompositionSensitivityResult phase_state_ln_fugacity_composition_sensitivity_cpp(
    double t,
    double p,
    std::vector<double> x,
    int phase,
    const add_args &cppargs
);

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
    std::vector<double> species_charges;
    int basis_rank = 0;
    std::string variable_model = "ascani_transformed_salt_pairs";
};

struct ElectrolyteTpdSearchNative {
    std::vector<StabilityTrialNative> trials;
    int multistart_count = 0;
    int polish_iterations = 0;
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

bool rachford_rice_beta(const std::vector<double>& feed, const std::vector<double>& k_values, double& beta, std::string& message);
std::pair<std::vector<double>, std::vector<double>> phase_compositions(
    const std::vector<double>& feed,
    const std::vector<double>& k_values,
    double beta,
    double min_composition
);

using namespace epcsaft::native::equilibrium;

struct EquilibriumBudgetExceeded : public std::runtime_error {
    explicit EquilibriumBudgetExceeded(std::string budget_trigger)
        : std::runtime_error(budget_trigger), trigger(std::move(budget_trigger)) {}

    std::string trigger;
};

struct ElectrolyteLLEBudget {
    explicit ElectrolyteLLEBudget(const EquilibriumOptionsNative& options)
        : timeout_seconds(options.timeout_seconds),
          max_seed_attempts(options.max_seed_attempts),
          max_density_failures(options.max_density_failures),
          max_total_objective_evaluations(options.max_total_objective_evaluations),
          start(std::chrono::steady_clock::now()) {}

    double timeout_seconds = 0.0;
    int max_seed_attempts = 0;
    int max_density_failures = 0;
    int max_total_objective_evaluations = 0;
    int objective_evaluations = 0;
    std::chrono::steady_clock::time_point start;

    double elapsed_seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    }

    void check_timeout(const std::string& trigger = "timeout_seconds") const {
        if (timeout_seconds > 0.0 && elapsed_seconds() >= timeout_seconds) {
            throw EquilibriumBudgetExceeded(trigger);
        }
    }

    void count_objective_evaluation() {
        ++objective_evaluations;
        if (max_total_objective_evaluations > 0 && objective_evaluations > max_total_objective_evaluations) {
            throw EquilibriumBudgetExceeded("max_total_objective_evaluations");
        }
        check_timeout();
    }

    void check_seed_attempt_count(int completed_seed_attempts) const {
        if (max_seed_attempts > 0 && completed_seed_attempts >= max_seed_attempts) {
            throw EquilibriumBudgetExceeded("max_seed_attempts");
        }
        check_timeout();
    }

    void check_density_failure_count(std::size_t density_failure_count) const {
        if (max_density_failures > 0 && density_failure_count >= static_cast<std::size_t>(max_density_failures)) {
            throw EquilibriumBudgetExceeded("max_density_failures");
        }
        check_timeout();
    }
};

void apply_electrolyte_lle_budget_diagnostics(
    EquilibriumResultNative& result,
    const EquilibriumOptionsNative& options,
    const ElectrolyteLLEBudget& budget,
    bool budget_exceeded,
    const std::string& budget_trigger
) {
    result.diagnostics_bool["budget_exceeded"] = budget_exceeded;
    result.diagnostics_string["budget_trigger"] = budget_trigger;
    result.diagnostics_double["elapsed_seconds"] = budget.elapsed_seconds();
    result.diagnostics_double["requested_timeout_seconds"] = options.timeout_seconds;
    result.diagnostics_int["objective_evaluation_count"] = budget.objective_evaluations;
    result.diagnostics_int["max_seed_attempts"] = options.max_seed_attempts;
    result.diagnostics_int["max_density_failures"] = options.max_density_failures;
    result.diagnostics_int["max_total_objective_evaluations"] = options.max_total_objective_evaluations;
    if (budget_exceeded) {
        if (result.split_detected) {
            result.diagnostics_string["message"] =
                "electrolyte LLE flash accepted a split after a configured exploratory work budget was exhausted";
        } else {
            result.diagnostics_string["message"] =
                "electrolyte LLE flash stopped after a configured work budget was exhausted";
            result.diagnostics_string["acceptance_gate"] = "predictive_budget_exhausted";
        }
    }
}

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

StabilityTrialNative fixed_composition_tpd_trial_from_state(
    const std::vector<double>& feed,
    const std::vector<double>& parent_ln_phi,
    const std::string& parent_phase,
    const std::string& trial_phase,
    const std::string& seed_name,
    const std::vector<double>& trial_composition,
    const std::vector<double>& trial_ln_phi,
    const EquilibriumOptionsNative& options,
    double threshold
) {
    std::vector<double> composition = clip_normalize(trial_composition, options.min_composition);
    double candidate_tpd = tpd_value(composition, trial_ln_phi, feed, parent_ln_phi);

    StabilityTrialNative out;
    out.parent_phase = parent_phase;
    out.trial_phase = trial_phase;
    out.seed_name = seed_name;
    out.composition = composition;
    out.tpd = candidate_tpd;
    out.iterations = 1;
    out.converged = true;
    out.unstable = candidate_tpd < -threshold;
    out.diagnostics_double["tpd_threshold"] = threshold;
    out.diagnostics_double["final_max_composition_delta"] = 0.0;
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

StabilityResultNative neutral_stability_result_from_trials(
    std::vector<StabilityTrialNative> trials,
    double threshold,
    bool fast_exit_used,
    const std::string& tpd_method
) {
    if (trials.empty()) {
        throw ValueError("stability analysis requires at least one parent and trial phase.");
    }
    auto best_iter = std::min_element(trials.begin(), trials.end(), [](const auto& a, const auto& b) {
        return a.tpd < b.tpd;
    });
    StabilityResultNative result;
    result.backend = "neutral_tpd";
    result.problem_kind = "stability";
    result.trials = std::move(trials);
    result.stable = best_iter->tpd >= -threshold;
    result.min_tpd = best_iter->tpd;
    result.parent_phase = best_iter->parent_phase;
    result.trial_phase = best_iter->trial_phase;
    result.trial_composition = best_iter->composition;
    result.diagnostics_string["stability_analysis"] = "neutral_tpd";
    result.diagnostics_bool["stable"] = result.stable;
    result.diagnostics_bool["fast_exit_used"] = fast_exit_used;
    result.diagnostics_double["tpd_threshold"] = threshold;
    result.diagnostics_double["min_tpd"] = result.min_tpd;
    result.diagnostics_int["trial_count"] = static_cast<int>(result.trials.size());
    result.diagnostics_string["min_seed_name"] = best_iter->seed_name;
    result.diagnostics_string["message"] = result.stable ? "no negative TPD trial found" : "unstable trial phase detected";
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_string["tpd_method"] = tpd_method;
    return result;
}

StabilityResultNative neutral_stability_native_impl(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& raw_feed,
    const EquilibriumOptionsNative& options,
    const std::vector<std::string>& parent_phases,
    const std::vector<std::string>& trial_phases,
    bool stop_on_first_unstable
) {
    if (mixture->has_ionic()) {
        throw ValueError("Neutral stability does not support ion-containing mixtures.");
    }
    std::vector<double> feed = normalize_feed(raw_feed, mixture->ncomp(), options.min_composition, "stability");
    std::vector<StabilityTrialNative> trials;
    double threshold = std::max(options.tolerance, 1.0e-8);
    const auto seeds = tpd_seeds(feed, options);
    bool fast_exit_used = false;
    for (const std::string& parent_phase : parent_phases) {
        PhaseStateNative parent = phase_state(mixture, t, p, feed, parent_phase);
        for (const std::string& trial_phase : trial_phases) {
            for (const auto& seed : seeds) {
                StabilityTrialNative trial = solve_tpd_trial(
                    mixture,
                    t,
                    p,
                    feed,
                    parent.ln_phi,
                    parent_phase,
                    trial_phase,
                    seed.first,
                    seed.second,
                    options,
                    threshold
                );
                trials.push_back(std::move(trial));
                if (stop_on_first_unstable && trials.back().unstable) {
                    fast_exit_used = true;
                    goto stability_done;
                }
            }
        }
    }

stability_done:
    return neutral_stability_result_from_trials(
        std::move(trials),
        threshold,
        fast_exit_used,
        fast_exit_used ? "native_tpd_fixed_point_fast_exit" : "native_tpd_fixed_point"
    );
}

StabilityResultNative neutral_split_stability_precheck_from_states(
    const std::vector<double>& raw_feed,
    const PhaseStateNative& parent_liq,
    const PhaseStateNative& parent_vap,
    const std::vector<double>& liquid_composition,
    const std::vector<double>& vapor_composition,
    const PhaseStateNative& liquid_state,
    const PhaseStateNative& vapor_state,
    const EquilibriumOptionsNative& options
) {
    std::vector<double> feed = raw_feed;
    double threshold = std::max(options.tolerance, 1.0e-8);
    std::vector<StabilityTrialNative> trials;

    trials.push_back(fixed_composition_tpd_trial_from_state(
        feed,
        parent_liq.ln_phi,
        "liq",
        "vap",
        "split_vapor_seed",
        vapor_composition,
        vapor_state.ln_phi,
        options,
        threshold
    ));
    if (!trials.back().unstable) {
        trials.push_back(fixed_composition_tpd_trial_from_state(
            feed,
            parent_vap.ln_phi,
            "vap",
            "liq",
            "split_liquid_seed",
            liquid_composition,
            liquid_state.ln_phi,
            options,
            threshold
        ));
    }
    return neutral_stability_result_from_trials(
        std::move(trials),
        threshold,
        true,
        "native_tpd_split_seed_check"
    );
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

std::vector<double> indices_to_double_vector(const std::vector<int>& indices) {
    std::vector<double> out;
    out.reserve(indices.size());
    for (int index : indices) {
        out.push_back(static_cast<double>(index));
    }
    return out;
}

std::vector<int> charged_indices_from_basis(const ElectrolyteBasisNative& basis) {
    std::vector<int> out = basis.cation_indices;
    out.insert(out.end(), basis.anion_indices.begin(), basis.anion_indices.end());
    return out;
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

std::vector<double> salt_pair_field_vector(const ElectrolyteBasisNative& basis, const std::string& field) {
    std::vector<double> out;
    out.reserve(basis.salt_pairs.size());
    for (const auto& pair : basis.salt_pairs) {
        if (field == "cation") {
            out.push_back(static_cast<double>(pair.cation));
        } else if (field == "anion") {
            out.push_back(static_cast<double>(pair.anion));
        } else if (field == "cation_stoich") {
            out.push_back(static_cast<double>(pair.cation_stoich));
        } else if (field == "anion_stoich") {
            out.push_back(static_cast<double>(pair.anion_stoich));
        }
    }
    return out;
}

void add_electrolyte_basis_diagnostics(
    EquilibriumResultNative& result,
    const ElectrolyteBasisNative& basis,
    std::size_t ncomp
) {
    result.diagnostics_string["basis_model"] = "charge_neutral_salt_pair_coordinates";
    result.diagnostics_string["basis_vector_model"] = "salt_pair_stoichiometry_rows_by_public_species";
    result.diagnostics_int["basis_rank"] = basis.basis_rank;
    result.diagnostics_int["explicit_species_count"] = static_cast<int>(ncomp);
    result.diagnostics_int["neutral_species_count"] = static_cast<int>(basis.neutral_indices.size());
    result.diagnostics_int["charged_species_count"] = static_cast<int>(basis.cation_indices.size() + basis.anion_indices.size());
    result.diagnostics_int["cation_species_count"] = static_cast<int>(basis.cation_indices.size());
    result.diagnostics_int["anion_species_count"] = static_cast<int>(basis.anion_indices.size());
    result.diagnostics_int["salt_pair_count"] = static_cast<int>(basis.salt_pairs.size());
    result.diagnostics_int["formula_variable_count"] = static_cast<int>(basis.formula_feed.size());
    result.diagnostics_int["transformed_variable_count"] = static_cast<int>(basis.formula_feed.size());
    result.diagnostics_int["basis_vector_rows"] = static_cast<int>(basis.salt_pairs.size());
    result.diagnostics_int["basis_vector_cols"] = static_cast<int>(ncomp);
    result.diagnostics_bool["phase_charge_enforced_by_basis"] = true;
    result.diagnostics_bool["material_balance_enforced_by_formula_transform"] = true;
    result.diagnostics_bool["formula_phase_positivity_enforced_by_transform"] = true;
    result.diagnostics_bool["explicit_public_species_reported"] = true;
    result.diagnostics_vector["neutral_species_indices"] = indices_to_double_vector(basis.neutral_indices);
    result.diagnostics_vector["cation_species_indices"] = indices_to_double_vector(basis.cation_indices);
    result.diagnostics_vector["anion_species_indices"] = indices_to_double_vector(basis.anion_indices);
    result.diagnostics_vector["charged_species_indices"] = indices_to_double_vector(charged_indices_from_basis(basis));
    result.diagnostics_vector["species_charge_vector"] = basis.species_charges;
    result.diagnostics_vector["formula_feed"] = basis.formula_feed;
    result.diagnostics_vector["salt_pair_cation_indices"] = salt_pair_field_vector(basis, "cation");
    result.diagnostics_vector["salt_pair_anion_indices"] = salt_pair_field_vector(basis, "anion");
    result.diagnostics_vector["salt_pair_cation_stoich"] = salt_pair_field_vector(basis, "cation_stoich");
    result.diagnostics_vector["salt_pair_anion_stoich"] = salt_pair_field_vector(basis, "anion_stoich");
    result.diagnostics_vector["basis_vectors_row_major"] = electrolyte_basis_vectors_row_major(basis, ncomp);
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
    auto fast_exit = stability.diagnostics_bool.find("fast_exit_used");
    result.diagnostics_bool["stability_fast_exit"] = fast_exit != stability.diagnostics_bool.end() && fast_exit->second;
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

std::vector<double> composition_logit_jacobian(const std::vector<double>& composition) {
    const std::size_t ncomp = composition.size();
    const std::size_t nvar = ncomp - 1;
    std::vector<double> jacobian(ncomp * nvar, 0.0);
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t var = 0; var < nvar; ++var) {
            const double delta = species == var ? 1.0 : 0.0;
            jacobian[species * nvar + var] = composition[species] * (delta - composition[var]);
        }
    }
    return jacobian;
}

std::vector<double> lle_residual_jacobian_row_major(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const LLECandidateNative& candidate
) {
    const std::size_t ncomp = feed.size();
    const std::size_t nlogit = ncomp - 1;
    const std::size_t nvar = 1 + 2 * nlogit;
    const std::size_t rows = 2 * ncomp;
    std::vector<double> jacobian(rows * nvar, 0.0);
    std::vector<double> comp1_logit_jac = composition_logit_jacobian(candidate.comp1);
    std::vector<double> comp2_logit_jac = composition_logit_jacobian(candidate.comp2);
    PhaseStateCompositionSensitivityResult sensitivity1 =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, candidate.comp1, phase_token_to_int("liq"), mixture->args());
    PhaseStateCompositionSensitivityResult sensitivity2 =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, candidate.comp2, phase_token_to_int("liq"), mixture->args());
    if (!sensitivity1.supported || !sensitivity2.supported) {
        throw SolutionError("Neutral LLE residual Jacobian requires supported phase-state fugacity composition sensitivities.");
    }

    const double dbeta_deta = candidate.beta * (1.0 - candidate.beta);
    for (std::size_t species = 0; species < ncomp; ++species) {
        const std::size_t phase_row = species;
        const std::size_t material_row = ncomp + species;
        jacobian[material_row * nvar] = dbeta_deta * (candidate.comp2[species] - candidate.comp1[species]);
        for (std::size_t var = 0; var < nlogit; ++var) {
            const std::size_t var1 = 1 + var;
            const std::size_t var2 = 1 + nlogit + var;
            double dlogf1 = phase_log_fugacity_derivative_for_species(
                candidate.comp1,
                comp1_logit_jac,
                sensitivity1.jacobian_row_major,
                species,
                var,
                ncomp,
                nlogit
            );
            double dlogf2 = phase_log_fugacity_derivative_for_species(
                candidate.comp2,
                comp2_logit_jac,
                sensitivity2.jacobian_row_major,
                species,
                var,
                ncomp,
                nlogit
            );
            jacobian[phase_row * nvar + var1] = -dlogf1;
            jacobian[phase_row * nvar + var2] = dlogf2;
            jacobian[material_row * nvar + var1] =
                (1.0 - candidate.beta) * comp1_logit_jac[species * nlogit + var];
            jacobian[material_row * nvar + var2] =
                candidate.beta * comp2_logit_jac[species * nlogit + var];
        }
    }
    return jacobian;
}

class NeutralLLEResidualProblem final : public NativeResidualProblem {
public:
    NeutralLLEResidualProblem(
        std::shared_ptr<ePCSAFTMixtureNative> mixture,
        double t,
        double p,
        std::vector<double> feed,
        EquilibriumOptionsNative options
    )
        : mixture_(std::move(mixture)),
          t_(t),
          p_(p),
          feed_(std::move(feed)),
          options_(std::move(options)) {}

    int variable_count() const override {
        return static_cast<int>(1 + 2 * (feed_.size() - 1));
    }

    int residual_count() const override {
        return static_cast<int>(2 * feed_.size());
    }

    NativeResidualEvaluation evaluate(const std::vector<double>& variables, bool need_jacobian) const override {
        LLECandidateNative candidate = evaluate_lle_variables(mixture_, t_, p_, feed_, variables, options_);
        NativeResidualEvaluation out;
        out.success = true;
        out.residual = candidate.residual;
        out.rows = static_cast<int>(candidate.residual.size());
        out.cols = static_cast<int>(variables.size());
        out.diagnostics_double["fugacity_residual_norm"] = candidate.fugacity_residual_norm;
        out.diagnostics_double["material_balance_norm"] = candidate.material_error;
        out.diagnostics_double["phase_distance"] = phase_distance(candidate.comp1, candidate.comp2);
        out.diagnostics_string["residual_surface"] = "native_neutral_lle_transformed_variables";
        out.diagnostics_string["jacobian_backend"] = "cppad_implicit";
        out.diagnostics_bool["jacobian_available"] = need_jacobian;
        if (need_jacobian) {
            out.jacobian_row_major = lle_residual_jacobian_row_major(mixture_, t_, p_, feed_, candidate);
        }
        return out;
    }

private:
    std::shared_ptr<ePCSAFTMixtureNative> mixture_;
    double t_ = 0.0;
    double p_ = 0.0;
    std::vector<double> feed_;
    EquilibriumOptionsNative options_;
};

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
    NativeResidualSolveOptions residual_options;
    residual_options.max_iterations = options.max_iterations;
    residual_options.residual_tolerance = options.tolerance;
    residual_options.function_tolerance = std::min(1.0e-12, std::max(1.0e-18, options.tolerance * options.tolerance));
    residual_options.gradient_tolerance = std::min(1.0e-10, std::max(1.0e-14, options.tolerance));
    residual_options.parameter_tolerance = std::min(1.0e-10, std::max(1.0e-14, options.tolerance));
    NeutralLLEResidualProblem residual_problem(mixture, t, p, feed, options);
    NativeResidualSolveResult solved = solve_native_residual_problem(residual_problem, variables, residual_options);
    LLECandidateNative best = evaluate_lle_variables(mixture, t, p, feed, solved.variables, options);
    best.iteration = solved.iterations;
    best.objective = solved.residual_norm_l2;
    attempt.candidate = best;
    if (solved.accepted && lle_converged(best, options)) {
        attempt.status = "converged";
        return attempt;
    }
    if (lle_degenerate(best, options)) {
        attempt.status = "degenerate";
        attempt.message = "no V2 LLE split found; best candidate collapsed to one liquid phase";
    } else {
        attempt.status = "failed";
        attempt.message = solved.rejection_reason.empty() ? "native residual solve was not accepted" : solved.rejection_reason;
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
    result.diagnostics_string["acceptance_gate"] = "residual_and_physical_gates";
    result.diagnostics_string["solver_backend"] = "ceres";
    result.diagnostics_string["selected_solver_backend"] = "ceres";
    result.diagnostics_string["solver_method"] = "ceres_trust_region_residual_solve";
    result.diagnostics_string["nonlinear_solver"] = "ceres_trust_region_residual_solve";
    result.diagnostics_string["solver_attempted"] = "ceres";
    result.diagnostics_string["solver_attempt_result"] = "accepted";
    result.diagnostics_string["accepted_solver_backend"] = "ceres";
    result.diagnostics_string["accepted_solver_method"] = "ceres_trust_region_residual_solve";
    result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
    result.diagnostics_string["jacobian_backend"] = "cppad_implicit";
    result.diagnostics_string["derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["accepted_derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["derivative_status"] = "residual_jacobian_available";
    result.diagnostics_bool["solution_accepted"] = true;
    result.diagnostics_bool["derivative_available"] = true;
    result.diagnostics_bool["jacobian_available"] = true;
    result.diagnostics_bool["jacobian_available_for_accepted_state"] = true;
    result.diagnostics_bool["derivative_available_for_accepted_state"] = true;
    result.diagnostics_bool["hessian_available"] = false;
    result.diagnostics_string["anti_trivial_solution_strategy"] = "phase_fraction_and_phase_distance_gate";
    result.diagnostics_string["association_solver_status"] = "not_active";
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
    result.diagnostics_string["acceptance_gate"] = "neutral_lle_not_accepted";
    result.diagnostics_string["solver_backend"] = "none";
    result.diagnostics_string["selected_solver_backend"] = "none";
    result.diagnostics_string["solver_method"] = "none";
    result.diagnostics_string["nonlinear_solver"] = "none";
    result.diagnostics_string["solver_attempted"] = attempt.candidate.residual.empty() ? "none" : "ceres";
    result.diagnostics_string["solver_attempt_result"] = "failed";
    result.diagnostics_string["accepted_solver_backend"] = "none";
    result.diagnostics_string["accepted_solver_method"] = "none";
    result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
    result.diagnostics_string["derivative_backend"] = "none";
    result.diagnostics_string["derivative_status"] = "no_accepted_state";
    result.diagnostics_string["accepted_derivative_backend"] = "none";
    result.diagnostics_string["attempted_jacobian_backend"] = attempt.candidate.residual.empty() ? "none" : "cppad_implicit";
    result.diagnostics_string["attempted_derivative_backend"] = attempt.candidate.residual.empty() ? "none" : "cppad_implicit";
    result.diagnostics_bool["solution_accepted"] = false;
    result.diagnostics_bool["derivative_available"] = false;
    result.diagnostics_bool["jacobian_available"] = false;
    result.diagnostics_bool["attempted_jacobian_available"] = !attempt.candidate.residual.empty();
    result.diagnostics_bool["attempted_derivative_available"] = !attempt.candidate.residual.empty();
    result.diagnostics_bool["jacobian_available_for_accepted_state"] = false;
    result.diagnostics_bool["derivative_available_for_accepted_state"] = false;
    result.diagnostics_bool["hessian_available"] = false;
    result.diagnostics_string["anti_trivial_solution_strategy"] = "phase_fraction_and_phase_distance_gate";
    result.diagnostics_string["association_solver_status"] = "not_active";
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
    return neutral_stability_native_impl(mixture, t, p, raw_feed, options, parent_phases, trial_phases, false);
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
        result.split_detected = false;
        if (options.stability_precheck) {
            EquilibriumOptionsNative stab_opts = precheck_options(options);
            StabilityResultNative stability = neutral_stability_native_impl(
                mixture,
                t,
                p,
                feed,
                stab_opts,
                {"liq", "vap"},
                {"liq", "vap"},
                false
            );
            merge_stability_diagnostics(result, stability, stab_opts);
            result.stable = stability.stable;
        } else {
            skipped_stability_diagnostics(result);
            result.stable = false;
        }
        result.diagnostics_int["iterations"] = 0;
        result.diagnostics_double["fugacity_residual_norm"] = 0.0;
        result.diagnostics_double["material_balance_error"] = 0.0;
        result.diagnostics_double["vapor_fraction"] = beta;
        result.diagnostics_string["message"] = no_split_message;
        result.diagnostics_bool["point_solver_split_detected"] = false;
        result.diagnostics_string["point_solver_message"] = no_split_message;
        result.diagnostics_string["solver_language"] = "c++";
        result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
        result.diagnostics_bool["neutral_fast_path"] = true;
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
            if (options.stability_precheck) {
                EquilibriumOptionsNative stab_opts = precheck_options(options);
                StabilityResultNative stability = neutral_split_stability_precheck_from_states(
                    feed,
                    liquid_seed,
                    vapor_seed,
                    best_x,
                    best_y,
                    best_liquid,
                    best_vapor,
                    stab_opts
                );
                if (stability.stable) {
                    stability = neutral_stability_native_impl(
                        mixture,
                        t,
                        p,
                        feed,
                        stab_opts,
                        {"liq", "vap"},
                        {"liq", "vap"},
                        false
                    );
                }
                merge_stability_diagnostics(result, stability, stab_opts);
            } else {
                skipped_stability_diagnostics(result);
            }
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
            result.diagnostics_bool["neutral_fast_path"] = true;
            return result;
        }
        for (std::size_t i = 0; i < feed.size(); ++i) {
            double ln_k_target = liquid.ln_phi[i] - vapor.ln_phi[i];
            ln_k[i] = (1.0 - options.damping) * ln_k[i] + options.damping * ln_k_target;
        }
    }
    throw SolutionError("neutral TP flash did not converge after native iterations.");
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
        stability = neutral_stability_native_impl(
            mixture,
            t,
            p,
            feed,
            stab_opts,
            {"liq"},
            {"liq"},
            has_initial_phases
        );
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

double residual_slice_max_abs(const std::vector<double>& residual, std::size_t begin, std::size_t end) {
    double out = 0.0;
    end = std::min(end, residual.size());
    for (std::size_t i = begin; i < end; ++i) {
        out = std::max(out, std::abs(residual[i]));
    }
    return out;
}

void add_electrolyte_candidate_state_diagnostics(
    EquilibriumResultNative& result,
    const ElectrolyteBasisNative& basis,
    const std::vector<double>& feed,
    const ElectrolyteCandidateNative& candidate,
    bool has_candidate,
    bool accepted_candidate
) {
    result.diagnostics_bool["transformed_variable_candidate_available"] = has_candidate;
    result.diagnostics_bool["accepted_transformed_variables_feasible"] = accepted_candidate;
    result.diagnostics_double["phase_charge_balance_feed"] = composition_charge(feed, basis.species_charges);
    if (!has_candidate) {
        return;
    }
    const double aq_charge = composition_charge(candidate.aq_comp, basis.species_charges);
    const double org_charge = composition_charge(candidate.org_comp, basis.species_charges);
    result.diagnostics_double["phase_charge_balance_aq"] = aq_charge;
    result.diagnostics_double["phase_charge_balance_org"] = org_charge;
    result.diagnostics_double["phase_charge_balance_max_abs"] = std::max(
        std::abs(result.diagnostics_double["phase_charge_balance_feed"]),
        std::max(std::abs(aq_charge), std::abs(org_charge))
    );
    const std::size_t neutral_rows = basis.neutral_indices.size();
    const std::size_t ionic_rows = basis.salt_pairs.size();
    result.diagnostics_double["neutral_fugacity_residual_norm"] =
        residual_slice_max_abs(candidate.residual, 0, neutral_rows);
    result.diagnostics_double["ionic_equilibrium_residual_norm"] =
        residual_slice_max_abs(candidate.residual, neutral_rows, neutral_rows + ionic_rows);
    result.diagnostics_double["phase_equilibrium_residual_norm"] = candidate.solver_residual_norm;
    result.diagnostics_double["material_balance_norm"] = candidate.material_error;
    result.diagnostics_double["phase_charge_balance_norm"] =
        result.diagnostics_double["phase_charge_balance_max_abs"];
    result.diagnostics_double["scaled_solver_residual_norm"] = candidate.solver_residual_norm;
    result.diagnostics_double["unscaled_solver_residual_norm"] = candidate.solver_residual_norm;
    const std::string prefix = accepted_candidate ? "accepted" : "best_candidate";
    result.diagnostics_double[prefix + "_beta_formula"] = candidate.beta_formula;
    result.diagnostics_double[prefix + "_beta_org"] = candidate.beta_org;
    result.diagnostics_vector[prefix + "_aq_formula"] = candidate.aq_formula;
    result.diagnostics_vector[prefix + "_org_formula"] = candidate.org_formula;
    result.diagnostics_vector["material_balance_residual"] = candidate.material_residual;
}

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

std::vector<double> electrolyte_full_residual_vector(const ElectrolyteCandidateNative& candidate) {
    std::vector<double> residual = candidate.residual;
    residual.insert(residual.end(), candidate.material_residual.begin(), candidate.material_residual.end());
    return residual;
}

double electrolyte_residual_cost(const std::vector<double>& residual) {
    double cost = 0.0;
    for (double value : residual) {
        cost += 0.5 * value * value;
    }
    return cost;
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
    result.diagnostics_string["solver_attempted"] = "ceres";
    result.diagnostics_string["solver_attempt_result"] = "failed";
    result.diagnostics_string["solver_backend"] = "none";
    result.diagnostics_string["selected_solver_backend"] = "none";
    result.diagnostics_string["accepted_solver_backend"] = "none";
    result.diagnostics_string["accepted_solver_method"] = "none";
    result.diagnostics_string["solver_seed_name"] = best_seed;
    result.diagnostics_string["solver_method"] = "none";
    result.diagnostics_string["attempted_solver_method"] = "ceres_trust_region_residual_solve";
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_string["jacobian_backend"] = "none";
    result.diagnostics_string["derivative_backend"] = "none";
    result.diagnostics_string["accepted_derivative_backend"] = "none";
    result.diagnostics_string["attempted_jacobian_backend"] = "cppad_implicit";
    result.diagnostics_string["attempted_derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["residual_surface_jacobian_backend"] = "cppad_implicit";
    result.diagnostics_string["residual_surface_derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["jacobian_scope"] = "transformed_variables_phase_state_implicit_density";
    result.diagnostics_string["tpd_method"] = "native_tpd_global_search";
    result.diagnostics_string["gibbs_seed_method"] = "native_golden_section";
    result.diagnostics_string["best_failure_reason"] = reason;
    result.diagnostics_string["message"] = "electrolyte LLE flash did not converge";
    result.diagnostics_string["density_diagnostics_mode"] = options.density_diagnostics;
    result.diagnostics_string["density_validity_gate"] = "not_evaluated";
    result.diagnostics_string["density_warm_start_source"] = "";
    result.diagnostics_bool["experimental_coupled_density_lle"] = options.experimental_coupled_density_lle;
    result.diagnostics_bool["coupled_density_lle_attempted"] = options.experimental_coupled_density_lle;
    result.diagnostics_bool["density_best_candidate_refinement_used"] = false;
    result.diagnostics_bool["return_best_effort"] = options.return_best_effort;
    result.diagnostics_int["density_failure_count"] = 0;
    result.diagnostics_bool["stability_checked"] = true;
    result.diagnostics_bool["stability_stable"] = stability.stable;
    result.diagnostics_bool["unstable_feed_collapsed_all_candidates"] = stability.min_tpd < -std::max(options.tolerance, 1.0e-8)
        && reason == "candidate collapsed to one phase";
    result.diagnostics_bool["solution_accepted"] = false;
    result.diagnostics_bool["ceres_accepted_solve"] = false;
    result.diagnostics_bool["jacobian_available"] = false;
    result.diagnostics_bool["derivative_available"] = false;
    result.diagnostics_bool["attempted_jacobian_available"] = has_best;
    result.diagnostics_bool["attempted_derivative_available"] = has_best;
    result.diagnostics_bool["jacobian_available_for_accepted_state"] = false;
    result.diagnostics_bool["derivative_available_for_accepted_state"] = false;
    add_electrolyte_basis_diagnostics(result, basis, feed.size());
    add_electrolyte_candidate_state_diagnostics(result, basis, feed, best, has_best, false);
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
    bool has_best_noncollapsed = has_best && best.phase_distance_value > std::max(1.0e-4, split_distance_tolerance(options));
    if (has_best_noncollapsed) {
        result.diagnostics_vector["best_noncollapsed_candidate_aq"] = best.aq_comp;
        result.diagnostics_vector["best_noncollapsed_candidate_org"] = best.org_comp;
        result.diagnostics_string["best_noncollapsed_candidate"] = "available";
    } else {
        result.diagnostics_string["best_noncollapsed_candidate"] = "none";
    }
    result.diagnostics_bool["best_effort_phases_returned"] = false;
    if (options.return_best_effort && has_best_noncollapsed && best.aq_comp.size() == feed.size() && best.org_comp.size() == feed.size()) {
        result.phases.push_back(phase_from_state("aq", best.aq_comp, 1.0 - best.beta_org, best.aq_state, options.include_phase_diagnostics));
        result.phases.push_back(phase_from_state("org", best.org_comp, best.beta_org, best.org_state, options.include_phase_diagnostics));
        result.diagnostics_bool["best_effort_phases_returned"] = true;
        result.diagnostics_string["message"] = "electrolyte LLE flash returned best-effort diagnostic phases without predictive acceptance";
    }
    result.attempt_diagnostics = attempts;
    return result;
}

std::vector<double> gibbs_seed_variables_from_trial(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& trial_composition,
    ElectrolyteLLEBudget* budget
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
        if (budget != nullptr) {
            budget->count_objective_evaluation();
        }
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
        if (budget != nullptr) {
            budget->check_timeout();
        }
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

struct ElectrolyteCeresAttemptNative {
    ElectrolyteCandidateNative candidate;
    std::vector<double> variables;
    bool solution_usable = false;
    bool budget_exceeded = false;
    std::string budget_trigger;
    std::string termination_type;
    std::string message;
    std::string trust_region_strategy = "ceres_internal_trust_region";
    std::string linear_solver = "dense_qr";
    double initial_cost = std::numeric_limits<double>::infinity();
    double final_cost = std::numeric_limits<double>::infinity();
    int status = 0;
    int iterations = 0;
    int residual_evaluations = 0;
    int jacobian_evaluations = 0;
};

#ifdef EPCSAFT_HAS_CERES
std::string ceres_termination_type_name(ceres::TerminationType type) {
    switch (type) {
        case ceres::CONVERGENCE:
            return "convergence";
        case ceres::NO_CONVERGENCE:
            return "no_convergence";
        case ceres::FAILURE:
            return "failure";
        case ceres::USER_SUCCESS:
            return "user_success";
        case ceres::USER_FAILURE:
            return "user_failure";
        default:
            return std::to_string(static_cast<int>(type));
    }
}

class ElectrolyteLLEResidualCeresCostFunction final : public ceres::CostFunction {
public:
    ElectrolyteLLEResidualCeresCostFunction(
        std::shared_ptr<ePCSAFTMixtureNative> mixture,
        double t,
        double p,
        std::vector<double> feed,
        ElectrolyteBasisNative basis,
        EquilibriumOptionsNative options,
        double gibbs_feed,
        std::vector<DensitySolveDiagnostics>* density_failures,
        ElectrolyteLLEBudget* budget
    )
        : mixture_(std::move(mixture)),
          t_(t),
          p_(p),
          feed_(std::move(feed)),
          basis_(std::move(basis)),
          options_(std::move(options)),
          gibbs_feed_(gibbs_feed),
          density_failures_(density_failures),
          budget_(budget)
    {
        const std::size_t residual_count = basis_.neutral_indices.size() + basis_.salt_pairs.size() + feed_.size();
        set_num_residuals(static_cast<int>(residual_count));
        mutable_parameter_block_sizes()->push_back(static_cast<int>(basis_.formula_feed.size()));
    }

    bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override {
        try {
            if (budget_ != nullptr) {
                budget_->count_objective_evaluation();
            }
            std::vector<double> variables(parameters[0], parameters[0] + basis_.formula_feed.size());
            ElectrolyteCandidateNative candidate = evaluate_predictive_electrolyte_variables(
                mixture_,
                t_,
                p_,
                feed_,
                basis_,
                variables,
                options_,
                gibbs_feed_
            );
            std::vector<double> residual = electrolyte_full_residual_vector(candidate);
            for (std::size_t row = 0; row < residual.size(); ++row) {
                residuals[row] = residual[row];
            }
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                std::vector<double> jacobian = electrolyte_residual_jacobian_row_major(
                    mixture_,
                    t_,
                    p_,
                    feed_,
                    basis_,
                    candidate
                );
                if (jacobian.size() != residual.size() * variables.size()) {
                    throw SolutionError("Electrolyte residual Jacobian shape does not match Ceres residual block.");
                }
                for (std::size_t i = 0; i < jacobian.size(); ++i) {
                    jacobians[0][i] = jacobian[i];
                }
            }
            return true;
        } catch (const EquilibriumBudgetExceeded& exc) {
            budget_exceeded_ = true;
            budget_trigger_ = exc.trigger;
            return false;
        } catch (const SolutionError&) {
            append_last_density_failure(mixture_, density_failures_);
            try {
                if (budget_ != nullptr && density_failures_ != nullptr) {
                    budget_->check_density_failure_count(density_failures_->size());
                }
            } catch (const EquilibriumBudgetExceeded& exc) {
                budget_exceeded_ = true;
                budget_trigger_ = exc.trigger;
            }
            return false;
        } catch (const std::exception&) {
            return false;
        }
    }

    bool budget_exceeded() const {
        return budget_exceeded_;
    }

    const std::string& budget_trigger() const {
        return budget_trigger_;
    }

private:
    std::shared_ptr<ePCSAFTMixtureNative> mixture_;
    double t_ = 0.0;
    double p_ = 0.0;
    std::vector<double> feed_;
    ElectrolyteBasisNative basis_;
    EquilibriumOptionsNative options_;
    double gibbs_feed_ = 0.0;
    std::vector<DensitySolveDiagnostics>* density_failures_ = nullptr;
    ElectrolyteLLEBudget* budget_ = nullptr;
    mutable bool budget_exceeded_ = false;
    mutable std::string budget_trigger_;
};
#endif

ElectrolyteCeresAttemptNative solve_predictive_electrolyte_ceres_attempt(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const ElectrolyteBasisNative& basis,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& seed_variables,
    double gibbs_feed,
    std::vector<DensitySolveDiagnostics>* density_failures,
    ElectrolyteLLEBudget* budget,
    bool* budget_exceeded,
    std::string* budget_trigger
) {
    ElectrolyteCeresAttemptNative out;
    out.variables = seed_variables;
    try {
        if (budget != nullptr) {
            budget->check_timeout();
        }
#ifndef EPCSAFT_HAS_CERES
        throw SolutionError("Ceres support is required for electrolyte LLE residual solving.");
#else
        ElectrolyteCandidateNative initial = evaluate_predictive_electrolyte_variables(
            mixture,
            t,
            p,
            feed,
            basis,
            out.variables,
            options,
            gibbs_feed
        );
        out.initial_cost = electrolyte_residual_cost(electrolyte_full_residual_vector(initial));

        ceres::Problem problem;
        auto *cost = new ElectrolyteLLEResidualCeresCostFunction(
            mixture,
            t,
            p,
            feed,
            basis,
            options,
            gibbs_feed,
            density_failures,
            budget
        );
        problem.AddResidualBlock(cost, nullptr, out.variables.data());
        for (int index = 0; index < static_cast<int>(out.variables.size()); ++index) {
            problem.SetParameterLowerBound(out.variables.data(), index, -100.0);
            problem.SetParameterUpperBound(out.variables.data(), index, 100.0);
        }

        ceres::Solver::Options ceres_options;
        ceres_options.linear_solver_type = ceres::DENSE_QR;
        ceres_options.max_num_iterations = options.max_iterations;
        ceres_options.minimizer_progress_to_stdout = false;
        ceres_options.logging_type = ceres::SILENT;
        ceres_options.function_tolerance = std::min(1.0e-12, std::max(1.0e-18, options.tolerance * options.tolerance));
        ceres_options.gradient_tolerance = std::min(1.0e-10, std::max(1.0e-14, options.tolerance));
        ceres_options.parameter_tolerance = std::min(1.0e-10, std::max(1.0e-14, options.tolerance));
        if (budget != nullptr && budget->timeout_seconds > 0.0) {
            double remaining = budget->timeout_seconds - budget->elapsed_seconds();
            if (remaining <= 0.0) {
                throw EquilibriumBudgetExceeded("timeout_seconds");
            }
            ceres_options.max_solver_time_in_seconds = remaining;
        }

        ceres::Solver::Summary summary;
        ceres::Solve(ceres_options, &problem, &summary);
        if (cost->budget_exceeded()) {
            out.budget_exceeded = true;
            out.budget_trigger = cost->budget_trigger();
        }
        out.status = static_cast<int>(summary.termination_type);
        out.termination_type = ceres_termination_type_name(summary.termination_type);
        out.message = summary.BriefReport();
        out.iterations = static_cast<int>(summary.iterations.size());
        out.residual_evaluations = static_cast<int>(summary.num_residual_evaluations);
        out.jacobian_evaluations = static_cast<int>(summary.num_jacobian_evaluations);
        out.solution_usable = summary.IsSolutionUsable();
        out.candidate = evaluate_predictive_electrolyte_variables(
            mixture,
            t,
            p,
            feed,
            basis,
            out.variables,
            options,
            gibbs_feed
        );
        out.candidate.iteration = out.iterations;
        std::vector<double> final_residual = electrolyte_full_residual_vector(out.candidate);
        out.final_cost = electrolyte_residual_cost(final_residual);
        out.candidate.objective = std::sqrt(std::max(0.0, 2.0 * out.final_cost));
        return out;
#endif
    } catch (const EquilibriumBudgetExceeded& exc) {
        if (budget_exceeded != nullptr) {
            *budget_exceeded = true;
        }
        if (budget_trigger != nullptr) {
            *budget_trigger = exc.trigger;
        }
        out.budget_exceeded = true;
        out.budget_trigger = exc.trigger;
        out.message = exc.what();
    } catch (const SolutionError& exc) {
        append_last_density_failure(mixture, density_failures);
        out.message = exc.what();
    }
    out.candidate.objective = 1.0e300;
    out.candidate.solver_residual_norm = 1.0e300;
    out.candidate.material_error = 1.0e300;
    out.candidate.charge_balance_error = 1.0e300;
    out.candidate.iteration = 0;
    out.candidate.residual.assign(std::max<std::size_t>(1, basis.formula_feed.size()), 1.0e300);
    return out;
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
    ElectrolyteLLEBudget budget(solver_options);
    bool budget_exceeded = false;
    std::string budget_trigger;
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
    try {
        budget.check_timeout();
        if (!stability.trial_composition.empty() && phase_distance(stability.trial_composition, feed) > split_distance_tolerance(solver_options)) {
            seed_variables.emplace_back("native_tpd_trial", pack_predictive_electrolyte_variables(0.5, explicit_to_formula(stability.trial_composition, basis)));
            try {
                seed_variables.emplace_back("native_gibbs_tpd_trial", gibbs_seed_variables_from_trial(
                    mixture,
                    t,
                    p,
                    feed,
                    basis,
                    solver_options,
                    stability.trial_composition,
                    &budget
                ));
            } catch (const SolutionError&) {
            }
        }
    } catch (const EquilibriumBudgetExceeded& exc) {
        budget_exceeded = true;
        budget_trigger = exc.trigger;
    }
    if (!budget_exceeded) {
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
    }
    ElectrolyteCandidateNative best;
    ElectrolyteCeresAttemptNative best_ceres;
    std::string best_seed;
    bool has_best = false;
    bool has_best_ceres = false;
    bool accepted = false;
    std::vector<EquilibriumAttemptDiagnosticsNative> attempts;
    std::vector<DensitySolveDiagnostics> density_failures;
    for (const auto& seed : seed_variables) {
        try {
            budget.check_seed_attempt_count(static_cast<int>(attempts.size()));
        } catch (const EquilibriumBudgetExceeded& exc) {
            budget_exceeded = true;
            budget_trigger = exc.trigger;
            break;
        }
        bool attempt_budget_exceeded = false;
        std::string attempt_budget_trigger;
        ElectrolyteCeresAttemptNative ceres_attempt = solve_predictive_electrolyte_ceres_attempt(
            mixture,
            t,
            p,
            feed,
            basis,
            solver_options,
            seed.second,
            gibbs_feed,
            &density_failures,
            &budget,
            &attempt_budget_exceeded,
            &attempt_budget_trigger
        );
        ElectrolyteCandidateNative candidate = ceres_attempt.candidate;
        bool candidate_accepted = ceres_attempt.solution_usable
            && ceres_attempt.termination_type == "convergence"
            && predictive_electrolyte_accepted(candidate, solver_options);
        attempts.push_back(electrolyte_attempt_diagnostics(seed.first, candidate, solver_options, candidate_accepted));
        if (!has_best || candidate.objective < best.objective) {
            best = candidate;
            best_ceres = ceres_attempt;
            best_seed = seed.first;
            has_best = true;
            has_best_ceres = true;
        }
        if (candidate_accepted) {
            best = candidate;
            best_ceres = ceres_attempt;
            best_seed = seed.first;
            accepted = true;
            has_best_ceres = true;
            break;
        }
        if (attempt_budget_exceeded || ceres_attempt.budget_exceeded) {
            budget_exceeded = true;
            budget_trigger = !attempt_budget_trigger.empty() ? attempt_budget_trigger : ceres_attempt.budget_trigger;
            break;
        }
        try {
            budget.check_seed_attempt_count(static_cast<int>(attempts.size()));
        } catch (const EquilibriumBudgetExceeded& exc) {
            budget_exceeded = true;
            budget_trigger = exc.trigger;
            break;
        }
    }
    if (!has_best || !accepted) {
        EquilibriumResultNative result = electrolyte_lle_failure_result(feed, basis, stability, best, best_seed, has_best, solver_options, attempts);
        result.density_diagnostics = density_failures;
        result.diagnostics_int["density_failure_count"] = static_cast<int>(density_failures.size());
        apply_electrolyte_lle_budget_diagnostics(result, solver_options, budget, budget_exceeded, budget_trigger);
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
    result.diagnostics_string["acceptance_gate"] = "ceres_residual_solve";
    result.diagnostics_string["solver_backend"] = "ceres";
    result.diagnostics_string["selected_solver_backend"] = "ceres";
    result.diagnostics_string["solver_attempted"] = "ceres";
    result.diagnostics_string["solver_attempt_result"] = "accepted";
    result.diagnostics_string["accepted_solver_backend"] = "ceres";
    result.diagnostics_string["accepted_solver_method"] = "ceres_trust_region_residual_solve";
    result.diagnostics_string["solver_seed_name"] = best_seed;
    result.diagnostics_string["solver_method"] = "ceres_trust_region_residual_solve";
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_equilibrium_native";
    result.diagnostics_string["jacobian_backend"] = "cppad_implicit";
    result.diagnostics_string["derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["accepted_derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["residual_surface_jacobian_backend"] = "cppad_implicit";
    result.diagnostics_string["residual_surface_derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["jacobian_scope"] = "transformed_variables_phase_state_implicit_density";
    result.diagnostics_string["ceres_trust_region_strategy"] = has_best_ceres ? best_ceres.trust_region_strategy : "";
    result.diagnostics_string["ceres_linear_solver"] = has_best_ceres ? best_ceres.linear_solver : "";
    result.diagnostics_string["ceres_termination_type"] = has_best_ceres ? best_ceres.termination_type : "";
    result.diagnostics_string["ceres_summary"] = has_best_ceres ? best_ceres.message : "";
    result.diagnostics_string["tpd_method"] = "native_tpd_global_search";
    result.diagnostics_string["gibbs_seed_method"] = "native_golden_section";
    result.diagnostics_string["phase_label_basis"] = phase_label_basis;
    result.diagnostics_string["density_diagnostics_mode"] = options.density_diagnostics;
    result.diagnostics_string["density_validity_gate"] = "passed";
    result.diagnostics_string["density_warm_start_source"] = mixture->last_density_diagnostics().warm_start_source;
    result.diagnostics_bool["phase_labels_swapped"] = labels_swapped;
    result.diagnostics_bool["experimental_coupled_density_lle"] = options.experimental_coupled_density_lle;
    result.diagnostics_bool["coupled_density_lle_attempted"] = options.experimental_coupled_density_lle;
    result.diagnostics_bool["density_best_candidate_refinement_used"] = false;
    result.diagnostics_bool["return_best_effort"] = options.return_best_effort;
    result.diagnostics_bool["best_effort_phases_returned"] = false;
    result.diagnostics_bool["stability_checked"] = true;
    result.diagnostics_bool["stability_stable"] = stability.stable;
    result.diagnostics_bool["solution_accepted"] = true;
    result.diagnostics_bool["ceres_accepted_solve"] = true;
    result.diagnostics_bool["ceres_solution_usable"] = has_best_ceres ? best_ceres.solution_usable : false;
    result.diagnostics_bool["ceres_converged"] = has_best_ceres && best_ceres.termination_type == "convergence";
    result.diagnostics_bool["jacobian_available"] = true;
    result.diagnostics_bool["derivative_available"] = true;
    result.diagnostics_bool["jacobian_available_for_accepted_state"] = true;
    result.diagnostics_bool["derivative_available_for_accepted_state"] = true;
    add_electrolyte_basis_diagnostics(result, basis, feed.size());
    add_electrolyte_candidate_state_diagnostics(result, basis, feed, best, true, true);
    result.diagnostics_int["basis_rank"] = basis.basis_rank;
    result.diagnostics_int["seed_attempt_count"] = static_cast<int>(attempts.size());
    result.diagnostics_int["density_failure_count"] = static_cast<int>(density_failures.size());
    result.diagnostics_int["iterations"] = best.iteration;
    result.diagnostics_int["repeated_stability_iterations"] = 1;
    result.diagnostics_int["requested_max_iterations"] = options.max_iterations;
    result.diagnostics_int["effective_max_iterations"] = solver_options.max_iterations;
    result.diagnostics_int["ceres_status"] = has_best_ceres ? best_ceres.status : 0;
    result.diagnostics_int["ceres_iteration_count"] = has_best_ceres ? best_ceres.iterations : 0;
    result.diagnostics_int["ceres_residual_evaluation_count"] = has_best_ceres ? best_ceres.residual_evaluations : 0;
    result.diagnostics_int["ceres_jacobian_evaluation_count"] = has_best_ceres ? best_ceres.jacobian_evaluations : 0;
    result.diagnostics_double["solver_residual_norm"] = best.solver_residual_norm;
    result.diagnostics_double["fugacity_residual_norm"] = best.solver_residual_norm;
    result.diagnostics_double["material_balance_error"] = best.material_error;
    result.diagnostics_double["charge_balance_error"] = best.charge_balance_error;
    result.diagnostics_double["gibbs_feed"] = best.gibbs_feed;
    result.diagnostics_double["gibbs_split"] = best.gibbs_split;
    result.diagnostics_double["gibbs_delta"] = best.gibbs_delta;
    result.diagnostics_double["phase_distance"] = best.phase_distance_value;
    result.diagnostics_double["stability_min_tpd"] = stability.min_tpd;
    result.diagnostics_double["ceres_initial_cost"] = has_best_ceres ? best_ceres.initial_cost : 0.0;
    result.diagnostics_double["ceres_final_cost"] = has_best_ceres ? best_ceres.final_cost : 0.0;
    result.diagnostics_vector["feed_composition"] = feed;
    result.diagnostics_vector["fugacity_residual"] = best.residual;
    result.diagnostics_string["best_noncollapsed_candidate"] = "accepted";
    result.attempt_diagnostics = attempts;
    result.density_diagnostics = density_failures;
    apply_electrolyte_lle_budget_diagnostics(result, solver_options, budget, budget_exceeded, budget_trigger);
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
