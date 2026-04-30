#include "epcsaft_equilibrium.h"

#include <Eigen/Dense>

#include <algorithm>
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

int phase_token_to_int(const std::string& phase) {
    if (phase == "liq" || phase == "liquid" || phase == "aq" || phase == "org" || phase == "liq1" || phase == "liq2") {
        return 0;
    }
    if (phase == "vap" || phase == "vapor" || phase == "gas") {
        return 1;
    }
    throw ValueError("phase must be 'liq' or 'vap'.");
}

std::vector<double> clip_normalize(const std::vector<double>& composition, double min_composition) {
    std::vector<double> out(composition.size(), min_composition);
    double total = 0.0;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        out[i] = std::max(composition[i], min_composition);
        total += out[i];
    }
    if (!std::isfinite(total) || total <= 0.0) {
        throw ValueError("composition must have a positive finite sum.");
    }
    for (double& value : out) {
        value /= total;
    }
    return out;
}

std::vector<double> normalize_feed(const std::vector<double>& feed, std::size_t ncomp, double min_composition, const std::string& kind) {
    if (feed.size() != ncomp) {
        std::ostringstream msg;
        msg << "Feed composition length (" << feed.size() << ") must match mixture component count (" << ncomp << ").";
        throw ValueError(msg.str());
    }
    double total = 0.0;
    for (double value : feed) {
        if (!std::isfinite(value)) {
            throw ValueError("Feed composition z must contain only finite values.");
        }
        if (value < 0.0) {
            throw ValueError("Feed composition z must be non-negative.");
        }
        total += value;
    }
    if (total <= 0.0) {
        throw ValueError("Feed composition z must have a positive sum.");
    }
    std::vector<double> out(feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        out[i] = feed[i] / total;
        if (out[i] < min_composition) {
            throw ValueError(kind + " requires each feed composition entry to be >= min_composition.");
        }
    }
    return out;
}

double max_abs(const std::vector<double>& values) {
    double out = 0.0;
    for (double value : values) {
        out = std::max(out, std::abs(value));
    }
    return out;
}

double phase_distance(const std::vector<double>& a, const std::vector<double>& b) {
    double out = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        out = std::max(out, std::abs(a[i] - b[i]));
    }
    return out;
}

double split_distance_tolerance(const EquilibriumOptionsNative& options) {
    return std::max(1.0e-8, 100.0 * options.min_composition);
}

PhaseStateNative phase_state(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& composition,
    const std::string& phase
) {
    PhaseStateNative out;
    out.state = mixture->state(t, composition, phase_token_to_int(phase), true, p, false, 0.0);
    out.ln_phi = out.state->ln_fugacity_coefficient();
    out.density = out.state->density();
    return out;
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

std::vector<double> composition_from_log_weights(const std::vector<double>& log_weights, double min_composition) {
    double largest = *std::max_element(log_weights.begin(), log_weights.end());
    std::vector<double> weights(log_weights.size(), 0.0);
    for (std::size_t i = 0; i < log_weights.size(); ++i) {
        weights[i] = std::exp(std::max(-700.0, std::min(700.0, log_weights[i] - largest)));
    }
    return clip_normalize(weights, min_composition);
}

std::vector<double> component_rich_composition(std::size_t ncomp, std::size_t rich_index, double min_composition) {
    std::vector<double> composition(ncomp, min_composition);
    composition[rich_index] = std::max(min_composition, 1.0 - min_composition * static_cast<double>(ncomp - 1));
    return clip_normalize(composition, min_composition);
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

double l2_norm(const std::vector<double>& values) {
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return std::sqrt(sum);
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

std::vector<double> damping_schedule(double damping) {
    double start = std::max(1.0e-6, std::min(1.0, damping));
    return {start, start * 0.5, start * 0.25, start * 0.1, start * 0.05, start * 0.01, start * 0.001};
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
