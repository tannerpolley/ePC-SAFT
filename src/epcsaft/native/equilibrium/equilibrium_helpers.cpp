#include "equilibrium_helpers.h"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace epcsaft::native::equilibrium {

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

double composition_charge(const std::vector<double>& composition, const std::vector<double>& charges) {
    double out = 0.0;
    for (std::size_t i = 0; i < composition.size(); ++i) {
        out += composition[i] * charges[i];
    }
    return out;
}

double l2_norm(const std::vector<double>& values) {
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return std::sqrt(sum);
}

std::vector<double> damping_schedule(double damping) {
    double start = std::max(1.0e-6, std::min(1.0, damping));
    return {start, start * 0.5, start * 0.25, start * 0.1, start * 0.05, start * 0.01, start * 0.001};
}

}  // namespace epcsaft::native::equilibrium
