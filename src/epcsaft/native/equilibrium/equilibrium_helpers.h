#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "../epcsaft_equilibrium.h"

namespace epcsaft::native::equilibrium {

int phase_token_to_int(const std::string& phase);
std::vector<double> clip_normalize(const std::vector<double>& composition, double min_composition);
std::vector<double> normalize_feed(const std::vector<double>& feed, std::size_t ncomp, double min_composition, const std::string& kind);
double max_abs(const std::vector<double>& values);
double phase_distance(const std::vector<double>& a, const std::vector<double>& b);
double split_distance_tolerance(const EquilibriumOptionsNative& options);
std::vector<double> composition_from_log_weights(const std::vector<double>& log_weights, double min_composition);
std::vector<double> component_rich_composition(std::size_t ncomp, std::size_t rich_index, double min_composition);
double composition_charge(const std::vector<double>& composition, const std::vector<double>& charges);
double l2_norm(const std::vector<double>& values);
std::vector<double> damping_schedule(double damping);

}  // namespace epcsaft::native::equilibrium
