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
double composition_charge(const std::vector<double>& composition, const std::vector<double>& charges);
double l2_norm(const std::vector<double>& values);

}  // namespace epcsaft::native::equilibrium
