#pragma once

#include <vector>

namespace epcsaft::native::equilibrium_nlp {

struct IdealReducedGibbsResult {
    double value = 0.0;
    std::vector<double> mole_fractions;
    std::vector<double> gradient;
    std::vector<double> hessian_row_major;
};

IdealReducedGibbsResult evaluate_ideal_reduced_gibbs(
    const std::vector<double>& amounts,
    const std::vector<double>& standard_mu_rt,
    bool include_hessian
);

}  // namespace epcsaft::native::equilibrium_nlp
