#pragma once

#include <vector>

namespace epcsaft::native::equilibrium_nlp {

struct IdealReactionQuotientResult {
    std::vector<double> log_q;
    std::vector<double> residuals;
    std::vector<double> jacobian_row_major;
};

IdealReactionQuotientResult evaluate_ideal_reaction_quotients(
    const std::vector<double>& amounts,
    int reaction_count,
    const std::vector<double>& stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants
);

std::vector<double> amounts_from_reaction_extents(
    const std::vector<double>& initial_amounts,
    int reaction_count,
    const std::vector<double>& stoichiometry_row_major,
    const std::vector<double>& extents
);

}  // namespace epcsaft::native::equilibrium_nlp
