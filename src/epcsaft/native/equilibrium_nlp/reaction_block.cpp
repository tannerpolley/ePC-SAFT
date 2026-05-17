#include "reaction_block.h"

#include "epcsaft_electrolyte.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace epcsaft::native::equilibrium_nlp {

namespace {

using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

void validate_reaction_shape(
    int species_count,
    int reaction_count,
    const std::vector<double>& stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants
) {
    if (species_count <= 0) {
        throw ValueError("Ideal reaction block requires at least one species.");
    }
    if (reaction_count <= 0) {
        throw ValueError("Ideal reaction block requires at least one reaction.");
    }
    const std::size_t expected_stoich_size =
        static_cast<std::size_t>(reaction_count) * static_cast<std::size_t>(species_count);
    if (stoichiometry_row_major.size() != expected_stoich_size) {
        throw ValueError("Ideal reaction stoichiometry must be a reaction-by-species row-major matrix.");
    }
    if (log_equilibrium_constants.size() != static_cast<std::size_t>(reaction_count)) {
        throw ValueError("Ideal reaction block requires one equilibrium constant per reaction.");
    }
    for (double coefficient : stoichiometry_row_major) {
        if (!std::isfinite(coefficient)) {
            throw ValueError("Ideal reaction stoichiometry must be finite.");
        }
    }
    for (double log_k : log_equilibrium_constants) {
        if (!std::isfinite(log_k)) {
            throw ValueError("Ideal reaction equilibrium constants must be finite in log form.");
        }
    }
}

void validate_reaction_inputs(
    const std::vector<double>& amounts,
    int reaction_count,
    const std::vector<double>& stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants
) {
    if (amounts.empty()) {
        throw ValueError("Ideal reaction block requires at least one species.");
    }
    const auto species = static_cast<int>(amounts.size());
    validate_reaction_shape(species, reaction_count, stoichiometry_row_major, log_equilibrium_constants);
    for (double amount : amounts) {
        if (!(std::isfinite(amount) && amount > 0.0)) {
            throw ValueError("Ideal reaction species amounts must be positive and finite.");
        }
    }
}

}  // namespace

std::vector<double> amounts_from_reaction_extents(
    const std::vector<double>& initial_amounts,
    int reaction_count,
    const std::vector<double>& stoichiometry_row_major,
    const std::vector<double>& extents
) {
    if (initial_amounts.empty()) {
        throw ValueError("Reaction extent mapping requires at least one species.");
    }
    if (reaction_count <= 0) {
        throw ValueError("Reaction extent mapping requires at least one reaction.");
    }
    const std::size_t species = initial_amounts.size();
    if (stoichiometry_row_major.size() != static_cast<std::size_t>(reaction_count) * species) {
        throw ValueError("Reaction extent stoichiometry must be a reaction-by-species row-major matrix.");
    }
    if (extents.size() != static_cast<std::size_t>(reaction_count)) {
        throw ValueError("Reaction extent vector length must match reaction count.");
    }
    std::vector<double> amounts = initial_amounts;
    for (double value : amounts) {
        if (!(std::isfinite(value) && value >= 0.0)) {
            throw ValueError("Initial species amounts must be finite and non-negative.");
        }
    }
    for (int reaction = 0; reaction < reaction_count; ++reaction) {
        const double extent = extents[static_cast<std::size_t>(reaction)];
        if (!std::isfinite(extent)) {
            throw ValueError("Reaction extents must be finite.");
        }
        for (std::size_t species_index = 0; species_index < species; ++species_index) {
            amounts[species_index] +=
                stoichiometry_row_major[static_cast<std::size_t>(reaction) * species + species_index] * extent;
        }
    }
    for (double value : amounts) {
        if (!(std::isfinite(value) && value > 0.0)) {
            throw ValueError("Reaction extent mapping produced non-positive species amounts.");
        }
    }
    return amounts;
}

std::vector<double> standard_mu_rt_from_reactions(
    int species_count,
    int reaction_count,
    const std::vector<double>& stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants
) {
    validate_reaction_shape(species_count, reaction_count, stoichiometry_row_major, log_equilibrium_constants);
    Eigen::Map<const RowMajorMatrix> stoich(
        stoichiometry_row_major.data(),
        reaction_count,
        species_count
    );
    Eigen::VectorXd rhs(reaction_count);
    for (int row = 0; row < reaction_count; ++row) {
        rhs[row] = -log_equilibrium_constants[static_cast<std::size_t>(row)];
    }
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> decomposition(stoich);
    const Eigen::VectorXd mu = decomposition.solve(rhs);
    const Eigen::VectorXd residual = stoich * mu - rhs;
    const double tolerance = 1.0e-10 * std::max(1.0, rhs.lpNorm<Eigen::Infinity>());
    if (residual.lpNorm<Eigen::Infinity>() > tolerance) {
        throw ValueError("Ideal Gibbs reaction constants are inconsistent with the stoichiometry matrix.");
    }
    return std::vector<double>(mu.data(), mu.data() + mu.size());
}

IdealReactionQuotientResult evaluate_ideal_reaction_quotients(
    const std::vector<double>& amounts,
    int reaction_count,
    const std::vector<double>& stoichiometry_row_major,
    const std::vector<double>& log_equilibrium_constants
) {
    validate_reaction_inputs(amounts, reaction_count, stoichiometry_row_major, log_equilibrium_constants);
    const std::size_t species = amounts.size();
    const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    if (!(std::isfinite(total) && total > 0.0)) {
        throw ValueError("Ideal reaction total amount must be positive and finite.");
    }

    std::vector<double> log_x(species);
    for (std::size_t index = 0; index < species; ++index) {
        log_x[index] = std::log(amounts[index] / total);
    }

    IdealReactionQuotientResult out;
    out.log_q.assign(static_cast<std::size_t>(reaction_count), 0.0);
    out.residuals.assign(static_cast<std::size_t>(reaction_count), 0.0);
    out.jacobian_row_major.assign(static_cast<std::size_t>(reaction_count) * species, 0.0);
    for (int reaction = 0; reaction < reaction_count; ++reaction) {
        double stoich_sum = 0.0;
        for (std::size_t species_index = 0; species_index < species; ++species_index) {
            const double coefficient = stoichiometry_row_major[static_cast<std::size_t>(reaction) * species + species_index];
            out.log_q[static_cast<std::size_t>(reaction)] += coefficient * log_x[species_index];
            stoich_sum += coefficient;
        }
        out.residuals[static_cast<std::size_t>(reaction)] =
            out.log_q[static_cast<std::size_t>(reaction)]
            - log_equilibrium_constants[static_cast<std::size_t>(reaction)];
        for (std::size_t species_index = 0; species_index < species; ++species_index) {
            const double coefficient = stoichiometry_row_major[static_cast<std::size_t>(reaction) * species + species_index];
            out.jacobian_row_major[static_cast<std::size_t>(reaction) * species + species_index] =
                coefficient / amounts[species_index] - stoich_sum / total;
        }
    }
    return out;
}

}  // namespace epcsaft::native::equilibrium_nlp
