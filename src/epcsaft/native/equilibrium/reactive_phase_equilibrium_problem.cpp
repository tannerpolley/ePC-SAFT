#include "epcsaft_equilibrium.h"

#include "equilibrium_helpers.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace {

using epcsaft::native::equilibrium::clip_normalize;
using epcsaft::native::equilibrium::composition_charge;
using epcsaft::native::equilibrium::l2_norm;
using epcsaft::native::equilibrium::max_abs;
using epcsaft::native::equilibrium::normalize_feed;
using epcsaft::native::equilibrium::phase_distance;

std::vector<double> matrix_vector_residual(
    const std::vector<double>& matrix_row_major,
    int rows,
    int cols,
    const std::vector<double>& values,
    const std::vector<double>& target
) {
    if (rows <= 0) {
        throw ValueError("reactive phase residual evaluation requires at least one balance row.");
    }
    if (matrix_row_major.size() != static_cast<std::size_t>(rows * cols)) {
        throw ValueError("balance_matrix has an invalid row-major size.");
    }
    if (target.size() != static_cast<std::size_t>(rows)) {
        throw ValueError("total_vector length must match balance row count.");
    }
    std::vector<double> residual(static_cast<std::size_t>(rows), 0.0);
    for (int r = 0; r < rows; ++r) {
        double value = 0.0;
        for (int c = 0; c < cols; ++c) {
            value += matrix_row_major[static_cast<std::size_t>(r * cols + c)] * values[static_cast<std::size_t>(c)];
        }
        residual[static_cast<std::size_t>(r)] = value - target[static_cast<std::size_t>(r)];
    }
    return residual;
}

std::vector<double> exp_amounts(const std::vector<double>& variables, std::size_t offset, std::size_t count) {
    std::vector<double> out(count, 0.0);
    for (std::size_t i = 0; i < count; ++i) {
        const double value = variables[offset + i];
        if (!std::isfinite(value)) {
            throw ValueError("reactive phase residual variables must be finite.");
        }
        out[i] = std::exp(std::max(-700.0, std::min(700.0, value)));
    }
    return out;
}

double sum_amounts(const std::vector<double>& values) {
    const double total = std::accumulate(values.begin(), values.end(), 0.0);
    if (!std::isfinite(total) || total <= 0.0) {
        throw ValueError("reactive phase residual variables produced invalid phase amounts.");
    }
    return total;
}

std::vector<double> composition_from_amounts(const std::vector<double>& amounts, double total, double floor) {
    std::vector<double> composition(amounts.size(), floor);
    double clipped = 0.0;
    for (std::size_t i = 0; i < amounts.size(); ++i) {
        composition[i] = std::max(amounts[i] / total, floor);
        clipped += composition[i];
    }
    if (!std::isfinite(clipped) || clipped <= 0.0) {
        throw ValueError("reactive phase residual composition normalization failed.");
    }
    for (double& item : composition) {
        item /= clipped;
    }
    return composition;
}

std::vector<double> default_variables_from_feed(const std::vector<double>& feed, double floor) {
    std::vector<double> out(2 * feed.size(), 0.0);
    for (std::size_t i = 0; i < feed.size(); ++i) {
        const double amount = std::max(0.5 * feed[i], floor);
        out[i] = std::log(amount);
        out[feed.size() + i] = std::log(amount);
    }
    return out;
}

std::vector<double> variables_from_initial_phases(
    const std::vector<double>& phase1,
    const std::vector<double>& phase2,
    double beta2,
    std::size_t ncomp,
    double floor
) {
    if (phase1.size() != ncomp || phase2.size() != ncomp) {
        throw ValueError("initial reactive phases must match mixture component count.");
    }
    if (!(std::isfinite(beta2) && beta2 > 0.0 && beta2 < 1.0)) {
        throw ValueError("initial reactive phase fraction must be > 0 and < 1.");
    }
    std::vector<double> x1 = clip_normalize(phase1, floor);
    std::vector<double> x2 = clip_normalize(phase2, floor);
    std::vector<double> out(2 * ncomp, 0.0);
    for (std::size_t i = 0; i < ncomp; ++i) {
        out[i] = std::log(std::max((1.0 - beta2) * x1[i], floor));
        out[ncomp + i] = std::log(std::max(beta2 * x2[i], floor));
    }
    return out;
}

struct ReactivePhaseState {
    std::vector<double> amounts;
    std::vector<double> composition;
    std::vector<double> ln_phi;
    double amount_total = 0.0;
    double density = 0.0;
};

ReactivePhaseState evaluate_phase_state(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& amounts,
    double floor,
    const std::string& density_scope
) {
    ReactivePhaseState out;
    out.amounts = amounts;
    out.amount_total = sum_amounts(amounts);
    out.composition = composition_from_amounts(amounts, out.amount_total, floor);
    out.density = mixture->solve_density_scoped(t, p, out.composition, 0, density_scope);
    std::shared_ptr<ePCSAFTStateNative> state = mixture->state(t, out.composition, 0, false, 0.0, true, out.density);
    out.ln_phi = state->ln_fugacity_coefficient();
    if (out.ln_phi.size() != out.composition.size()) {
        throw ValueError("reactive phase residual fugacity payload length mismatch.");
    }
    return out;
}

std::vector<double> ln_activities(const ReactivePhaseState& state, double floor) {
    std::vector<double> out(state.composition.size(), 0.0);
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = std::log(std::max(state.composition[i], floor)) + state.ln_phi[i];
    }
    return out;
}

std::vector<double> reaction_residuals(
    const std::vector<double>& stoichiometry_row_major,
    int reaction_rows,
    int ncomp,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<double>& ln_activity
) {
    if (reaction_rows < 0) {
        throw ValueError("reaction row count must be non-negative.");
    }
    if (stoichiometry_row_major.size() != static_cast<std::size_t>(reaction_rows * ncomp)) {
        throw ValueError("reaction_stoichiometry has an invalid row-major size.");
    }
    if (log_equilibrium_constants.size() != static_cast<std::size_t>(reaction_rows)) {
        throw ValueError("log equilibrium constant length must match reaction row count.");
    }
    std::vector<double> residual(static_cast<std::size_t>(reaction_rows), 0.0);
    for (int r = 0; r < reaction_rows; ++r) {
        double value = 0.0;
        for (int i = 0; i < ncomp; ++i) {
            value += stoichiometry_row_major[static_cast<std::size_t>(r * ncomp + i)]
                * ln_activity[static_cast<std::size_t>(i)];
        }
        residual[static_cast<std::size_t>(r)] = value - log_equilibrium_constants[static_cast<std::size_t>(r)];
    }
    return residual;
}

std::vector<double> neutral_phase_residuals(
    const std::vector<double>& charges,
    const std::vector<double>& ln_activity1,
    const std::vector<double>& ln_activity2
) {
    std::vector<double> residual;
    for (std::size_t i = 0; i < charges.size(); ++i) {
        if (std::abs(charges[i]) <= 1.0e-12) {
            residual.push_back(ln_activity1[i] - ln_activity2[i]);
        }
    }
    return residual;
}

std::vector<double> ionic_phase_residuals(
    const std::vector<double>& charges,
    const std::vector<double>& ln_activity1,
    const std::vector<double>& ln_activity2
) {
    std::vector<int> cations;
    std::vector<int> anions;
    for (std::size_t i = 0; i < charges.size(); ++i) {
        if (charges[i] > 1.0e-12) {
            cations.push_back(static_cast<int>(i));
        } else if (charges[i] < -1.0e-12) {
            anions.push_back(static_cast<int>(i));
        }
    }
    std::vector<double> residual;
    for (int cation : cations) {
        for (int anion : anions) {
            const double cation_weight = std::abs(charges[static_cast<std::size_t>(anion)]);
            const double anion_weight = std::abs(charges[static_cast<std::size_t>(cation)]);
            residual.push_back(
                cation_weight * (ln_activity1[static_cast<std::size_t>(cation)] - ln_activity2[static_cast<std::size_t>(cation)])
                + anion_weight * (ln_activity1[static_cast<std::size_t>(anion)] - ln_activity2[static_cast<std::size_t>(anion)])
            );
        }
    }
    return residual;
}

void append_block(std::vector<double>& target, const std::vector<double>& block) {
    target.insert(target.end(), block.begin(), block.end());
}

}  // namespace

ReactivePhaseResidualEvaluationNative evaluate_reactive_phase_equilibrium_residual_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& raw_feed,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& variables,
    bool has_variables,
    const std::vector<double>& initial_phase1,
    const std::vector<double>& initial_phase2,
    double initial_phase_fraction_phase2,
    bool has_initial_phases
) {
    const std::size_t ncomp = mixture->ncomp();
    if (reaction_standard_states.size() != static_cast<std::size_t>(reaction_rows)) {
        throw ValueError("reaction standard state length must match reaction row count.");
    }
    std::vector<double> feed = normalize_feed(raw_feed, ncomp, options.min_composition, "reactive_phase_equilibrium");
    std::vector<double> eval_variables = variables;
    if (has_variables) {
        if (variables.size() != 2 * ncomp) {
            throw ValueError("reactive phase residual variables length must be 2 * mixture component count.");
        }
    } else if (has_initial_phases) {
        eval_variables = variables_from_initial_phases(
            initial_phase1,
            initial_phase2,
            initial_phase_fraction_phase2,
            ncomp,
            options.min_composition
        );
    } else {
        eval_variables = default_variables_from_feed(feed, options.min_composition);
    }

    ReactivePhaseState phase1 = evaluate_phase_state(
        mixture,
        t,
        p,
        exp_amounts(eval_variables, 0, ncomp),
        options.min_composition,
        "reactive_phase_equilibrium_phase1"
    );
    ReactivePhaseState phase2 = evaluate_phase_state(
        mixture,
        t,
        p,
        exp_amounts(eval_variables, ncomp, ncomp),
        options.min_composition,
        "reactive_phase_equilibrium_phase2"
    );
    std::vector<double> total_amounts(ncomp, 0.0);
    for (std::size_t i = 0; i < ncomp; ++i) {
        total_amounts[i] = phase1.amounts[i] + phase2.amounts[i];
    }
    std::vector<double> ln_activity1 = ln_activities(phase1, options.min_composition);
    std::vector<double> ln_activity2 = ln_activities(phase2, options.min_composition);
    const int ncomp_int = static_cast<int>(ncomp);

    ReactivePhaseResidualEvaluationNative out;
    out.variable_model = "log_phase_species_amounts";
    out.variables = eval_variables;
    out.lower_bounds.assign(eval_variables.size(), std::log(options.min_composition));
    out.upper_bounds.assign(eval_variables.size(), 50.0);
    out.phase1_amounts = phase1.amounts;
    out.phase2_amounts = phase2.amounts;
    out.phase1_composition = phase1.composition;
    out.phase2_composition = phase2.composition;
    out.phase1_ln_fugacity_coefficient = phase1.ln_phi;
    out.phase2_ln_fugacity_coefficient = phase2.ln_phi;
    out.phase1_density = phase1.density;
    out.phase2_density = phase2.density;
    const double total_phase_amount = phase1.amount_total + phase2.amount_total;
    out.phase_fraction_phase2 = phase2.amount_total / total_phase_amount;
    out.element_balance_residuals = matrix_vector_residual(
        balance_matrix_row_major,
        balance_rows,
        ncomp_int,
        total_amounts,
        total_vector
    );
    out.reaction_residuals_phase1 = reaction_residuals(
        reaction_stoichiometry_row_major,
        reaction_rows,
        ncomp_int,
        log_equilibrium_constants,
        ln_activity1
    );
    out.reaction_residuals_phase2 = reaction_residuals(
        reaction_stoichiometry_row_major,
        reaction_rows,
        ncomp_int,
        log_equilibrium_constants,
        ln_activity2
    );
    const std::vector<double>& charges = mixture->args().z;
    out.phase_charge_residuals = {
        composition_charge(phase1.amounts, charges),
        composition_charge(phase2.amounts, charges),
    };
    out.neutral_phase_equilibrium_residuals = neutral_phase_residuals(charges, ln_activity1, ln_activity2);
    out.ionic_equilibrium_residuals = ionic_phase_residuals(charges, ln_activity1, ln_activity2);
    append_block(out.residual, out.element_balance_residuals);
    append_block(out.residual, out.reaction_residuals_phase1);
    append_block(out.residual, out.reaction_residuals_phase2);
    append_block(out.residual, out.neutral_phase_equilibrium_residuals);
    append_block(out.residual, out.ionic_equilibrium_residuals);
    append_block(out.residual, out.phase_charge_residuals);
    for (double value : out.residual) {
        out.objective += 0.5 * value * value;
    }
    out.gradient.assign(eval_variables.size(), 0.0);
    out.jacobian_rows = static_cast<int>(out.residual.size());
    out.jacobian_cols = static_cast<int>(eval_variables.size());
    out.phase_distance = phase_distance(phase1.composition, phase2.composition);

    out.diagnostics_string["residual_surface"] = "native_reactive_phase_equilibrium_coupled_state";
    out.diagnostics_string["variable_model"] = out.variable_model;
    out.diagnostics_string["residual_blocks"] = "element_balance,reaction_equilibrium,neutral_phase_equilibrium,ionic_equilibrium,phase_charge";
    out.diagnostics_string["solver_backend"] = "residual_surface_only";
    out.diagnostics_string["solver_method"] = "not_started_until_ceres_slice";
    out.diagnostics_string["jacobian_backend"] = "not_available";
    out.diagnostics_string["derivative_backend"] = "not_available";
    out.diagnostics_string["coupling_level"] = "single_native_residual_state";
    out.diagnostics_string["phase_model"] = "two_liquid_phases";
    out.diagnostics_bool["jacobian_available"] = false;
    out.diagnostics_bool["derivative_available"] = false;
    out.diagnostics_bool["reaction_and_phase_residuals_share_state"] = true;
    out.diagnostics_bool["nonnegative_amounts_enforced_by_transform"] = true;
    out.diagnostics_bool["composition_normalization_enforced_by_transform"] = true;
    out.diagnostics_bool["ceres_accepted_solve"] = false;
    out.diagnostics_int["phase_count"] = 2;
    out.diagnostics_int["component_count"] = static_cast<int>(ncomp);
    out.diagnostics_int["reaction_count"] = reaction_rows;
    out.diagnostics_int["balance_row_count"] = balance_rows;
    out.diagnostics_int["variable_count"] = static_cast<int>(eval_variables.size());
    out.diagnostics_int["residual_size"] = static_cast<int>(out.residual.size());
    out.diagnostics_int["element_balance_residual_size"] = static_cast<int>(out.element_balance_residuals.size());
    out.diagnostics_int["reaction_residual_size_per_phase"] = static_cast<int>(out.reaction_residuals_phase1.size());
    out.diagnostics_int["neutral_phase_equilibrium_residual_size"] = static_cast<int>(out.neutral_phase_equilibrium_residuals.size());
    out.diagnostics_int["ionic_equilibrium_residual_size"] = static_cast<int>(out.ionic_equilibrium_residuals.size());
    out.diagnostics_int["phase_charge_residual_size"] = static_cast<int>(out.phase_charge_residuals.size());
    out.diagnostics_double["element_balance_norm"] = max_abs(out.element_balance_residuals);
    out.diagnostics_double["reaction_residual_norm"] = std::max(
        max_abs(out.reaction_residuals_phase1),
        max_abs(out.reaction_residuals_phase2)
    );
    out.diagnostics_double["phase_equilibrium_residual_norm"] = std::max(
        max_abs(out.neutral_phase_equilibrium_residuals),
        max_abs(out.ionic_equilibrium_residuals)
    );
    out.diagnostics_double["phase_charge_balance_norm"] = max_abs(out.phase_charge_residuals);
    out.diagnostics_double["residual_norm"] = max_abs(out.residual);
    out.diagnostics_double["residual_l2_norm"] = l2_norm(out.residual);
    out.diagnostics_double["objective"] = out.objective;
    out.diagnostics_double["phase_distance"] = out.phase_distance;
    out.diagnostics_double["phase_fraction_phase2"] = out.phase_fraction_phase2;
    out.diagnostics_vector["feed_composition"] = feed;
    out.diagnostics_vector["total_phase_amounts"] = total_amounts;
    out.diagnostics_vector["element_balance_residual"] = out.element_balance_residuals;
    out.diagnostics_vector["reaction_residual_phase1"] = out.reaction_residuals_phase1;
    out.diagnostics_vector["reaction_residual_phase2"] = out.reaction_residuals_phase2;
    out.diagnostics_vector["neutral_phase_equilibrium_residual"] = out.neutral_phase_equilibrium_residuals;
    out.diagnostics_vector["ionic_equilibrium_residual"] = out.ionic_equilibrium_residuals;
    out.diagnostics_vector["phase_charge_residual"] = out.phase_charge_residuals;
    out.diagnostics_vector["reaction_standard_states"] = std::vector<double>(
        reaction_standard_states.begin(),
        reaction_standard_states.end()
    );
    return out;
}
