#include "epcsaft_equilibrium.h"

#include "equilibrium_helpers.h"
#include "../equilibrium_nlp/ipopt_adapter.h"
#include "../equilibrium_nlp/nlp_problem.h"
#include "../equilibrium_nlp/second_order.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>

PhaseStateCompositionSensitivityResult phase_state_ln_fugacity_composition_sensitivity_cpp(
    double t,
    double p,
    std::vector<double> x,
    int phase,
    const add_args& cppargs
);

namespace {

using epcsaft::native::equilibrium::clip_normalize;
using epcsaft::native::equilibrium::composition_charge;
using epcsaft::native::equilibrium::l2_norm;
using epcsaft::native::equilibrium::max_abs;
using epcsaft::native::equilibrium::normalize_feed;
using epcsaft::native::equilibrium::phase_distance;

constexpr int STANDARD_STATE_MOLE_FRACTION_ACTIVITY = 0;
constexpr int STANDARD_STATE_IDEAL_MOLE_FRACTION = 1;
constexpr int STANDARD_STATE_CONCENTRATION = 2;

std::string reaction_standard_state_name(int standard_state) {
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return "mole_fraction_activity";
    }
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return "ideal_mole_fraction";
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        return "concentration";
    }
    return "unsupported";
}

std::string reaction_standard_state_summary(const std::vector<int>& reaction_standard_states) {
    if (reaction_standard_states.empty()) {
        return "none";
    }
    std::ostringstream out;
    for (std::size_t i = 0; i < reaction_standard_states.size(); ++i) {
        if (i > 0) {
            out << ",";
        }
        out << reaction_standard_state_name(reaction_standard_states[i]);
    }
    return out.str();
}

void validate_reaction_standard_states(const std::vector<int>& reaction_standard_states, int reaction_rows) {
    if (reaction_standard_states.size() != static_cast<std::size_t>(reaction_rows)) {
        throw ValueError("reaction_standard_states length must match reaction row count.");
    }
    for (int standard_state : reaction_standard_states) {
        if (standard_state != STANDARD_STATE_MOLE_FRACTION_ACTIVITY
            && standard_state != STANDARD_STATE_IDEAL_MOLE_FRACTION
            && standard_state != STANDARD_STATE_CONCENTRATION) {
            throw ValueError("reaction standard state code is outside the native reactive-phase contract.");
        }
    }
}

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

std::vector<double> species_amount_upper_bounds(
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    std::size_t ncomp,
    double floor
) {
    if (balance_matrix_row_major.size() != static_cast<std::size_t>(balance_rows) * ncomp) {
        throw ValueError("reactive phase upper-bound builder received an invalid balance matrix.");
    }
    if (total_vector.size() != static_cast<std::size_t>(balance_rows)) {
        throw ValueError("reactive phase upper-bound builder received an invalid total vector.");
    }

    double default_upper = 0.0;
    for (double total : total_vector) {
        if (std::isfinite(total) && total > 0.0) {
            default_upper += total;
        }
    }
    default_upper = std::max(default_upper, 1.0);

    std::vector<double> upper(ncomp, default_upper);
    for (std::size_t species = 0; species < ncomp; ++species) {
        double limit = std::numeric_limits<double>::infinity();
        for (int row = 0; row < balance_rows; ++row) {
            const double coefficient = balance_matrix_row_major[static_cast<std::size_t>(row) * ncomp + species];
            const double total = total_vector[static_cast<std::size_t>(row)];
            if (std::isfinite(coefficient) && coefficient > 0.0 && std::isfinite(total) && total >= 0.0) {
                limit = std::min(limit, total / coefficient);
            }
        }
        if (std::isfinite(limit)) {
            upper[species] = std::min(upper[species], limit);
        }
        upper[species] = std::max(upper[species], 10.0 * floor);
    }
    return upper;
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
    std::vector<double> model_composition;
    std::vector<double> ln_phi;
    std::vector<int> global_indices;
    std::size_t global_component_count = 0;
    double amount_total = 0.0;
    double density = 0.0;
};

std::vector<int> full_component_indices(std::size_t ncomp) {
    std::vector<int> indices(ncomp, 0);
    for (std::size_t i = 0; i < ncomp; ++i) {
        indices[i] = static_cast<int>(i);
    }
    return indices;
}

bool contains_index(const std::vector<int>& indices, std::size_t index) {
    return std::find(indices.begin(), indices.end(), static_cast<int>(index)) != indices.end();
}

std::vector<double> selected_amounts(
    const std::vector<double>& amounts,
    const std::vector<int>& indices,
    std::size_t global_ncomp
) {
    std::vector<double> out;
    out.reserve(indices.size());
    for (int index : indices) {
        if (index < 0 || static_cast<std::size_t>(index) >= global_ncomp) {
            throw ValueError("phase model component index is outside the global reactive mixture.");
        }
        out.push_back(amounts[static_cast<std::size_t>(index)]);
    }
    return out;
}

void validate_phase_model_indices(
    const std::shared_ptr<ePCSAFTMixtureNative>& phase_mixture,
    const std::vector<int>& indices,
    std::size_t global_ncomp,
    const std::string& label
) {
    if (!phase_mixture) {
        throw ValueError("reactive electrolyte phase_models require native aqueous and organic mixtures.");
    }
    if (indices.size() != phase_mixture->ncomp()) {
        throw ValueError(label + " phase model species count must match its global index map.");
    }
    std::vector<int> seen;
    seen.reserve(indices.size());
    for (int index : indices) {
        if (index < 0 || static_cast<std::size_t>(index) >= global_ncomp) {
            throw ValueError(label + " phase model index is outside the global reactive mixture.");
        }
        if (std::find(seen.begin(), seen.end(), index) != seen.end()) {
            throw ValueError(label + " phase model index map contains duplicates.");
        }
        seen.push_back(index);
    }
}

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
    out.model_composition = out.composition;
    out.global_component_count = out.composition.size();
    out.global_indices = full_component_indices(out.composition.size());
    out.density = mixture->solve_density_scoped(t, p, out.composition, 0, density_scope);
    std::shared_ptr<ePCSAFTStateNative> state = mixture->state(t, out.composition, 0, false, 0.0, true, out.density);
    out.ln_phi = state->ln_fugacity_coefficient();
    if (out.ln_phi.size() != out.composition.size()) {
        throw ValueError("reactive phase residual fugacity payload length mismatch.");
    }
    return out;
}

ReactivePhaseState evaluate_scoped_phase_state(
    const std::shared_ptr<ePCSAFTMixtureNative>& phase_mixture,
    double t,
    double p,
    const std::vector<double>& amounts,
    double floor,
    const std::string& density_scope,
    const std::vector<int>& global_indices,
    std::size_t global_ncomp
);

std::vector<double> ln_activities(const ReactivePhaseState& state, double floor) {
    std::vector<double> out(state.composition.size(), 0.0);
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = std::log(std::max(state.composition[i], floor)) + state.ln_phi[i];
    }
    return out;
}

std::vector<double> reaction_standard_state_log_terms(
    const ReactivePhaseState& state,
    int standard_state,
    double floor
) {
    std::vector<double> out(state.composition.size(), 0.0);
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return ln_activities(state, floor);
    }
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = std::log(std::max(state.composition[i], floor));
        }
        return out;
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        if (!(std::isfinite(state.density) && state.density > 0.0)) {
            throw ValueError("concentration reaction standard state requires a finite positive molar density.");
        }
        for (std::size_t i = 0; i < out.size(); ++i) {
            out[i] = std::log(std::max(state.composition[i] * state.density, floor));
        }
        return out;
    }
    throw ValueError("reaction standard state code is outside the native reactive-phase contract.");
}

std::vector<double> reaction_residuals(
    const std::vector<double>& stoichiometry_row_major,
    int reaction_rows,
    int ncomp,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const ReactivePhaseState& phase,
    double floor
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
    validate_reaction_standard_states(reaction_standard_states, reaction_rows);
    std::vector<double> residual(static_cast<std::size_t>(reaction_rows), 0.0);
    for (int r = 0; r < reaction_rows; ++r) {
        const std::vector<double> log_terms = reaction_standard_state_log_terms(
            phase,
            reaction_standard_states[static_cast<std::size_t>(r)],
            floor
        );
        double value = 0.0;
        for (int i = 0; i < ncomp; ++i) {
            value += stoichiometry_row_major[static_cast<std::size_t>(r * ncomp + i)]
                * log_terms[static_cast<std::size_t>(i)];
        }
        residual[static_cast<std::size_t>(r)] = value - log_equilibrium_constants[static_cast<std::size_t>(r)];
    }
    return residual;
}

bool has_phase_tagged_reaction_stoichiometry(
    const std::vector<double>& phase_stoichiometry_row_major,
    int reaction_rows,
    int ncomp
) {
    if (phase_stoichiometry_row_major.empty()) {
        return false;
    }
    if (reaction_rows < 0) {
        throw ValueError("reaction row count must be non-negative.");
    }
    if (phase_stoichiometry_row_major.size() != static_cast<std::size_t>(reaction_rows * 2 * ncomp)) {
        throw ValueError("reaction_phase_stoichiometry has an invalid row-major size.");
    }
    return true;
}

std::vector<double> cross_phase_reaction_residuals(
    const std::vector<double>& phase_stoichiometry_row_major,
    int reaction_rows,
    int ncomp,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const ReactivePhaseState& phase1,
    const ReactivePhaseState& phase2,
    double floor
) {
    if (!has_phase_tagged_reaction_stoichiometry(phase_stoichiometry_row_major, reaction_rows, ncomp)) {
        return {};
    }
    if (log_equilibrium_constants.size() != static_cast<std::size_t>(reaction_rows)) {
        throw ValueError("log equilibrium constant length must match reaction row count.");
    }
    validate_reaction_standard_states(reaction_standard_states, reaction_rows);
    std::vector<double> residual(static_cast<std::size_t>(reaction_rows), 0.0);
    for (int r = 0; r < reaction_rows; ++r) {
        const std::vector<double> log_terms1 = reaction_standard_state_log_terms(
            phase1,
            reaction_standard_states[static_cast<std::size_t>(r)],
            floor
        );
        const std::vector<double> log_terms2 = reaction_standard_state_log_terms(
            phase2,
            reaction_standard_states[static_cast<std::size_t>(r)],
            floor
        );
        double value = 0.0;
        const std::size_t row_offset = static_cast<std::size_t>(r * 2 * ncomp);
        for (int i = 0; i < ncomp; ++i) {
            value += phase_stoichiometry_row_major[row_offset + static_cast<std::size_t>(i)]
                * log_terms1[static_cast<std::size_t>(i)];
            value += phase_stoichiometry_row_major[row_offset + static_cast<std::size_t>(ncomp + i)]
                * log_terms2[static_cast<std::size_t>(i)];
        }
        residual[static_cast<std::size_t>(r)] = value - log_equilibrium_constants[static_cast<std::size_t>(r)];
    }
    return residual;
}

double reaction_residual_norm(const ReactivePhaseResidualEvaluationNative& eval) {
    if (!eval.reaction_residuals_cross_phase.empty()) {
        return max_abs(eval.reaction_residuals_cross_phase);
    }
    return std::max(max_abs(eval.reaction_residuals_phase1), max_abs(eval.reaction_residuals_phase2));
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

std::vector<double> phase_log_mole_fraction_log_amount_jacobian(const ReactivePhaseState& state) {
    const std::size_t ncomp = state.global_component_count > 0
        ? state.global_component_count
        : state.composition.size();
    const std::vector<int> indices = state.global_indices.empty()
        ? full_component_indices(ncomp)
        : state.global_indices;
    const std::vector<double>& composition = state.model_composition.empty()
        ? state.composition
        : state.model_composition;
    std::vector<double> jacobian(ncomp * ncomp, 0.0);
    for (std::size_t species = 0; species < indices.size(); ++species) {
        const std::size_t global_species = static_cast<std::size_t>(indices[species]);
        for (std::size_t variable = 0; variable < indices.size(); ++variable) {
            const std::size_t global_variable = static_cast<std::size_t>(indices[variable]);
            jacobian[global_species * ncomp + global_variable] =
                (species == variable ? 1.0 : 0.0) - composition[variable];
        }
    }
    return jacobian;
}

std::vector<double> phase_composition_log_amount_hessian_tensor(const ReactivePhaseState& state) {
    const std::size_t ncomp = state.global_component_count > 0
        ? state.global_component_count
        : state.composition.size();
    const std::vector<int> indices = state.global_indices.empty()
        ? full_component_indices(ncomp)
        : state.global_indices;
    const std::vector<double>& composition = state.model_composition.empty()
        ? state.composition
        : state.model_composition;
    if (composition.size() != indices.size()) {
        throw ValueError("reactive phase composition Hessian shape does not match phase indices.");
    }
    std::vector<double> hessian(ncomp * ncomp * ncomp, 0.0);
    for (std::size_t local_species = 0; local_species < indices.size(); ++local_species) {
        const std::size_t global_species = static_cast<std::size_t>(indices[local_species]);
        const double x_species = composition[local_species];
        for (std::size_t first = 0; first < indices.size(); ++first) {
            const std::size_t global_first = static_cast<std::size_t>(indices[first]);
            const double first_delta = local_species == first ? 1.0 : 0.0;
            const double x_first = composition[first];
            for (std::size_t second = 0; second < indices.size(); ++second) {
                const std::size_t global_second = static_cast<std::size_t>(indices[second]);
                const double second_delta = local_species == second ? 1.0 : 0.0;
                const double cross_delta = first == second ? 1.0 : 0.0;
                const double x_second = composition[second];
                hessian[global_species * ncomp * ncomp + global_first * ncomp + global_second] =
                    x_species * (second_delta - x_second) * (first_delta - x_first)
                    - x_species * x_first * (cross_delta - x_second);
            }
        }
    }
    return hessian;
}

std::vector<double> phase_log_mole_fraction_log_amount_hessian_tensor(const ReactivePhaseState& state) {
    const std::size_t ncomp = state.global_component_count > 0
        ? state.global_component_count
        : state.composition.size();
    const std::vector<int> indices = state.global_indices.empty()
        ? full_component_indices(ncomp)
        : state.global_indices;
    const std::vector<double>& composition = state.model_composition.empty()
        ? state.composition
        : state.model_composition;
    if (composition.size() != indices.size()) {
        throw ValueError("reactive phase log-composition Hessian shape does not match phase indices.");
    }
    std::vector<double> hessian(ncomp * ncomp * ncomp, 0.0);
    for (std::size_t local_species = 0; local_species < indices.size(); ++local_species) {
        const std::size_t global_species = static_cast<std::size_t>(indices[local_species]);
        for (std::size_t first = 0; first < indices.size(); ++first) {
            const std::size_t global_first = static_cast<std::size_t>(indices[first]);
            const double x_first = composition[first];
            for (std::size_t second = 0; second < indices.size(); ++second) {
                const std::size_t global_second = static_cast<std::size_t>(indices[second]);
                const double cross_delta = first == second ? 1.0 : 0.0;
                hessian[global_species * ncomp * ncomp + global_first * ncomp + global_second] =
                    -x_first * (cross_delta - composition[second]);
            }
        }
    }
    return hessian;
}

std::vector<double> phase_ln_activity_log_amount_jacobian(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ReactivePhaseState& state
) {
    const std::size_t ncomp = state.global_component_count > 0
        ? state.global_component_count
        : state.composition.size();
    const std::vector<int> indices = state.global_indices.empty()
        ? full_component_indices(ncomp)
        : state.global_indices;
    const std::vector<double>& composition = state.model_composition.empty()
        ? state.composition
        : state.model_composition;
    const std::size_t local_ncomp = composition.size();
    PhaseStateCompositionSensitivityResult sensitivity =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, composition, 0, mixture->args());
    if (!sensitivity.supported || sensitivity.jacobian_row_major.size() != local_ncomp * local_ncomp) {
        throw ValueError("reactive phase residual Jacobian requires supported phase-state fugacity sensitivities.");
    }
    std::vector<double> jacobian(ncomp * ncomp, 0.0);
    for (std::size_t species = 0; species < indices.size(); ++species) {
        const std::size_t global_species = static_cast<std::size_t>(indices[species]);
        for (std::size_t variable = 0; variable < indices.size(); ++variable) {
            const std::size_t global_variable = static_cast<std::size_t>(indices[variable]);
            const double dlogx = (species == variable ? 1.0 : 0.0) - composition[variable];
            double dlnphi = 0.0;
            for (std::size_t k = 0; k < local_ncomp; ++k) {
                const double dxk_dlogn =
                    composition[k] * ((k == variable ? 1.0 : 0.0) - composition[variable]);
                dlnphi += sensitivity.jacobian_row_major[species * local_ncomp + k] * dxk_dlogn;
            }
            jacobian[global_species * ncomp + global_variable] = dlogx + dlnphi;
        }
    }
    return jacobian;
}

std::vector<double> phase_ln_activity_log_amount_hessian_tensor(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ReactivePhaseState& state
) {
    const std::size_t ncomp = state.global_component_count > 0
        ? state.global_component_count
        : state.composition.size();
    const std::vector<int> indices = state.global_indices.empty()
        ? full_component_indices(ncomp)
        : state.global_indices;
    const std::vector<double>& composition = state.model_composition.empty()
        ? state.composition
        : state.model_composition;
    const std::size_t local_ncomp = composition.size();
    PhaseStateCompositionSensitivityResult sensitivity =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, composition, 0, mixture->args());
    if (!sensitivity.supported
        || sensitivity.jacobian_row_major.size() != local_ncomp * local_ncomp
        || sensitivity.hessian_tensor_row_major.size() != local_ncomp * local_ncomp * local_ncomp) {
        throw ValueError("reactive phase residual Hessian requires supported phase-state fugacity Hessians.");
    }
    std::vector<double> hessian = phase_log_mole_fraction_log_amount_hessian_tensor(state);
    auto dx = [&](std::size_t species, std::size_t variable) {
        return composition[species] * ((species == variable ? 1.0 : 0.0) - composition[variable]);
    };
    auto d2x = [&](std::size_t species, std::size_t first, std::size_t second) {
        return composition[species]
            * ((species == second ? 1.0 : 0.0) - composition[second])
            * ((species == first ? 1.0 : 0.0) - composition[first])
            - composition[species] * composition[first] * ((first == second ? 1.0 : 0.0) - composition[second]);
    };
    for (std::size_t local_species = 0; local_species < indices.size(); ++local_species) {
        const std::size_t global_species = static_cast<std::size_t>(indices[local_species]);
        for (std::size_t first = 0; first < indices.size(); ++first) {
            const std::size_t global_first = static_cast<std::size_t>(indices[first]);
            for (std::size_t second = 0; second < indices.size(); ++second) {
                const std::size_t global_second = static_cast<std::size_t>(indices[second]);
                double value = hessian[global_species * ncomp * ncomp + global_first * ncomp + global_second];
                for (std::size_t k = 0; k < local_ncomp; ++k) {
                    value += sensitivity.jacobian_row_major[local_species * local_ncomp + k]
                        * d2x(k, first, second);
                    for (std::size_t l = 0; l < local_ncomp; ++l) {
                        value += sensitivity.hessian_tensor_row_major[
                            local_species * local_ncomp * local_ncomp + k * local_ncomp + l
                        ] * dx(k, first) * dx(l, second);
                    }
                }
                hessian[global_species * ncomp * ncomp + global_first * ncomp + global_second] = value;
            }
        }
    }
    return hessian;
}

std::vector<double> phase_log_concentration_log_amount_jacobian(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ReactivePhaseState& state
) {
    const std::size_t ncomp = state.global_component_count > 0
        ? state.global_component_count
        : state.composition.size();
    const std::vector<int> indices = state.global_indices.empty()
        ? full_component_indices(ncomp)
        : state.global_indices;
    const std::vector<double>& composition = state.model_composition.empty()
        ? state.composition
        : state.model_composition;
    const std::size_t local_ncomp = composition.size();
    PhaseStateCompositionSensitivityResult sensitivity =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, composition, 0, mixture->args());
    if (!sensitivity.supported || sensitivity.density_composition_derivative.size() != local_ncomp) {
        throw ValueError("concentration reaction residual Jacobian requires supported density composition sensitivities.");
    }
    const double density = sensitivity.density > 0.0 ? sensitivity.density : state.density;
    if (!(std::isfinite(density) && density > 0.0)) {
        throw ValueError("concentration reaction residual Jacobian requires a finite positive molar density.");
    }
    std::vector<double> jacobian = phase_log_mole_fraction_log_amount_jacobian(state);
    std::vector<double> dlogrho_dlogn(local_ncomp, 0.0);
    for (std::size_t variable = 0; variable < local_ncomp; ++variable) {
        double drho_dlogn = 0.0;
        for (std::size_t k = 0; k < local_ncomp; ++k) {
            const double dxk_dlogn =
                composition[k] * ((k == variable ? 1.0 : 0.0) - composition[variable]);
            drho_dlogn += sensitivity.density_composition_derivative[k] * dxk_dlogn;
        }
        dlogrho_dlogn[variable] = drho_dlogn / density;
    }
    for (std::size_t species = 0; species < indices.size(); ++species) {
        const std::size_t global_species = static_cast<std::size_t>(indices[species]);
        for (std::size_t variable = 0; variable < indices.size(); ++variable) {
            const std::size_t global_variable = static_cast<std::size_t>(indices[variable]);
            jacobian[global_species * ncomp + global_variable] += dlogrho_dlogn[variable];
        }
    }
    return jacobian;
}

std::vector<double> phase_log_concentration_log_amount_hessian_tensor(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ReactivePhaseState& state
) {
    const std::size_t ncomp = state.global_component_count > 0
        ? state.global_component_count
        : state.composition.size();
    const std::vector<int> indices = state.global_indices.empty()
        ? full_component_indices(ncomp)
        : state.global_indices;
    const std::vector<double>& composition = state.model_composition.empty()
        ? state.composition
        : state.model_composition;
    const std::size_t local_ncomp = composition.size();
    PhaseStateCompositionSensitivityResult sensitivity =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, composition, 0, mixture->args());
    if (!sensitivity.supported
        || sensitivity.density_composition_derivative.size() != local_ncomp
        || sensitivity.density_composition_hessian_row_major.size() != local_ncomp * local_ncomp) {
        throw ValueError("concentration reaction residual Hessian requires supported density composition Hessians.");
    }
    const double density = sensitivity.density > 0.0 ? sensitivity.density : state.density;
    if (!(std::isfinite(density) && density > 0.0)) {
        throw ValueError("concentration reaction residual Hessian requires a finite positive molar density.");
    }
    std::vector<double> hessian = phase_log_mole_fraction_log_amount_hessian_tensor(state);
    auto dx = [&](std::size_t species, std::size_t variable) {
        return composition[species] * ((species == variable ? 1.0 : 0.0) - composition[variable]);
    };
    auto d2x = [&](std::size_t species, std::size_t first, std::size_t second) {
        return composition[species]
            * ((species == second ? 1.0 : 0.0) - composition[second])
            * ((species == first ? 1.0 : 0.0) - composition[first])
            - composition[species] * composition[first] * ((first == second ? 1.0 : 0.0) - composition[second]);
    };
    std::vector<double> dlogrho(local_ncomp, 0.0);
    for (std::size_t variable = 0; variable < local_ncomp; ++variable) {
        for (std::size_t k = 0; k < local_ncomp; ++k) {
            dlogrho[variable] += sensitivity.density_composition_derivative[k] * dx(k, variable);
        }
        dlogrho[variable] /= density;
    }
    std::vector<double> d2logrho(local_ncomp * local_ncomp, 0.0);
    for (std::size_t first = 0; first < local_ncomp; ++first) {
        for (std::size_t second = 0; second < local_ncomp; ++second) {
            double d2rho = 0.0;
            for (std::size_t k = 0; k < local_ncomp; ++k) {
                d2rho += sensitivity.density_composition_derivative[k] * d2x(k, first, second);
                for (std::size_t l = 0; l < local_ncomp; ++l) {
                    d2rho += sensitivity.density_composition_hessian_row_major[k * local_ncomp + l]
                        * dx(k, first) * dx(l, second);
                }
            }
            d2logrho[first * local_ncomp + second] =
                d2rho / density - dlogrho[first] * dlogrho[second];
        }
    }
    for (std::size_t local_species = 0; local_species < indices.size(); ++local_species) {
        const std::size_t global_species = static_cast<std::size_t>(indices[local_species]);
        for (std::size_t first = 0; first < indices.size(); ++first) {
            const std::size_t global_first = static_cast<std::size_t>(indices[first]);
            for (std::size_t second = 0; second < indices.size(); ++second) {
                const std::size_t global_second = static_cast<std::size_t>(indices[second]);
                hessian[global_species * ncomp * ncomp + global_first * ncomp + global_second] +=
                    d2logrho[first * local_ncomp + second];
            }
        }
    }
    return hessian;
}

std::vector<double> phase_reaction_standard_state_log_amount_jacobian(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ReactivePhaseState& state,
    int standard_state
) {
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return phase_ln_activity_log_amount_jacobian(mixture, t, p, state);
    }
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return phase_log_mole_fraction_log_amount_jacobian(state);
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        return phase_log_concentration_log_amount_jacobian(mixture, t, p, state);
    }
    throw ValueError("reaction standard state code is outside the native reactive-phase contract.");
}

std::vector<double> phase_reaction_standard_state_log_amount_hessian_tensor(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ReactivePhaseState& state,
    int standard_state
) {
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return phase_ln_activity_log_amount_hessian_tensor(mixture, t, p, state);
    }
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return phase_log_mole_fraction_log_amount_hessian_tensor(state);
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        return phase_log_concentration_log_amount_hessian_tensor(mixture, t, p, state);
    }
    throw ValueError("reaction standard state code is outside the native reactive-phase contract.");
}

std::vector<double> reactive_phase_residual_jacobian_row_major(
    const std::shared_ptr<ePCSAFTMixtureNative>& phase1_mixture,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase2_mixture,
    double t,
    double p,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    const ReactivePhaseState& phase1,
    const ReactivePhaseState& phase2,
    const std::vector<double>& charges,
    const ReactivePhaseResidualEvaluationNative& residual_eval
) {
    const std::size_t ncomp = phase1.composition.size();
    const std::size_t nvars = 2 * ncomp;
    std::vector<double> jacobian(residual_eval.residual.size() * nvars, 0.0);
    std::vector<double> activity_jac1 = phase_ln_activity_log_amount_jacobian(phase1_mixture, t, p, phase1);
    std::vector<double> activity_jac2 = phase_ln_activity_log_amount_jacobian(phase2_mixture, t, p, phase2);
    validate_reaction_standard_states(reaction_standard_states, reaction_rows);
    const bool phase_tagged_reactions = has_phase_tagged_reaction_stoichiometry(
        reaction_phase_stoichiometry_row_major,
        reaction_rows,
        static_cast<int>(ncomp)
    );
    std::size_t row = 0;

    for (int balance = 0; balance < balance_rows; ++balance) {
        for (std::size_t species = 0; species < ncomp; ++species) {
            const double coefficient = balance_matrix_row_major[static_cast<std::size_t>(balance) * ncomp + species];
            jacobian[row * nvars + species] = coefficient * phase1.amounts[species];
            jacobian[row * nvars + ncomp + species] = coefficient * phase2.amounts[species];
        }
        ++row;
    }

    if (phase_tagged_reactions) {
        for (int reaction = 0; reaction < reaction_rows; ++reaction) {
            const std::vector<double> reaction_jac1 = phase_reaction_standard_state_log_amount_jacobian(
                phase1_mixture,
                t,
                p,
                phase1,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            const std::vector<double> reaction_jac2 = phase_reaction_standard_state_log_amount_jacobian(
                phase2_mixture,
                t,
                p,
                phase2,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            const std::size_t reaction_offset = static_cast<std::size_t>(reaction) * 2 * ncomp;
            for (std::size_t variable = 0; variable < ncomp; ++variable) {
                double phase1_value = 0.0;
                double phase2_value = 0.0;
                for (std::size_t species = 0; species < ncomp; ++species) {
                    phase1_value += reaction_phase_stoichiometry_row_major[reaction_offset + species]
                        * reaction_jac1[species * ncomp + variable];
                    phase2_value += reaction_phase_stoichiometry_row_major[reaction_offset + ncomp + species]
                        * reaction_jac2[species * ncomp + variable];
                }
                jacobian[row * nvars + variable] = phase1_value;
                jacobian[row * nvars + ncomp + variable] = phase2_value;
            }
            ++row;
        }
    } else {
        for (int reaction = 0; reaction < reaction_rows; ++reaction) {
            const std::vector<double> reaction_jac1 = phase_reaction_standard_state_log_amount_jacobian(
                phase1_mixture,
                t,
                p,
                phase1,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            for (std::size_t variable = 0; variable < ncomp; ++variable) {
                double value = 0.0;
                for (std::size_t species = 0; species < ncomp; ++species) {
                    value += reaction_stoichiometry_row_major[static_cast<std::size_t>(reaction) * ncomp + species]
                        * reaction_jac1[species * ncomp + variable];
                }
                jacobian[row * nvars + variable] = value;
            }
            ++row;
        }

        for (int reaction = 0; reaction < reaction_rows; ++reaction) {
            const std::vector<double> reaction_jac2 = phase_reaction_standard_state_log_amount_jacobian(
                phase2_mixture,
                t,
                p,
                phase2,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            for (std::size_t variable = 0; variable < ncomp; ++variable) {
                double value = 0.0;
                for (std::size_t species = 0; species < ncomp; ++species) {
                    value += reaction_stoichiometry_row_major[static_cast<std::size_t>(reaction) * ncomp + species]
                        * reaction_jac2[species * ncomp + variable];
                }
                jacobian[row * nvars + ncomp + variable] = value;
            }
            ++row;
        }
    }

    if (!phase_tagged_reactions) {
        for (std::size_t species = 0; species < ncomp; ++species) {
            if (std::abs(charges[species]) > 1.0e-12) {
                continue;
            }
            for (std::size_t variable = 0; variable < ncomp; ++variable) {
                jacobian[row * nvars + variable] = activity_jac1[species * ncomp + variable];
                jacobian[row * nvars + ncomp + variable] = -activity_jac2[species * ncomp + variable];
            }
            ++row;
        }

        std::vector<int> cations;
        std::vector<int> anions;
        for (std::size_t species = 0; species < ncomp; ++species) {
            if (charges[species] > 1.0e-12) {
                cations.push_back(static_cast<int>(species));
            } else if (charges[species] < -1.0e-12) {
                anions.push_back(static_cast<int>(species));
            }
        }
        for (int cation : cations) {
            for (int anion : anions) {
                const std::size_t c = static_cast<std::size_t>(cation);
                const std::size_t a = static_cast<std::size_t>(anion);
                const double cation_weight = std::abs(charges[a]);
                const double anion_weight = std::abs(charges[c]);
                for (std::size_t variable = 0; variable < ncomp; ++variable) {
                    jacobian[row * nvars + variable] =
                        cation_weight * activity_jac1[c * ncomp + variable]
                        + anion_weight * activity_jac1[a * ncomp + variable];
                    jacobian[row * nvars + ncomp + variable] =
                        -cation_weight * activity_jac2[c * ncomp + variable]
                        - anion_weight * activity_jac2[a * ncomp + variable];
                }
                ++row;
            }
        }
    }

    for (std::size_t species = 0; species < ncomp; ++species) {
        jacobian[row * nvars + species] = charges[species] * phase1.amounts[species];
    }
    ++row;
    for (std::size_t species = 0; species < ncomp; ++species) {
        jacobian[row * nvars + ncomp + species] = charges[species] * phase2.amounts[species];
    }
    return jacobian;
}

void add_phase_species_hessian(
    std::vector<double>& target,
    std::size_t residual_row,
    std::size_t variable_count,
    std::size_t phase_offset,
    std::size_t species_count,
    double coefficient,
    const std::vector<double>& species_hessian_tensor,
    std::size_t species
) {
    if (coefficient == 0.0) {
        return;
    }
    if (species_hessian_tensor.size() != species_count * species_count * species_count) {
        throw ValueError("reactive residual Hessian tensor shape does not match species count.");
    }
    for (std::size_t first = 0; first < species_count; ++first) {
        for (std::size_t second = 0; second < species_count; ++second) {
            target[residual_row * variable_count * variable_count
                + (phase_offset + first) * variable_count
                + phase_offset + second] += coefficient
                * species_hessian_tensor[species * species_count * species_count + first * species_count + second];
        }
    }
}

void add_phase_amount_diagonal_hessian(
    std::vector<double>& target,
    std::size_t residual_row,
    std::size_t variable_count,
    std::size_t phase_offset,
    const std::vector<double>& coefficients,
    const std::vector<double>& amounts
) {
    if (coefficients.size() != amounts.size()) {
        throw ValueError("reactive amount Hessian coefficients must match phase amounts.");
    }
    for (std::size_t species = 0; species < amounts.size(); ++species) {
        target[residual_row * variable_count * variable_count
            + (phase_offset + species) * variable_count
            + phase_offset + species] += coefficients[species] * amounts[species];
    }
}

std::vector<double> reactive_phase_residual_hessian_tensor_row_major(
    const std::shared_ptr<ePCSAFTMixtureNative>& phase1_mixture,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase2_mixture,
    double t,
    double p,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    const ReactivePhaseState& phase1,
    const ReactivePhaseState& phase2,
    const std::vector<double>& charges,
    const ReactivePhaseResidualEvaluationNative& residual_eval
) {
    const std::size_t ncomp = phase1.composition.size();
    const std::size_t nvars = 2 * ncomp;
    std::vector<double> tensor(residual_eval.residual.size() * nvars * nvars, 0.0);
    const std::vector<double> activity_hess1 =
        phase_ln_activity_log_amount_hessian_tensor(phase1_mixture, t, p, phase1);
    const std::vector<double> activity_hess2 =
        phase_ln_activity_log_amount_hessian_tensor(phase2_mixture, t, p, phase2);
    validate_reaction_standard_states(reaction_standard_states, reaction_rows);
    const bool phase_tagged_reactions = has_phase_tagged_reaction_stoichiometry(
        reaction_phase_stoichiometry_row_major,
        reaction_rows,
        static_cast<int>(ncomp)
    );
    std::size_t row = 0;

    for (int balance = 0; balance < balance_rows; ++balance) {
        std::vector<double> coefficients(ncomp, 0.0);
        for (std::size_t species = 0; species < ncomp; ++species) {
            coefficients[species] = balance_matrix_row_major[static_cast<std::size_t>(balance) * ncomp + species];
        }
        add_phase_amount_diagonal_hessian(tensor, row, nvars, 0, coefficients, phase1.amounts);
        add_phase_amount_diagonal_hessian(tensor, row, nvars, ncomp, coefficients, phase2.amounts);
        ++row;
    }

    if (phase_tagged_reactions) {
        for (int reaction = 0; reaction < reaction_rows; ++reaction) {
            const std::vector<double> reaction_hess1 = phase_reaction_standard_state_log_amount_hessian_tensor(
                phase1_mixture,
                t,
                p,
                phase1,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            const std::vector<double> reaction_hess2 = phase_reaction_standard_state_log_amount_hessian_tensor(
                phase2_mixture,
                t,
                p,
                phase2,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            const std::size_t reaction_offset = static_cast<std::size_t>(reaction) * 2 * ncomp;
            for (std::size_t species = 0; species < ncomp; ++species) {
                add_phase_species_hessian(
                    tensor,
                    row,
                    nvars,
                    0,
                    ncomp,
                    reaction_phase_stoichiometry_row_major[reaction_offset + species],
                    reaction_hess1,
                    species
                );
                add_phase_species_hessian(
                    tensor,
                    row,
                    nvars,
                    ncomp,
                    ncomp,
                    reaction_phase_stoichiometry_row_major[reaction_offset + ncomp + species],
                    reaction_hess2,
                    species
                );
            }
            ++row;
        }
    } else {
        for (int reaction = 0; reaction < reaction_rows; ++reaction) {
            const std::vector<double> reaction_hess1 = phase_reaction_standard_state_log_amount_hessian_tensor(
                phase1_mixture,
                t,
                p,
                phase1,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            for (std::size_t species = 0; species < ncomp; ++species) {
                add_phase_species_hessian(
                    tensor,
                    row,
                    nvars,
                    0,
                    ncomp,
                    reaction_stoichiometry_row_major[static_cast<std::size_t>(reaction) * ncomp + species],
                    reaction_hess1,
                    species
                );
            }
            ++row;
        }

        for (int reaction = 0; reaction < reaction_rows; ++reaction) {
            const std::vector<double> reaction_hess2 = phase_reaction_standard_state_log_amount_hessian_tensor(
                phase2_mixture,
                t,
                p,
                phase2,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            for (std::size_t species = 0; species < ncomp; ++species) {
                add_phase_species_hessian(
                    tensor,
                    row,
                    nvars,
                    ncomp,
                    ncomp,
                    reaction_stoichiometry_row_major[static_cast<std::size_t>(reaction) * ncomp + species],
                    reaction_hess2,
                    species
                );
            }
            ++row;
        }
    }

    if (!phase_tagged_reactions) {
        for (std::size_t species = 0; species < ncomp; ++species) {
            if (std::abs(charges[species]) > 1.0e-12) {
                continue;
            }
            add_phase_species_hessian(tensor, row, nvars, 0, ncomp, 1.0, activity_hess1, species);
            add_phase_species_hessian(tensor, row, nvars, ncomp, ncomp, -1.0, activity_hess2, species);
            ++row;
        }

        std::vector<int> cations;
        std::vector<int> anions;
        for (std::size_t species = 0; species < ncomp; ++species) {
            if (charges[species] > 1.0e-12) {
                cations.push_back(static_cast<int>(species));
            } else if (charges[species] < -1.0e-12) {
                anions.push_back(static_cast<int>(species));
            }
        }
        for (int cation : cations) {
            for (int anion : anions) {
                const std::size_t c = static_cast<std::size_t>(cation);
                const std::size_t a = static_cast<std::size_t>(anion);
                const double cation_weight = std::abs(charges[a]);
                const double anion_weight = std::abs(charges[c]);
                add_phase_species_hessian(tensor, row, nvars, 0, ncomp, cation_weight, activity_hess1, c);
                add_phase_species_hessian(tensor, row, nvars, 0, ncomp, anion_weight, activity_hess1, a);
                add_phase_species_hessian(tensor, row, nvars, ncomp, ncomp, -cation_weight, activity_hess2, c);
                add_phase_species_hessian(tensor, row, nvars, ncomp, ncomp, -anion_weight, activity_hess2, a);
                ++row;
            }
        }
    }

    add_phase_amount_diagonal_hessian(tensor, row, nvars, 0, charges, phase1.amounts);
    ++row;
    add_phase_amount_diagonal_hessian(tensor, row, nvars, ncomp, charges, phase2.amounts);
    ++row;
    if (row != residual_eval.residual.size()) {
        throw ValueError("reactive residual Hessian row count did not match the residual vector.");
    }
    return tensor;
}

void validate_jacobian_backend(const EquilibriumOptionsNative& options) {
    if (options.jacobian_backend == "auto"
        || options.jacobian_backend == "analytic"
        || options.jacobian_backend == "cppad") {
        return;
    }
    throw ValueError("reactive phase jacobian_backend must be auto, analytic, or cppad.");
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
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    const std::vector<double>& variables,
    bool has_variables,
    const std::vector<double>& initial_phase1,
    const std::vector<double>& initial_phase2,
    double initial_phase_fraction_phase2,
    bool has_initial_phases,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase1_mixture,
    const std::vector<int>& phase1_global_indices,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase2_mixture,
    const std::vector<int>& phase2_global_indices
) {
    const std::size_t ncomp = mixture->ncomp();
    const bool use_phase_models = phase1_mixture || phase2_mixture;
    if (use_phase_models) {
        validate_phase_model_indices(phase1_mixture, phase1_global_indices, ncomp, "phase1");
        validate_phase_model_indices(phase2_mixture, phase2_global_indices, ncomp, "phase2");
    }
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

    ReactivePhaseState phase1 = use_phase_models
        ? evaluate_scoped_phase_state(
            phase1_mixture,
            t,
            p,
            exp_amounts(eval_variables, 0, ncomp),
            options.min_composition,
            "reactive_phase_equilibrium_phase1",
            phase1_global_indices,
            ncomp
        )
        : evaluate_phase_state(
            mixture,
            t,
            p,
            exp_amounts(eval_variables, 0, ncomp),
            options.min_composition,
            "reactive_phase_equilibrium_phase1"
        );
    ReactivePhaseState phase2 = use_phase_models
        ? evaluate_scoped_phase_state(
            phase2_mixture,
            t,
            p,
            exp_amounts(eval_variables, ncomp, ncomp),
            options.min_composition,
            "reactive_phase_equilibrium_phase2",
            phase2_global_indices,
            ncomp
        )
        : evaluate_phase_state(
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
    const bool phase_tagged_reactions = has_phase_tagged_reaction_stoichiometry(
        reaction_phase_stoichiometry_row_major,
        reaction_rows,
        ncomp_int
    );

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
    if (phase_tagged_reactions) {
        out.reaction_residuals_cross_phase = cross_phase_reaction_residuals(
            reaction_phase_stoichiometry_row_major,
            reaction_rows,
            ncomp_int,
            log_equilibrium_constants,
            reaction_standard_states,
            phase1,
            phase2,
            options.min_composition
        );
    } else {
        out.reaction_residuals_phase1 = reaction_residuals(
            reaction_stoichiometry_row_major,
            reaction_rows,
            ncomp_int,
            log_equilibrium_constants,
            reaction_standard_states,
            phase1,
            options.min_composition
        );
        out.reaction_residuals_phase2 = reaction_residuals(
            reaction_stoichiometry_row_major,
            reaction_rows,
            ncomp_int,
            log_equilibrium_constants,
            reaction_standard_states,
            phase2,
            options.min_composition
        );
    }
    std::vector<double> charges = mixture->args().z;
    if (charges.size() != ncomp) {
        charges.assign(ncomp, 0.0);
    }
    out.phase_charge_residuals = {
        composition_charge(phase1.amounts, charges),
        composition_charge(phase2.amounts, charges),
    };
    if (!phase_tagged_reactions) {
        out.neutral_phase_equilibrium_residuals = neutral_phase_residuals(charges, ln_activity1, ln_activity2);
        out.ionic_equilibrium_residuals = ionic_phase_residuals(charges, ln_activity1, ln_activity2);
    }
    append_block(out.residual, out.element_balance_residuals);
    if (phase_tagged_reactions) {
        append_block(out.residual, out.reaction_residuals_cross_phase);
    } else {
        append_block(out.residual, out.reaction_residuals_phase1);
        append_block(out.residual, out.reaction_residuals_phase2);
    }
    append_block(out.residual, out.neutral_phase_equilibrium_residuals);
    append_block(out.residual, out.ionic_equilibrium_residuals);
    append_block(out.residual, out.phase_charge_residuals);
    for (double value : out.residual) {
        out.objective += 0.5 * value * value;
    }
    out.jacobian_rows = static_cast<int>(out.residual.size());
    out.jacobian_cols = static_cast<int>(eval_variables.size());
    out.phase_distance = phase_distance(phase1.composition, phase2.composition);
    validate_jacobian_backend(options);
    out.jacobian_row_major = reactive_phase_residual_jacobian_row_major(
        use_phase_models ? phase1_mixture : mixture,
        use_phase_models ? phase2_mixture : mixture,
        t,
        p,
        balance_matrix_row_major,
        balance_rows,
        reaction_stoichiometry_row_major,
        reaction_rows,
        reaction_standard_states,
        reaction_phase_stoichiometry_row_major,
        phase1,
        phase2,
        charges,
        out
    );
    out.residual_hessian_tensor_row_major = reactive_phase_residual_hessian_tensor_row_major(
        use_phase_models ? phase1_mixture : mixture,
        use_phase_models ? phase2_mixture : mixture,
        t,
        p,
        balance_matrix_row_major,
        balance_rows,
        reaction_stoichiometry_row_major,
        reaction_rows,
        reaction_standard_states,
        reaction_phase_stoichiometry_row_major,
        phase1,
        phase2,
        charges,
        out
    );
    out.gradient.assign(eval_variables.size(), 0.0);
    for (std::size_t row = 0; row < out.residual.size(); ++row) {
        for (std::size_t col = 0; col < eval_variables.size(); ++col) {
            out.gradient[col] += out.jacobian_row_major[row * eval_variables.size() + col] * out.residual[row];
        }
    }

    out.diagnostics_string["residual_surface"] = "native_reactive_phase_equilibrium_coupled_state";
    out.diagnostics_string["variable_model"] = out.variable_model;
    out.diagnostics_string["residual_blocks"] = phase_tagged_reactions
        ? "element_balance,phase_tagged_reaction_equilibrium,neutral_phase_equilibrium,ionic_equilibrium,phase_charge"
        : "element_balance,reaction_equilibrium,neutral_phase_equilibrium,ionic_equilibrium,phase_charge";
    out.diagnostics_string["jacobian_backend"] = "cppad_implicit";
    out.diagnostics_string["derivative_backend"] = "cppad_implicit";
    out.diagnostics_string["coupling_level"] = use_phase_models
        ? "phase_tagged_distinct_native_residual_states"
        : "single_native_residual_state";
    out.diagnostics_string["phase_model"] = use_phase_models
        ? "phase_models_aq_org"
        : "two_liquid_phases";
    out.diagnostics_string["reaction_residual_basis"] = reaction_standard_state_summary(reaction_standard_states);
    out.diagnostics_string["reaction_phase_scope"] = phase_tagged_reactions
        ? "phase_tagged_cross_phase"
        : "per_phase_same_stoichiometry";
    out.diagnostics_bool["jacobian_available"] = true;
    out.diagnostics_bool["derivative_available"] = true;
    out.diagnostics_bool["solved_state_sensitivity_available"] = true;
    out.diagnostics_bool["reaction_and_phase_residuals_share_state"] = true;
    out.diagnostics_bool["reaction_residual_standard_state_applied"] = true;
    out.diagnostics_bool["phase_tagged_reaction_stoichiometry"] = phase_tagged_reactions;
    out.diagnostics_bool["cross_phase_reaction_residuals"] = phase_tagged_reactions;
    out.diagnostics_bool["phase_models_applied"] = use_phase_models;
    out.diagnostics_bool["nonnegative_amounts_enforced_by_transform"] = true;
    out.diagnostics_bool["composition_normalization_enforced_by_transform"] = true;
    out.diagnostics_int["phase_count"] = 2;
    out.diagnostics_int["component_count"] = static_cast<int>(ncomp);
    out.diagnostics_int["reaction_count"] = reaction_rows;
    out.diagnostics_int["balance_row_count"] = balance_rows;
    out.diagnostics_int["variable_count"] = static_cast<int>(eval_variables.size());
    out.diagnostics_int["residual_size"] = static_cast<int>(out.residual.size());
    out.diagnostics_int["element_balance_residual_size"] = static_cast<int>(out.element_balance_residuals.size());
    out.diagnostics_int["reaction_residual_size_per_phase"] = static_cast<int>(out.reaction_residuals_phase1.size());
    out.diagnostics_int["cross_phase_reaction_residual_size"] = static_cast<int>(out.reaction_residuals_cross_phase.size());
    out.diagnostics_int["reaction_residual_size"] = phase_tagged_reactions
        ? static_cast<int>(out.reaction_residuals_cross_phase.size())
        : static_cast<int>(out.reaction_residuals_phase1.size() + out.reaction_residuals_phase2.size());
    out.diagnostics_int["neutral_phase_equilibrium_residual_size"] = static_cast<int>(out.neutral_phase_equilibrium_residuals.size());
    out.diagnostics_int["ionic_equilibrium_residual_size"] = static_cast<int>(out.ionic_equilibrium_residuals.size());
    out.diagnostics_int["phase_charge_residual_size"] = static_cast<int>(out.phase_charge_residuals.size());
    out.diagnostics_double["element_balance_norm"] = max_abs(out.element_balance_residuals);
    out.diagnostics_double["reaction_residual_norm"] = reaction_residual_norm(out);
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
    out.diagnostics_vector["reaction_residual_cross_phase"] = out.reaction_residuals_cross_phase;
    out.diagnostics_vector["neutral_phase_equilibrium_residual"] = out.neutral_phase_equilibrium_residuals;
    out.diagnostics_vector["ionic_equilibrium_residual"] = out.ionic_equilibrium_residuals;
    out.diagnostics_vector["phase_charge_residual"] = out.phase_charge_residuals;
    out.diagnostics_vector["reaction_standard_states"] = std::vector<double>(
        reaction_standard_states.begin(),
        reaction_standard_states.end()
    );
    out.diagnostics_vector["reaction_standard_state_codes"] = out.diagnostics_vector["reaction_standard_states"];
    return out;
}

namespace {

double clamped_log_amount(double amount, double floor) {
    return std::log(std::max(amount, floor));
}

std::vector<double> build_reactive_liquid_root_initial_variables(
    const std::vector<double>& feed,
    const std::vector<double>& charges,
    double floor,
    double shift_sign = 1.0
) {
    const std::size_t ncomp = feed.size();
    std::vector<double> phase1 = feed;
    std::vector<double> phase2 = feed;
    std::vector<std::size_t> neutral_indices;
    for (std::size_t i = 0; i < ncomp; ++i) {
        if (charges.size() == ncomp && std::abs(charges[i]) > 1.0e-12) {
            continue;
        }
        neutral_indices.push_back(i);
    }
    if (neutral_indices.size() >= 2) {
        const std::size_t first = neutral_indices.front();
        const std::size_t second = neutral_indices.back();
        const double shift = 0.20 * std::min(feed[first], feed[second]);
        if (shift > floor && feed[first] > shift + floor && feed[second] > shift + floor) {
            phase1[first] += shift_sign * shift;
            phase1[second] -= shift_sign * shift;
            phase2[first] -= shift_sign * shift;
            phase2[second] += shift_sign * shift;
        }
    }

    const double beta2 = 0.5;
    std::vector<double> out(2 * ncomp, 0.0);
    for (std::size_t i = 0; i < ncomp; ++i) {
        out[i] = clamped_log_amount((1.0 - beta2) * phase1[i], floor);
        out[ncomp + i] = clamped_log_amount(beta2 * phase2[i], floor);
    }
    return out;
}

ReactivePhaseState evaluate_scoped_phase_state(
    const std::shared_ptr<ePCSAFTMixtureNative>& phase_mixture,
    double t,
    double p,
    const std::vector<double>& amounts,
    double floor,
    const std::string& density_scope,
    const std::vector<int>& global_indices,
    std::size_t global_ncomp
) {
    validate_phase_model_indices(phase_mixture, global_indices, global_ncomp, density_scope);
    ReactivePhaseState out;
    out.amounts = amounts;
    out.amount_total = sum_amounts(amounts);
    out.global_component_count = global_ncomp;
    out.global_indices = global_indices;
    const std::vector<double> local_amounts = selected_amounts(amounts, global_indices, global_ncomp);
    const double local_total = sum_amounts(local_amounts);
    out.model_composition = composition_from_amounts(local_amounts, local_total, floor);
    out.composition.assign(global_ncomp, 0.0);
    for (std::size_t local = 0; local < global_indices.size(); ++local) {
        out.composition[static_cast<std::size_t>(global_indices[local])] = out.model_composition[local];
    }
    out.density = phase_mixture->solve_density_scoped(t, p, out.model_composition, 0, density_scope);
    std::shared_ptr<ePCSAFTStateNative> state =
        phase_mixture->state(t, out.model_composition, 0, false, 0.0, true, out.density);
    const std::vector<double> local_ln_phi = state->ln_fugacity_coefficient();
    if (local_ln_phi.size() != out.model_composition.size()) {
        throw ValueError("reactive phase model fugacity payload length mismatch.");
    }
    out.ln_phi.assign(global_ncomp, 0.0);
    for (std::size_t local = 0; local < global_indices.size(); ++local) {
        out.ln_phi[static_cast<std::size_t>(global_indices[local])] = local_ln_phi[local];
    }
    return out;
}

std::vector<double> build_reactive_phase_model_initial_variables(
    const std::vector<double>& feed,
    const std::vector<int>& phase1_indices,
    const std::vector<int>& phase2_indices,
    const std::vector<double>& charges,
    double floor,
    double shift_sign = 1.0
) {
    const std::size_t ncomp = feed.size();
    std::vector<double> out(2 * ncomp, std::log(floor));
    for (std::size_t species = 0; species < ncomp; ++species) {
        if (contains_index(phase1_indices, species)) {
            out[species] = clamped_log_amount(feed[species], floor);
        }
        if (contains_index(phase2_indices, species)) {
            out[ncomp + species] = clamped_log_amount(feed[species], floor);
        }
    }
    std::vector<std::size_t> shared_neutral_indices;
    for (std::size_t species = 0; species < ncomp; ++species) {
        const bool shared = contains_index(phase1_indices, species) && contains_index(phase2_indices, species);
        const bool charged = charges.size() == ncomp && std::abs(charges[species]) > 1.0e-12;
        if (shared && !charged) {
            shared_neutral_indices.push_back(species);
        }
    }
    if (shared_neutral_indices.size() >= 2) {
        const std::size_t first = shared_neutral_indices.front();
        const std::size_t second = shared_neutral_indices.back();
        const double shift = 0.20 * std::min(feed[first], feed[second]);
        if (shift > floor && feed[first] > shift + floor && feed[second] > shift + floor) {
            out[first] = clamped_log_amount(feed[first] + shift_sign * shift, floor);
            out[second] = clamped_log_amount(feed[second] - shift_sign * shift, floor);
            out[ncomp + first] = clamped_log_amount(feed[first] - shift_sign * shift, floor);
            out[ncomp + second] = clamped_log_amount(feed[second] + shift_sign * shift, floor);
        }
    }
    return out;
}

struct NamedInitialVariables {
    std::string seed_name;
    std::vector<double> variables;
};

std::vector<NamedInitialVariables> reactive_liquid_root_seed_candidates(
    const std::vector<double>& feed,
    const std::vector<double>& charges,
    double floor
) {
    return {
        {"canonical_shifted_feed", build_reactive_liquid_root_initial_variables(feed, charges, floor, 1.0)},
        {"mirrored_shifted_feed", build_reactive_liquid_root_initial_variables(feed, charges, floor, -1.0)},
    };
}

std::vector<NamedInitialVariables> reactive_phase_model_seed_candidates(
    const std::vector<double>& feed,
    const std::vector<int>& phase1_indices,
    const std::vector<int>& phase2_indices,
    const std::vector<double>& charges,
    double floor
) {
    std::vector<NamedInitialVariables> out;
    out.push_back({
        "canonical_shifted_feed",
        build_reactive_phase_model_initial_variables(feed, phase1_indices, phase2_indices, charges, floor, 1.0)
    });
    const std::vector<double> mirrored =
        build_reactive_phase_model_initial_variables(feed, phase1_indices, phase2_indices, charges, floor, -1.0);
    if (mirrored != out.front().variables) {
        out.push_back({"mirrored_shifted_feed", mirrored});
    }
    return out;
}

class ReactiveLiquidRootTwoPhaseProblem final : public epcsaft::native::equilibrium_nlp::NlpProblem {
public:
    ReactiveLiquidRootTwoPhaseProblem(
        std::shared_ptr<ePCSAFTMixtureNative> mixture,
        double temperature,
        double target_pressure,
        std::vector<double> raw_feed,
        EquilibriumOptionsNative options,
        std::vector<double> balance_matrix_row_major,
        int balance_rows,
        std::vector<double> total_vector,
        std::vector<double> reaction_stoichiometry_row_major,
        int reaction_rows,
        std::vector<double> log_equilibrium_constants,
        std::vector<int> reaction_standard_states,
        std::vector<double> reaction_phase_stoichiometry_row_major,
        double phase_distance_tolerance,
        std::shared_ptr<ePCSAFTMixtureNative> phase1_mixture = nullptr,
        std::vector<int> phase1_global_indices = {},
        std::shared_ptr<ePCSAFTMixtureNative> phase2_mixture = nullptr,
        std::vector<int> phase2_global_indices = {}
    )
        : mixture_(std::move(mixture)),
          phase1_mixture_(std::move(phase1_mixture)),
          phase2_mixture_(std::move(phase2_mixture)),
          temperature_(temperature),
          target_pressure_(target_pressure),
          options_(std::move(options)),
          balance_matrix_row_major_(std::move(balance_matrix_row_major)),
          balance_rows_(balance_rows),
          total_vector_(std::move(total_vector)),
          reaction_stoichiometry_row_major_(std::move(reaction_stoichiometry_row_major)),
          reaction_rows_(reaction_rows),
          log_equilibrium_constants_(std::move(log_equilibrium_constants)),
          reaction_standard_states_(std::move(reaction_standard_states)),
          reaction_phase_stoichiometry_row_major_(std::move(reaction_phase_stoichiometry_row_major)),
          minimum_phase_distance_(std::max(phase_distance_tolerance, 1.0e-8)) {
        if (!mixture_) {
            throw ValueError("Reactive liquid-root LLE NLP requires a native mixture.");
        }
        if (!phase1_mixture_) {
            phase1_mixture_ = mixture_;
        }
        if (!phase2_mixture_) {
            phase2_mixture_ = mixture_;
        }
        if (!std::isfinite(temperature_) || temperature_ <= 0.0 || !std::isfinite(target_pressure_) || target_pressure_ <= 0.0) {
            throw ValueError("Reactive liquid-root LLE NLP received invalid T/P specifications.");
        }
        feed_ = normalize_feed(raw_feed, mixture_->ncomp(), options_.min_composition, "reactive_phase_equilibrium");
        phase1_global_indices_ = phase1_global_indices.empty()
            ? full_component_indices(feed_.size())
            : std::move(phase1_global_indices);
        phase2_global_indices_ = phase2_global_indices.empty()
            ? full_component_indices(feed_.size())
            : std::move(phase2_global_indices);
        const bool distinct_phase_models = phase1_mixture_ != mixture_ || phase2_mixture_ != mixture_;
        validate_phase_model_indices(phase1_mixture_, phase1_global_indices_, feed_.size(), "phase1");
        validate_phase_model_indices(phase2_mixture_, phase2_global_indices_, feed_.size(), "phase2");
        validate_reaction_standard_states(reaction_standard_states_, reaction_rows_);
        const std::size_t ncomp = feed_.size();
        if (balance_rows_ <= 0) {
            throw ValueError("Reactive liquid-root LLE NLP requires at least one balance row.");
        }
        if (reaction_rows_ <= 0) {
            throw ValueError("Reactive liquid-root LLE NLP requires at least one reaction row.");
        }
        if (balance_matrix_row_major_.size() != static_cast<std::size_t>(balance_rows_) * ncomp) {
            throw ValueError("Reactive liquid-root LLE balance matrix has an invalid row-major size.");
        }
        if (total_vector_.size() != static_cast<std::size_t>(balance_rows_)) {
            throw ValueError("Reactive liquid-root LLE total vector length must match balance rows.");
        }
        if (reaction_stoichiometry_row_major_.size() != static_cast<std::size_t>(reaction_rows_) * ncomp) {
            throw ValueError("Reactive liquid-root LLE reaction matrix has an invalid row-major size.");
        }
        if (log_equilibrium_constants_.size() != static_cast<std::size_t>(reaction_rows_)) {
            throw ValueError("Reactive liquid-root LLE log equilibrium constants length must match reaction rows.");
        }
        if (!reaction_phase_stoichiometry_row_major_.empty()) {
            has_phase_tagged_reaction_stoichiometry(
                reaction_phase_stoichiometry_row_major_,
                reaction_rows_,
                static_cast<int>(ncomp)
            );
        }
        initial_variables_ = distinct_phase_models
            ? build_reactive_phase_model_initial_variables(
                feed_,
                phase1_global_indices_,
                phase2_global_indices_,
                mixture_->args().z,
                options_.min_composition
            )
            : build_reactive_liquid_root_initial_variables(
                feed_,
                mixture_->args().z,
                options_.min_composition
            );
        const std::vector<double> species_upper = species_amount_upper_bounds(
            balance_matrix_row_major_,
            balance_rows_,
            total_vector_,
            ncomp,
            options_.min_composition
        );
        variable_upper_bounds_.assign(2 * ncomp, std::log(1.0));
        for (std::size_t species = 0; species < ncomp; ++species) {
            const double upper = std::log(species_upper[species]);
            variable_upper_bounds_[species] = upper;
            variable_upper_bounds_[ncomp + species] = upper;
            if (distinct_phase_models && !contains_index(phase1_global_indices_, species)) {
                variable_upper_bounds_[species] = std::log(10.0 * options_.min_composition);
            }
            if (distinct_phase_models && !contains_index(phase2_global_indices_, species)) {
                variable_upper_bounds_[ncomp + species] = std::log(10.0 * options_.min_composition);
            }
        }
        if (minimum_phase_distance_ > 0.0) {
            const std::vector<double> first_amounts = exp_amounts(initial_variables_, 0, ncomp);
            const std::vector<double> second_amounts = exp_amounts(initial_variables_, ncomp, ncomp);
            const std::vector<double> first = composition_from_amounts(
                first_amounts,
                sum_amounts(first_amounts),
                options_.min_composition
            );
            const std::vector<double> second = composition_from_amounts(
                second_amounts,
                sum_amounts(second_amounts),
                options_.min_composition
            );
            double max_distance = 0.0;
            for (std::size_t species = 0; species < ncomp; ++species) {
                const double diff = first[species] - second[species];
                if (std::abs(diff) > max_distance) {
                    max_distance = std::abs(diff);
                    separation_species_index_ = static_cast<int>(species);
                    separation_sign_ = diff >= 0.0 ? 1.0 : -1.0;
                }
            }
            if (max_distance <= 0.0) {
                throw ValueError("Reactive liquid-root LLE NLP requires distinct initial phases for phase-separation gating.");
            }
        }
    }

    std::string name() const override {
        return "reactive_liquid_root_eos";
    }

    int variable_count() const override {
        return static_cast<int>(initial_variables_.size());
    }

    int constraint_count() const override {
        return balance_rows_
            + phase_tagged_reaction_constraint_count()
            + phase_charge_constraint_count()
            + separation_constraint_count();
    }

    int jacobian_nonzero_count() const override {
        int balance_nonzeros = 0;
        for (double coefficient : balance_matrix_row_major_) {
            if (coefficient != 0.0) {
                balance_nonzeros += 2;
            }
        }
        const int reaction_nonzeros = phase_tagged_reaction_constraint_count() * variable_count();
        int charge_nonzeros = 0;
        for (double charge : mixture_->args().z) {
            if (charge != 0.0) {
                charge_nonzeros += 2;
            }
        }
        const int separation_nonzeros = separation_constraint_count() > 0 ? 2 * static_cast<int>(feed_.size()) : 0;
        return balance_nonzeros + reaction_nonzeros + charge_nonzeros + separation_nonzeros;
    }

    epcsaft::native::equilibrium_nlp::NlpBounds bounds() const override {
        epcsaft::native::equilibrium_nlp::NlpBounds out;
        out.variable_lower.assign(initial_variables_.size(), std::log(options_.min_composition));
        out.variable_upper = variable_upper_bounds_;
        out.constraint_lower.assign(static_cast<std::size_t>(balance_rows_), 0.0);
        out.constraint_upper.assign(static_cast<std::size_t>(balance_rows_), 0.0);
        for (int row = 0; row < phase_tagged_reaction_constraint_count(); ++row) {
            out.constraint_lower.push_back(0.0);
            out.constraint_upper.push_back(0.0);
        }
        for (int row = 0; row < phase_charge_constraint_count(); ++row) {
            out.constraint_lower.push_back(0.0);
            out.constraint_upper.push_back(0.0);
        }
        if (separation_constraint_count() > 0) {
            out.constraint_lower.push_back(minimum_phase_distance_);
            out.constraint_upper.push_back(1.0e12);
        }
        return out;
    }

    std::vector<double> initial_point() const override {
        return initial_variables_;
    }

    double objective(const std::vector<double>& variables) const override {
        try {
            return evaluate_cached(variables).objective;
        } catch (const std::exception&) {
            return 1.0e100;
        }
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        try {
            return evaluate_cached(variables).gradient;
        } catch (const std::exception&) {
            return std::vector<double>(variables.size(), 0.0);
        }
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        const ReactivePhaseResidualEvaluationNative& eval = evaluate_cached(variables);
        std::vector<double> out = eval.element_balance_residuals;
        if (phase_tagged_reaction_constraint_count() > 0) {
            out.insert(out.end(), eval.reaction_residuals_cross_phase.begin(), eval.reaction_residuals_cross_phase.end());
        }
        if (phase_charge_constraint_count() > 0) {
            out.insert(out.end(), eval.phase_charge_residuals.begin(), eval.phase_charge_residuals.end());
        }
        if (separation_constraint_count() > 0) {
            out.push_back(phase_separation(variables));
        }
        return out;
    }

    epcsaft::native::equilibrium_nlp::NlpJacobianStructure jacobian_structure() const override {
        epcsaft::native::equilibrium_nlp::NlpJacobianStructure out;
        out.rows.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        out.cols.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < balance_rows_; ++row) {
            for (std::size_t species = 0; species < feed_.size(); ++species) {
                const double coefficient =
                    balance_matrix_row_major_[static_cast<std::size_t>(row) * feed_.size() + species];
                if (coefficient == 0.0) {
                    continue;
                }
                out.rows.push_back(row);
                out.cols.push_back(static_cast<int>(species));
                out.rows.push_back(row);
                out.cols.push_back(static_cast<int>(feed_.size() + species));
            }
        }
        for (int reaction_row = 0; reaction_row < phase_tagged_reaction_constraint_count(); ++reaction_row) {
            const int row = balance_rows_ + reaction_row;
            for (int col = 0; col < variable_count(); ++col) {
                out.rows.push_back(row);
                out.cols.push_back(col);
            }
        }
        for (int charge_row = 0; charge_row < phase_charge_constraint_count(); ++charge_row) {
            const int row = balance_rows_ + phase_tagged_reaction_constraint_count() + charge_row;
            const bool phase1 = charge_row == 0;
            const int phase_offset = phase1 ? 0 : static_cast<int>(feed_.size());
            for (std::size_t species = 0; species < feed_.size(); ++species) {
                const double charge = mixture_->args().z[species];
                if (charge == 0.0) {
                    continue;
                }
                out.rows.push_back(row);
                out.cols.push_back(phase_offset + static_cast<int>(species));
            }
        }
        if (separation_constraint_count() > 0) {
            const int row = balance_rows_
                + phase_tagged_reaction_constraint_count()
                + phase_charge_constraint_count();
            for (std::size_t species = 0; species < feed_.size(); ++species) {
                out.rows.push_back(row);
                out.cols.push_back(static_cast<int>(species));
                out.rows.push_back(row);
                out.cols.push_back(static_cast<int>(feed_.size() + species));
            }
        }
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        const std::size_t ncomp = feed_.size();
        const ReactivePhaseResidualEvaluationNative& eval = evaluate_cached(variables);
        std::vector<double> values;
        values.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < balance_rows_; ++row) {
            for (std::size_t species = 0; species < ncomp; ++species) {
                const double coefficient =
                    balance_matrix_row_major_[static_cast<std::size_t>(row) * ncomp + species];
                if (coefficient == 0.0) {
                    continue;
                }
                values.push_back(coefficient * eval.phase1_amounts[species]);
                values.push_back(coefficient * eval.phase2_amounts[species]);
            }
        }
        if (phase_tagged_reaction_constraint_count() > 0) {
            const std::size_t nvars = static_cast<std::size_t>(variable_count());
            for (int reaction_row = 0; reaction_row < phase_tagged_reaction_constraint_count(); ++reaction_row) {
                const std::size_t residual_row = static_cast<std::size_t>(balance_rows_ + reaction_row);
                values.insert(
                    values.end(),
                    eval.jacobian_row_major.begin() + static_cast<std::ptrdiff_t>(residual_row * nvars),
                    eval.jacobian_row_major.begin() + static_cast<std::ptrdiff_t>((residual_row + 1) * nvars)
                );
            }
        }
        if (phase_charge_constraint_count() > 0) {
            std::vector<double> charges = mixture_->args().z;
            if (charges.size() != ncomp) {
                charges.assign(ncomp, 0.0);
            }
            for (std::size_t species = 0; species < ncomp; ++species) {
                if (charges[species] != 0.0) {
                    values.push_back(charges[species] * eval.phase1_amounts[species]);
                }
            }
            for (std::size_t species = 0; species < ncomp; ++species) {
                if (charges[species] != 0.0) {
                    values.push_back(charges[species] * eval.phase2_amounts[species]);
                }
            }
        }
        if (separation_constraint_count() > 0) {
            const std::vector<double> separation_jacobian = phase_separation_jacobian(variables);
            for (std::size_t species = 0; species < ncomp; ++species) {
                values.push_back(separation_jacobian[species]);
                values.push_back(separation_jacobian[ncomp + species]);
            }
        }
        return values;
    }

    bool has_exact_hessian() const override {
        return true;
    }

    int hessian_nonzero_count() const override {
        return epcsaft::native::equilibrium_nlp::LagrangianHessianAssembler(variable_count()).nonzero_count();
    }

    epcsaft::native::equilibrium_nlp::NlpHessianStructure hessian_structure() const override {
        return epcsaft::native::equilibrium_nlp::LagrangianHessianAssembler(variable_count()).structure();
    }

    std::vector<double> hessian_values(
        const std::vector<double>& variables,
        double objective_factor,
        const std::vector<double>& constraint_multipliers
    ) const override {
        const ReactivePhaseResidualEvaluationNative& eval = evaluate_cached(variables);
        if (constraint_multipliers.size() != static_cast<std::size_t>(constraint_count())) {
            throw ValueError("Reactive liquid-root route Hessian multiplier vector size does not match constraints.");
        }
        epcsaft::native::equilibrium_nlp::ResidualSecondOrderData residuals;
        residuals.residual_count = eval.jacobian_rows;
        residuals.variable_count = eval.jacobian_cols;
        residuals.residuals = eval.residual;
        residuals.jacobian_row_major = eval.jacobian_row_major;
        residuals.residual_hessian_tensor_row_major = eval.residual_hessian_tensor_row_major;
        residuals.backend = "cppad_implicit_reactive_residual";
        epcsaft::native::equilibrium_nlp::ObjectiveSecondOrderData objective =
            epcsaft::native::equilibrium_nlp::least_squares_objective_second_order(residuals);

        const int nvars = variable_count();
        const std::size_t n = static_cast<std::size_t>(nvars);
        epcsaft::native::equilibrium_nlp::ConstraintSecondOrderData constraints_data;
        constraints_data.constraint_count = constraint_count();
        constraints_data.variable_count = nvars;
        constraints_data.values = constraints(variables);
        constraints_data.hessian_tensor_row_major.assign(
            static_cast<std::size_t>(constraints_data.constraint_count) * n * n,
            0.0
        );
        constraints_data.has_hessian.assign(static_cast<std::size_t>(constraints_data.constraint_count), false);
        constraints_data.backend = "cppad_implicit_reactive_residual";
        const auto copy_residual_hessian = [&](int residual_row, int constraint_row) {
            constraints_data.has_hessian[static_cast<std::size_t>(constraint_row)] = true;
            const std::size_t source_offset = static_cast<std::size_t>(residual_row) * n * n;
            const std::size_t target_offset = static_cast<std::size_t>(constraint_row) * n * n;
            std::copy(
                eval.residual_hessian_tensor_row_major.begin() + static_cast<std::ptrdiff_t>(source_offset),
                eval.residual_hessian_tensor_row_major.begin() + static_cast<std::ptrdiff_t>(source_offset + n * n),
                constraints_data.hessian_tensor_row_major.begin() + static_cast<std::ptrdiff_t>(target_offset)
            );
        };
        for (int row = 0; row < balance_rows_; ++row) {
            copy_residual_hessian(row, row);
        }
        for (int reaction_row = 0; reaction_row < phase_tagged_reaction_constraint_count(); ++reaction_row) {
            copy_residual_hessian(balance_rows_ + reaction_row, balance_rows_ + reaction_row);
        }
        if (phase_charge_constraint_count() > 0) {
            const int charge_residual_start = eval.jacobian_rows - 2;
            const int charge_constraint_start = balance_rows_ + phase_tagged_reaction_constraint_count();
            copy_residual_hessian(charge_residual_start, charge_constraint_start);
            copy_residual_hessian(charge_residual_start + 1, charge_constraint_start + 1);
        }
        if (separation_constraint_count() > 0) {
            const int row = constraint_count() - 1;
            constraints_data.has_hessian[static_cast<std::size_t>(row)] = true;
            const std::vector<double> separation = phase_separation_hessian(variables);
            const std::size_t offset = static_cast<std::size_t>(row) * n * n;
            std::copy(
                separation.begin(),
                separation.end(),
                constraints_data.hessian_tensor_row_major.begin() + static_cast<std::ptrdiff_t>(offset)
            );
        }
        return epcsaft::native::equilibrium_nlp::LagrangianHessianAssembler(nvars).values(
            objective_factor,
            objective,
            constraints_data,
            constraint_multipliers
        );
    }

    std::string hessian_backend() const override {
        return "cppad_implicit_reactive_residual";
    }

    epcsaft::native::equilibrium_nlp::NlpScaling scaling() const override {
        epcsaft::native::equilibrium_nlp::NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(initial_variables_.size(), 1.0);
        out.constraints.assign(static_cast<std::size_t>(constraint_count()), 1.0);
        return out;
    }

    std::map<std::string, std::string> diagnostics() const override {
        return {
            {"derivative_backend", "cppad_implicit"},
            {"density_backend", "liquid_pressure_root"},
            {"variable_model", "log_phase_species_amounts"},
        };
    }

    ReactivePhaseResidualEvaluationNative evaluate(const std::vector<double>& variables) const {
        return evaluate_cached(variables);
    }

    const ReactivePhaseResidualEvaluationNative& evaluate_cached(const std::vector<double>& variables) const {
        if (cached_evaluation_valid_ && cached_variables_ == variables) {
            return cached_evaluation_;
        }
        cached_evaluation_ = evaluate_reactive_phase_equilibrium_residual_native(
            mixture_,
            temperature_,
            target_pressure_,
            feed_,
            options_,
            balance_matrix_row_major_,
            balance_rows_,
            total_vector_,
            reaction_stoichiometry_row_major_,
            reaction_rows_,
            log_equilibrium_constants_,
            reaction_standard_states_,
            reaction_phase_stoichiometry_row_major_,
            variables,
            true,
            {},
            {},
            0.5,
            false,
            phase1_mixture_,
            phase1_global_indices_,
            phase2_mixture_,
            phase2_global_indices_
        );
        cached_variables_ = variables;
        cached_evaluation_valid_ = true;
        return cached_evaluation_;
    }

    const std::vector<double>& feed() const {
        return feed_;
    }

    int balance_row_count() const {
        return balance_rows_;
    }

    int reaction_count() const {
        return reaction_rows_;
    }

    double minimum_phase_distance() const {
        return minimum_phase_distance_;
    }

private:
    int phase_tagged_reaction_constraint_count() const {
        if (reaction_phase_stoichiometry_row_major_.empty()) {
            return 0;
        }
        return reaction_rows_;
    }

    int phase_charge_constraint_count() const {
        std::vector<double> charges = mixture_->args().z;
        if (charges.size() != feed_.size()) {
            return 0;
        }
        for (double charge : charges) {
            if (std::abs(charge) > 1.0e-12) {
                return 2;
            }
        }
        return 0;
    }

    int separation_constraint_count() const {
        return minimum_phase_distance_ > 0.0 ? 1 : 0;
    }

    double phase_separation(const std::vector<double>& variables) const {
        const std::size_t ncomp = feed_.size();
        const std::vector<double> first_amounts = exp_amounts(variables, 0, ncomp);
        const std::vector<double> second_amounts = exp_amounts(variables, ncomp, ncomp);
        const std::vector<double> first = composition_from_amounts(
            first_amounts,
            sum_amounts(first_amounts),
            options_.min_composition
        );
        const std::vector<double> second = composition_from_amounts(
            second_amounts,
            sum_amounts(second_amounts),
            options_.min_composition
        );
        const std::size_t selected = static_cast<std::size_t>(separation_species_index_);
        return separation_sign_ * (first[selected] - second[selected]);
    }

    std::vector<double> phase_separation_jacobian(const std::vector<double>& variables) const {
        const std::size_t ncomp = feed_.size();
        const std::vector<double> first_amounts = exp_amounts(variables, 0, ncomp);
        const std::vector<double> second_amounts = exp_amounts(variables, ncomp, ncomp);
        const std::vector<double> first = composition_from_amounts(
            first_amounts,
            sum_amounts(first_amounts),
            options_.min_composition
        );
        const std::vector<double> second = composition_from_amounts(
            second_amounts,
            sum_amounts(second_amounts),
            options_.min_composition
        );

        std::vector<double> row(static_cast<std::size_t>(variable_count()), 0.0);
        const std::size_t selected = static_cast<std::size_t>(separation_species_index_);
        for (std::size_t species = 0; species < ncomp; ++species) {
            const double first_indicator = species == selected ? 1.0 : 0.0;
            row[species] = separation_sign_ * first[selected] * (first_indicator - first[species]);

            const std::size_t second_offset = ncomp + species;
            const double second_indicator = species == selected ? 1.0 : 0.0;
            row[second_offset] = -separation_sign_ * second[selected] * (second_indicator - second[species]);
        }
        return row;
    }

    std::vector<double> phase_separation_hessian(const std::vector<double>& variables) const {
        const std::size_t ncomp = feed_.size();
        const std::size_t nvars = static_cast<std::size_t>(variable_count());
        const std::vector<double> first_amounts = exp_amounts(variables, 0, ncomp);
        const std::vector<double> second_amounts = exp_amounts(variables, ncomp, ncomp);
        const std::vector<double> first = composition_from_amounts(
            first_amounts,
            sum_amounts(first_amounts),
            options_.min_composition
        );
        const std::vector<double> second = composition_from_amounts(
            second_amounts,
            sum_amounts(second_amounts),
            options_.min_composition
        );
        std::vector<double> out(nvars * nvars, 0.0);
        const std::size_t selected = static_cast<std::size_t>(separation_species_index_);
        const auto add_block = [&](std::size_t offset, const std::vector<double>& composition, double scale) {
            const double x_selected = composition[selected];
            for (std::size_t first_index = 0; first_index < ncomp; ++first_index) {
                const double first_delta = selected == first_index ? 1.0 : 0.0;
                const double x_first = composition[first_index];
                for (std::size_t second_index = 0; second_index < ncomp; ++second_index) {
                    const double second_delta = selected == second_index ? 1.0 : 0.0;
                    const double cross_delta = first_index == second_index ? 1.0 : 0.0;
                    const double x_second = composition[second_index];
                    const double value =
                        x_selected * (second_delta - x_second) * (first_delta - x_first)
                        - x_selected * x_first * (cross_delta - x_second);
                    out[(offset + first_index) * nvars + offset + second_index] += scale * value;
                }
            }
        };
        add_block(0, first, separation_sign_);
        add_block(ncomp, second, -separation_sign_);
        return out;
    }

    std::shared_ptr<ePCSAFTMixtureNative> mixture_;
    std::shared_ptr<ePCSAFTMixtureNative> phase1_mixture_;
    std::shared_ptr<ePCSAFTMixtureNative> phase2_mixture_;
    double temperature_ = 0.0;
    double target_pressure_ = 0.0;
    EquilibriumOptionsNative options_;
    std::vector<double> feed_;
    std::vector<double> balance_matrix_row_major_;
    int balance_rows_ = 0;
    std::vector<double> total_vector_;
    std::vector<double> reaction_stoichiometry_row_major_;
    int reaction_rows_ = 0;
    std::vector<double> log_equilibrium_constants_;
    std::vector<int> reaction_standard_states_;
    std::vector<double> reaction_phase_stoichiometry_row_major_;
    std::vector<int> phase1_global_indices_;
    std::vector<int> phase2_global_indices_;
    std::vector<double> initial_variables_;
    std::vector<double> variable_upper_bounds_;
    mutable bool cached_evaluation_valid_ = false;
    mutable std::vector<double> cached_variables_;
    mutable ReactivePhaseResidualEvaluationNative cached_evaluation_;
    double minimum_phase_distance_ = 1.0e-8;
    int separation_species_index_ = 0;
    double separation_sign_ = 1.0;
};

epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosNlpContract reactive_liquid_root_contract_from_problem(
    const ReactiveLiquidRootTwoPhaseProblem& problem
) {
    epcsaft::native::equilibrium_nlp::validate_nlp_problem_shape(problem);
    const std::vector<double> initial = problem.initial_point();
    const epcsaft::native::equilibrium_nlp::NlpBounds bounds = problem.bounds();
    const epcsaft::native::equilibrium_nlp::NlpJacobianStructure structure = problem.jacobian_structure();
    const ReactivePhaseResidualEvaluationNative eval = problem.evaluate(initial);
    const std::vector<double> constraints = problem.constraints(initial);
    const std::vector<double> jacobian = problem.jacobian_values(initial);

    epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosNlpContract out;
    out.problem_name = problem.name();
    out.derivative_backend = "cppad_implicit";
    out.variable_model = "log_phase_species_amounts";
    out.density_backend = "liquid_pressure_root";
    out.phase_count = 2;
    out.species_count = static_cast<int>(problem.feed().size());
    out.balance_row_count = problem.balance_row_count();
    out.reaction_count = problem.reaction_count();
    out.variable_count = problem.variable_count();
    out.constraint_count = problem.constraint_count();
    out.jacobian_nonzero_count = problem.jacobian_nonzero_count();
    out.initial_point = initial;
    out.variable_lower_bounds = bounds.variable_lower;
    out.variable_upper_bounds = bounds.variable_upper;
    out.constraint_lower_bounds = bounds.constraint_lower;
    out.constraint_upper_bounds = bounds.constraint_upper;
    out.objective_at_initial = eval.objective;
    out.gradient_at_initial = eval.gradient;
    out.constraints_at_initial = constraints;
    out.jacobian_rows = structure.rows;
    out.jacobian_cols = structure.cols;
    out.jacobian_values_at_initial = jacobian;
    return out;
}

std::vector<std::vector<double>> phase_amounts_from_reactive_eval(
    const ReactivePhaseResidualEvaluationNative& eval
) {
    return {eval.phase1_amounts, eval.phase2_amounts};
}

epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosPostsolve reactive_liquid_root_postsolve_from_eval(
    const ReactiveLiquidRootTwoPhaseProblem& problem,
    const ReactivePhaseResidualEvaluationNative& eval,
    double conserved_balance_tolerance,
    double reaction_tolerance,
    double phase_equilibrium_tolerance,
    double phase_distance_tolerance
) {
    epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosPostsolve out;
    out.derivative_backend = "cppad_implicit";
    out.density_backend = "liquid_pressure_root";
    out.phase_count = 2;
    out.species_count = static_cast<int>(problem.feed().size());
    out.balance_row_count = problem.balance_row_count();
    out.reaction_count = problem.reaction_count();
    out.objective = eval.objective;
    out.constraints = eval.residual;
    out.reaction_stationarity_residuals = eval.reaction_residuals_cross_phase.empty()
        ? eval.reaction_residuals_phase1
        : eval.reaction_residuals_cross_phase;
    if (eval.reaction_residuals_cross_phase.empty()) {
        out.reaction_stationarity_residuals.insert(
            out.reaction_stationarity_residuals.end(),
            eval.reaction_residuals_phase2.begin(),
            eval.reaction_residuals_phase2.end()
        );
    }
    out.phase_equilibrium_residuals = eval.neutral_phase_equilibrium_residuals;
    out.phase_equilibrium_residuals.insert(
        out.phase_equilibrium_residuals.end(),
        eval.ionic_equilibrium_residuals.begin(),
        eval.ionic_equilibrium_residuals.end()
    );
    out.phase_charge_residuals = eval.phase_charge_residuals;
    out.phase_amount_totals = {
        std::accumulate(eval.phase1_amounts.begin(), eval.phase1_amounts.end(), 0.0),
        std::accumulate(eval.phase2_amounts.begin(), eval.phase2_amounts.end(), 0.0),
    };
    out.phase_densities = {eval.phase1_density, eval.phase2_density};
    out.phase_volumes = {
        out.phase_amount_totals[0] / eval.phase1_density,
        out.phase_amount_totals[1] / eval.phase2_density,
    };
    out.phase_compositions = {eval.phase1_composition, eval.phase2_composition};
    out.phase_ln_fugacity_coefficients = {
        eval.phase1_ln_fugacity_coefficient,
        eval.phase2_ln_fugacity_coefficient,
    };
    out.conserved_balance_norm = max_abs(eval.element_balance_residuals);
    out.charge_balance_norm = max_abs(eval.phase_charge_residuals);
    out.phase_equilibrium_norm = max_abs(out.phase_equilibrium_residuals);
    out.pressure_consistency_norm = 0.0;
    out.reaction_stationarity_norm = reaction_residual_norm(eval);
    out.phase_distance = eval.phase_distance;
    const double effective_phase_distance = std::max(
        phase_distance_tolerance,
        problem.minimum_phase_distance()
    );
    const bool finite_liquid_densities = std::isfinite(eval.phase1_density)
        && std::isfinite(eval.phase2_density)
        && eval.phase1_density > 0.0
        && eval.phase2_density > 0.0;
    const bool retained_phase_split = out.phase_amount_totals[0] > 0.0 && out.phase_amount_totals[1] > 0.0;
    out.accepted = out.conserved_balance_norm <= conserved_balance_tolerance
        && out.charge_balance_norm <= conserved_balance_tolerance
        && out.reaction_stationarity_norm <= reaction_tolerance
        && out.phase_equilibrium_norm <= phase_equilibrium_tolerance
        && out.phase_distance >= effective_phase_distance
        && finite_liquid_densities
        && retained_phase_split;
    if (out.accepted) {
        out.rejection_reason = "accepted";
    } else if (out.conserved_balance_norm > conserved_balance_tolerance) {
        out.rejection_reason = "conserved_balance";
    } else if (out.charge_balance_norm > conserved_balance_tolerance) {
        out.rejection_reason = "charge_balance";
    } else if (out.reaction_stationarity_norm > reaction_tolerance) {
        out.rejection_reason = "reaction_equilibrium";
    } else if (out.phase_equilibrium_norm > phase_equilibrium_tolerance) {
        out.rejection_reason = "phase_equilibrium";
    } else if (!finite_liquid_densities) {
        out.rejection_reason = "liquid_density";
    } else if (!retained_phase_split) {
        out.rejection_reason = "phase_fraction";
    } else {
        out.rejection_reason = "phase_distance";
    }
    return out;
}

int reactive_route_quality(const epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult& result) {
    if (result.accepted) {
        return 3;
    }
    if (result.solver_accepted) {
        return 2;
    }
    if (result.ran) {
        return 1;
    }
    return 0;
}

epcsaft::native::equilibrium_nlp::RouteSeedAttempt reactive_seed_attempt_from_result(
    const epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult& result
) {
    epcsaft::native::equilibrium_nlp::RouteSeedAttempt out;
    out.seed_name = result.seed_name;
    out.status = result.status;
    out.solver_status = result.solver_status;
    out.application_status = result.application_status;
    out.solver_accepted = result.solver_accepted;
    out.accepted = result.accepted;
    out.iteration_count = result.iteration_count;
    out.objective = result.objective;
    out.phase_distance = result.postsolve.phase_distance;
    out.conserved_balance_norm = result.postsolve.conserved_balance_norm;
    out.charge_balance_norm = result.postsolve.charge_balance_norm;
    out.phase_equilibrium_norm = result.postsolve.phase_equilibrium_norm;
    out.reaction_stationarity_norm = result.postsolve.reaction_stationarity_norm;
    return out;
}

}  // namespace

epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosNlpContract evaluate_reactive_phase_liquid_root_nlp_contract_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& options,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    double phase_distance_tolerance
) {
    const ReactiveLiquidRootTwoPhaseProblem problem(
        mixture,
        t,
        p,
        feed,
        options,
        balance_matrix_row_major,
        balance_rows,
        total_vector,
        reaction_stoichiometry_row_major,
        reaction_rows,
        log_equilibrium_constants,
        reaction_standard_states,
        reaction_phase_stoichiometry_row_major,
        phase_distance_tolerance
    );
    return reactive_liquid_root_contract_from_problem(problem);
}

epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult solve_reactive_phase_liquid_root_route_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& equilibrium_options,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    const epcsaft::native::equilibrium_nlp::IpoptSolveOptions& solve_options,
    double conserved_balance_tolerance,
    double reaction_tolerance,
    double phase_equilibrium_tolerance,
    double phase_distance_tolerance
) {
    const epcsaft::native::equilibrium_nlp::IpoptAdapterInfo adapter =
        epcsaft::native::equilibrium_nlp::native_ipopt_adapter_info();
    epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult best;
    best.compiled = adapter.compiled;
    best.adapter_available = adapter.adapter_available;
    best.adapter_kind = adapter.adapter_kind;
    best.exact_gradient_required = adapter.exact_gradient_required;
    best.exact_jacobian_required = adapter.exact_jacobian_required;
    best.problem_name = "reactive_liquid_root_eos";
    best.derivative_backend = "cppad_implicit";
    best.postsolve.derivative_backend = "cppad_implicit";
    best.postsolve.density_backend = "liquid_pressure_root";
    if (!adapter.compiled) {
        best.status = "ipopt_dependency_required";
        return best;
    }

    const ReactiveLiquidRootTwoPhaseProblem problem(
        mixture,
        t,
        p,
        feed,
        equilibrium_options,
        balance_matrix_row_major,
        balance_rows,
        total_vector,
        reaction_stoichiometry_row_major,
        reaction_rows,
        log_equilibrium_constants,
        reaction_standard_states,
        reaction_phase_stoichiometry_row_major,
        phase_distance_tolerance
    );
    best.phase_count = 2;
    best.species_count = static_cast<int>(problem.feed().size());
    best.balance_row_count = problem.balance_row_count();
    best.reaction_count = problem.reaction_count();
    const std::vector<NamedInitialVariables> seeds = reactive_liquid_root_seed_candidates(
        problem.feed(),
        mixture->args().z,
        equilibrium_options.min_composition
    );
    bool have_best = false;
    std::vector<epcsaft::native::equilibrium_nlp::RouteSeedAttempt> attempts;
    attempts.reserve(seeds.size() + (solve_options.initial_variables.empty() ? 0 : 1));

    auto run_attempt = [&](
        const std::string& seed_name,
        const epcsaft::native::equilibrium_nlp::IpoptSolveOptions& attempt_options
    ) {
        const epcsaft::native::equilibrium_nlp::IpoptSolveResult solve =
            epcsaft::native::equilibrium_nlp::solve_ipopt_nlp(problem, attempt_options);
        epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult result;
        result.compiled = adapter.compiled;
        result.adapter_available = adapter.adapter_available;
        result.adapter_kind = adapter.adapter_kind;
        result.exact_gradient_required = adapter.exact_gradient_required;
        result.exact_jacobian_required = adapter.exact_jacobian_required;
        result.problem_name = "reactive_liquid_root_eos";
        result.derivative_backend = "cppad_implicit";
        result.postsolve.derivative_backend = "cppad_implicit";
        result.postsolve.density_backend = "liquid_pressure_root";
        result.initial_point_strategy = "deterministic_multistart";
        result.seed_name = seed_name;
        result.phase_count = 2;
        result.species_count = static_cast<int>(problem.feed().size());
        result.balance_row_count = problem.balance_row_count();
        result.reaction_count = problem.reaction_count();
        result.ran = solve.solver_ran;
        result.solver_accepted = solve.accepted;
        result.accepted = solve.accepted;
        result.solver_status = solve.solver_status;
        result.application_status = solve.application_status;
        epcsaft::native::equilibrium_nlp::apply_ipopt_solve_metadata(result, solve);
        result.objective = solve.objective;
        result.variables = solve.variables;
        result.constraints = solve.constraints;
        result.status = solve.accepted ? "accepted" : "solver_rejected";
        if (!solve.accepted) {
            attempts.push_back(reactive_seed_attempt_from_result(result));
            if (!have_best || reactive_route_quality(result) > reactive_route_quality(best)) {
                best = result;
                have_best = true;
            }
            return result;
        }

        const ReactivePhaseResidualEvaluationNative eval = problem.evaluate(solve.variables);
        result.phase_amounts = phase_amounts_from_reactive_eval(eval);
        result.postsolve = reactive_liquid_root_postsolve_from_eval(
            problem,
            eval,
            conserved_balance_tolerance,
            reaction_tolerance,
            phase_equilibrium_tolerance,
            phase_distance_tolerance
        );
        result.phase_volumes = result.postsolve.phase_volumes;
        result.accepted = result.postsolve.accepted;
        result.status = result.accepted ? "accepted" : "postsolve_rejected";
        attempts.push_back(reactive_seed_attempt_from_result(result));
        if (!have_best || reactive_route_quality(result) > reactive_route_quality(best)) {
            best = result;
            have_best = true;
        }
        return result;
    };

    if (!solve_options.initial_variables.empty()) {
        const auto continuation = run_attempt("continuation_state", solve_options);
        if (continuation.accepted) {
            best.seed_attempts = attempts;
            return best;
        }
    }

    for (const auto& seed : seeds) {
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions attempt_options = solve_options;
        attempt_options.initial_variables = seed.variables;
        attempt_options.initial_bound_lower_multipliers.clear();
        attempt_options.initial_bound_upper_multipliers.clear();
        attempt_options.initial_constraint_multipliers.clear();
        const auto attempt = run_attempt(seed.seed_name, attempt_options);
        if (attempt.accepted) {
            break;
        }
    }

    best.initial_point_strategy = "deterministic_multistart";
    best.seed_attempts = attempts;
    return best;
}

epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult solve_reactive_phase_liquid_root_phase_model_route_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase1_mixture,
    const std::shared_ptr<ePCSAFTMixtureNative>& phase2_mixture,
    const std::vector<int>& phase1_global_indices,
    const std::vector<int>& phase2_global_indices,
    double t,
    double p,
    const std::vector<double>& feed,
    const EquilibriumOptionsNative& equilibrium_options,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const std::vector<double>& reaction_phase_stoichiometry_row_major,
    const epcsaft::native::equilibrium_nlp::IpoptSolveOptions& solve_options,
    double conserved_balance_tolerance,
    double reaction_tolerance,
    double phase_equilibrium_tolerance,
    double phase_distance_tolerance
) {
    const epcsaft::native::equilibrium_nlp::IpoptAdapterInfo adapter =
        epcsaft::native::equilibrium_nlp::native_ipopt_adapter_info();
    epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult best;
    best.compiled = adapter.compiled;
    best.adapter_available = adapter.adapter_available;
    best.adapter_kind = adapter.adapter_kind;
    best.exact_gradient_required = adapter.exact_gradient_required;
    best.exact_jacobian_required = adapter.exact_jacobian_required;
    best.problem_name = "reactive_liquid_root_eos";
    best.derivative_backend = "cppad_implicit";
    best.postsolve.derivative_backend = "cppad_implicit";
    best.postsolve.density_backend = "liquid_pressure_root";
    if (!adapter.compiled) {
        best.status = "ipopt_dependency_required";
        return best;
    }

    const ReactiveLiquidRootTwoPhaseProblem problem(
        mixture,
        t,
        p,
        feed,
        equilibrium_options,
        balance_matrix_row_major,
        balance_rows,
        total_vector,
        reaction_stoichiometry_row_major,
        reaction_rows,
        log_equilibrium_constants,
        reaction_standard_states,
        reaction_phase_stoichiometry_row_major,
        phase_distance_tolerance,
        phase1_mixture,
        phase1_global_indices,
        phase2_mixture,
        phase2_global_indices
    );
    best.phase_count = 2;
    best.species_count = static_cast<int>(problem.feed().size());
    best.balance_row_count = problem.balance_row_count();
    best.reaction_count = problem.reaction_count();
    const std::vector<NamedInitialVariables> seeds = reactive_phase_model_seed_candidates(
        problem.feed(),
        phase1_global_indices,
        phase2_global_indices,
        mixture->args().z,
        equilibrium_options.min_composition
    );
    bool have_best = false;
    std::vector<epcsaft::native::equilibrium_nlp::RouteSeedAttempt> attempts;
    attempts.reserve(seeds.size() + (solve_options.initial_variables.empty() ? 0 : 1));

    auto run_attempt = [&](
        const std::string& seed_name,
        const epcsaft::native::equilibrium_nlp::IpoptSolveOptions& attempt_options
    ) {
        const epcsaft::native::equilibrium_nlp::IpoptSolveResult solve =
            epcsaft::native::equilibrium_nlp::solve_ipopt_nlp(problem, attempt_options);
        epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult result;
        result.compiled = adapter.compiled;
        result.adapter_available = adapter.adapter_available;
        result.adapter_kind = adapter.adapter_kind;
        result.exact_gradient_required = adapter.exact_gradient_required;
        result.exact_jacobian_required = adapter.exact_jacobian_required;
        result.problem_name = "reactive_liquid_root_eos";
        result.derivative_backend = "cppad_implicit";
        result.postsolve.derivative_backend = "cppad_implicit";
        result.postsolve.density_backend = "liquid_pressure_root";
        result.initial_point_strategy = "deterministic_multistart";
        result.seed_name = seed_name;
        result.phase_count = 2;
        result.species_count = static_cast<int>(problem.feed().size());
        result.balance_row_count = problem.balance_row_count();
        result.reaction_count = problem.reaction_count();
        result.ran = solve.solver_ran;
        result.solver_accepted = solve.accepted;
        result.accepted = solve.accepted;
        result.solver_status = solve.solver_status;
        result.application_status = solve.application_status;
        epcsaft::native::equilibrium_nlp::apply_ipopt_solve_metadata(result, solve);
        result.objective = solve.objective;
        result.variables = solve.variables;
        result.constraints = solve.constraints;
        result.status = solve.accepted ? "accepted" : "solver_rejected";
        if (!solve.accepted) {
            attempts.push_back(reactive_seed_attempt_from_result(result));
            if (!have_best || reactive_route_quality(result) > reactive_route_quality(best)) {
                best = result;
                have_best = true;
            }
            return result;
        }

        const ReactivePhaseResidualEvaluationNative eval = problem.evaluate(solve.variables);
        result.phase_amounts = phase_amounts_from_reactive_eval(eval);
        result.postsolve = reactive_liquid_root_postsolve_from_eval(
            problem,
            eval,
            conserved_balance_tolerance,
            reaction_tolerance,
            phase_equilibrium_tolerance,
            phase_distance_tolerance
        );
        result.phase_volumes = result.postsolve.phase_volumes;
        result.accepted = result.postsolve.accepted;
        result.status = result.accepted ? "accepted" : "postsolve_rejected";
        attempts.push_back(reactive_seed_attempt_from_result(result));
        if (!have_best || reactive_route_quality(result) > reactive_route_quality(best)) {
            best = result;
            have_best = true;
        }
        return result;
    };

    if (!solve_options.initial_variables.empty()) {
        const auto continuation = run_attempt("continuation_state", solve_options);
        if (continuation.accepted) {
            best.seed_attempts = attempts;
            return best;
        }
    }

    for (const auto& seed : seeds) {
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions attempt_options = solve_options;
        attempt_options.initial_variables = seed.variables;
        attempt_options.initial_bound_lower_multipliers.clear();
        attempt_options.initial_bound_upper_multipliers.clear();
        attempt_options.initial_constraint_multipliers.clear();
        const auto attempt = run_attempt(seed.seed_name, attempt_options);
        if (attempt.accepted) {
            break;
        }
    }

    best.initial_point_strategy = "deterministic_multistart";
    best.seed_attempts = attempts;
    return best;
}
