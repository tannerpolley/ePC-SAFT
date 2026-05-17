#include "epcsaft_equilibrium.h"

#include "equilibrium_helpers.h"

#include <algorithm>
#include <cmath>
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
    const std::size_t ncomp = state.composition.size();
    std::vector<double> jacobian(ncomp * ncomp, 0.0);
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t variable = 0; variable < ncomp; ++variable) {
            jacobian[species * ncomp + variable] =
                (species == variable ? 1.0 : 0.0) - state.composition[variable];
        }
    }
    return jacobian;
}

std::vector<double> phase_ln_activity_log_amount_jacobian(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ReactivePhaseState& state
) {
    const std::size_t ncomp = state.composition.size();
    PhaseStateCompositionSensitivityResult sensitivity =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, state.composition, 0, mixture->args());
    if (!sensitivity.supported || sensitivity.jacobian_row_major.size() != ncomp * ncomp) {
        throw ValueError("reactive phase residual Jacobian requires supported phase-state fugacity sensitivities.");
    }
    std::vector<double> jacobian(ncomp * ncomp, 0.0);
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t variable = 0; variable < ncomp; ++variable) {
            const double dlogx = (species == variable ? 1.0 : 0.0) - state.composition[variable];
            double dlnphi = 0.0;
            for (std::size_t k = 0; k < ncomp; ++k) {
                const double dxk_dlogn =
                    state.composition[k] * ((k == variable ? 1.0 : 0.0) - state.composition[variable]);
                dlnphi += sensitivity.jacobian_row_major[species * ncomp + k] * dxk_dlogn;
            }
            jacobian[species * ncomp + variable] = dlogx + dlnphi;
        }
    }
    return jacobian;
}

std::vector<double> phase_log_concentration_log_amount_jacobian(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ReactivePhaseState& state
) {
    const std::size_t ncomp = state.composition.size();
    PhaseStateCompositionSensitivityResult sensitivity =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, state.composition, 0, mixture->args());
    if (!sensitivity.supported || sensitivity.density_composition_derivative.size() != ncomp) {
        throw ValueError("concentration reaction residual Jacobian requires supported density composition sensitivities.");
    }
    const double density = sensitivity.density > 0.0 ? sensitivity.density : state.density;
    if (!(std::isfinite(density) && density > 0.0)) {
        throw ValueError("concentration reaction residual Jacobian requires a finite positive molar density.");
    }
    std::vector<double> jacobian = phase_log_mole_fraction_log_amount_jacobian(state);
    std::vector<double> dlogrho_dlogn(ncomp, 0.0);
    for (std::size_t variable = 0; variable < ncomp; ++variable) {
        double drho_dlogn = 0.0;
        for (std::size_t k = 0; k < ncomp; ++k) {
            const double dxk_dlogn =
                state.composition[k] * ((k == variable ? 1.0 : 0.0) - state.composition[variable]);
            drho_dlogn += sensitivity.density_composition_derivative[k] * dxk_dlogn;
        }
        dlogrho_dlogn[variable] = drho_dlogn / density;
    }
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t variable = 0; variable < ncomp; ++variable) {
            jacobian[species * ncomp + variable] += dlogrho_dlogn[variable];
        }
    }
    return jacobian;
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

std::vector<double> reactive_phase_residual_jacobian_row_major(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
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
    std::vector<double> activity_jac1 = phase_ln_activity_log_amount_jacobian(mixture, t, p, phase1);
    std::vector<double> activity_jac2 = phase_ln_activity_log_amount_jacobian(mixture, t, p, phase2);
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
                mixture,
                t,
                p,
                phase1,
                reaction_standard_states[static_cast<std::size_t>(reaction)]
            );
            const std::vector<double> reaction_jac2 = phase_reaction_standard_state_log_amount_jacobian(
                mixture,
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
                mixture,
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
                mixture,
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

    for (std::size_t species = 0; species < ncomp; ++species) {
        jacobian[row * nvars + species] = charges[species] * phase1.amounts[species];
    }
    ++row;
    for (std::size_t species = 0; species < ncomp; ++species) {
        jacobian[row * nvars + ncomp + species] = charges[species] * phase2.amounts[species];
    }
    return jacobian;
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
    const std::vector<double>& charges = mixture->args().z;
    out.phase_charge_residuals = {
        composition_charge(phase1.amounts, charges),
        composition_charge(phase2.amounts, charges),
    };
    out.neutral_phase_equilibrium_residuals = neutral_phase_residuals(charges, ln_activity1, ln_activity2);
    out.ionic_equilibrium_residuals = ionic_phase_residuals(charges, ln_activity1, ln_activity2);
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
        mixture,
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
    out.diagnostics_string["coupling_level"] = "single_native_residual_state";
    out.diagnostics_string["phase_model"] = "two_liquid_phases";
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
