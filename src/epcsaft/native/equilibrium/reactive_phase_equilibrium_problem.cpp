#include "epcsaft_equilibrium.h"

#include "equilibrium_helpers.h"

#ifdef EPCSAFT_HAS_CERES
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
#endif

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
using epcsaft::native::equilibrium::split_distance_tolerance;

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

std::vector<double> reactive_phase_residual_jacobian_row_major(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
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
    std::size_t row = 0;

    for (int balance = 0; balance < balance_rows; ++balance) {
        for (std::size_t species = 0; species < ncomp; ++species) {
            const double coefficient = balance_matrix_row_major[static_cast<std::size_t>(balance) * ncomp + species];
            jacobian[row * nvars + species] = coefficient * phase1.amounts[species];
            jacobian[row * nvars + ncomp + species] = coefficient * phase2.amounts[species];
        }
        ++row;
    }

    for (int reaction = 0; reaction < reaction_rows; ++reaction) {
        for (std::size_t variable = 0; variable < ncomp; ++variable) {
            double value = 0.0;
            for (std::size_t species = 0; species < ncomp; ++species) {
                value += reaction_stoichiometry_row_major[static_cast<std::size_t>(reaction) * ncomp + species]
                    * activity_jac1[species * ncomp + variable];
            }
            jacobian[row * nvars + variable] = value;
        }
        ++row;
    }

    for (int reaction = 0; reaction < reaction_rows; ++reaction) {
        for (std::size_t variable = 0; variable < ncomp; ++variable) {
            double value = 0.0;
            for (std::size_t species = 0; species < ncomp; ++species) {
                value += reaction_stoichiometry_row_major[static_cast<std::size_t>(reaction) * ncomp + species]
                    * activity_jac2[species * ncomp + variable];
            }
            jacobian[row * nvars + ncomp + variable] = value;
        }
        ++row;
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

bool should_compute_jacobian(const EquilibriumOptionsNative& options) {
    return options.jacobian_backend != "auto";
}

#ifdef EPCSAFT_HAS_CERES
std::string ceres_termination_type_name_reactive(ceres::TerminationType type) {
    switch (type) {
        case ceres::CONVERGENCE:
            return "convergence";
        case ceres::NO_CONVERGENCE:
            return "no_convergence";
        case ceres::FAILURE:
            return "failure";
        case ceres::USER_SUCCESS:
            return "user_success";
        case ceres::USER_FAILURE:
            return "user_failure";
        default:
            return "unknown";
    }
}

class ReactivePhaseEquilibriumCostFunction final : public ceres::CostFunction {
public:
    ReactivePhaseEquilibriumCostFunction(
        std::shared_ptr<ePCSAFTMixtureNative> mixture,
        double t,
        double p,
        std::vector<double> feed,
        EquilibriumOptionsNative options,
        std::vector<double> balance_matrix,
        int balance_rows,
        std::vector<double> total_vector,
        std::vector<double> reaction_stoichiometry,
        int reaction_rows,
        std::vector<double> log_equilibrium_constants,
        std::vector<int> reaction_standard_states,
        int residual_size,
        int variable_count
    )
        : mixture_(std::move(mixture)),
          t_(t),
          p_(p),
          feed_(std::move(feed)),
          options_(std::move(options)),
          balance_matrix_(std::move(balance_matrix)),
          balance_rows_(balance_rows),
          total_vector_(std::move(total_vector)),
          reaction_stoichiometry_(std::move(reaction_stoichiometry)),
          reaction_rows_(reaction_rows),
          log_equilibrium_constants_(std::move(log_equilibrium_constants)),
          reaction_standard_states_(std::move(reaction_standard_states)) {
        set_num_residuals(residual_size);
        mutable_parameter_block_sizes()->push_back(variable_count);
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const override {
        std::vector<double> variables(parameters[0], parameters[0] + parameter_block_sizes()[0]);
        ReactivePhaseResidualEvaluationNative eval = evaluate_reactive_phase_equilibrium_residual_native(
            mixture_,
            t_,
            p_,
            feed_,
            options_,
            balance_matrix_,
            balance_rows_,
            total_vector_,
            reaction_stoichiometry_,
            reaction_rows_,
            log_equilibrium_constants_,
            reaction_standard_states_,
            variables,
            true
        );
        for (std::size_t i = 0; i < eval.residual.size(); ++i) {
            residuals[i] = eval.residual[i];
        }
        if (jacobians != nullptr && jacobians[0] != nullptr) {
            if (eval.jacobian_row_major.size() != eval.residual.size() * variables.size()) {
                return false;
            }
            std::copy(eval.jacobian_row_major.begin(), eval.jacobian_row_major.end(), jacobians[0]);
        }
        return true;
    }

private:
    std::shared_ptr<ePCSAFTMixtureNative> mixture_;
    double t_;
    double p_;
    std::vector<double> feed_;
    EquilibriumOptionsNative options_;
    std::vector<double> balance_matrix_;
    int balance_rows_;
    std::vector<double> total_vector_;
    std::vector<double> reaction_stoichiometry_;
    int reaction_rows_;
    std::vector<double> log_equilibrium_constants_;
    std::vector<int> reaction_standard_states_;
};
#endif

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
    out.jacobian_rows = static_cast<int>(out.residual.size());
    out.jacobian_cols = static_cast<int>(eval_variables.size());
    out.phase_distance = phase_distance(phase1.composition, phase2.composition);
    if (should_compute_jacobian(options)) {
        out.jacobian_row_major = reactive_phase_residual_jacobian_row_major(
            mixture,
            t,
            p,
            balance_matrix_row_major,
            balance_rows,
            reaction_stoichiometry_row_major,
            reaction_rows,
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
    } else {
        out.gradient.assign(eval_variables.size(), 0.0);
    }

    out.diagnostics_string["residual_surface"] = "native_reactive_phase_equilibrium_coupled_state";
    out.diagnostics_string["variable_model"] = out.variable_model;
    out.diagnostics_string["residual_blocks"] = "element_balance,reaction_equilibrium,neutral_phase_equilibrium,ionic_equilibrium,phase_charge";
    out.diagnostics_string["solver_backend"] = "residual_surface_only";
    out.diagnostics_string["solver_method"] = "not_started_until_ceres_slice";
    out.diagnostics_string["jacobian_backend"] = should_compute_jacobian(options) ? "cppad_implicit" : "not_available";
    out.diagnostics_string["derivative_backend"] = should_compute_jacobian(options) ? "cppad_implicit" : "not_available";
    out.diagnostics_string["coupling_level"] = "single_native_residual_state";
    out.diagnostics_string["phase_model"] = "two_liquid_phases";
    out.diagnostics_bool["jacobian_available"] = should_compute_jacobian(options);
    out.diagnostics_bool["derivative_available"] = should_compute_jacobian(options);
    out.diagnostics_bool["solved_state_sensitivity_available"] = should_compute_jacobian(options);
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

EquilibriumResultNative reactive_phase_equilibrium_native(
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
    const std::vector<double>& initial_phase1,
    const std::vector<double>& initial_phase2,
    double initial_phase_fraction_phase2,
    bool has_initial_phases
) {
#ifndef EPCSAFT_HAS_CERES
    (void)mixture;
    (void)t;
    (void)p;
    (void)raw_feed;
    (void)options;
    (void)balance_matrix_row_major;
    (void)balance_rows;
    (void)total_vector;
    (void)reaction_stoichiometry_row_major;
    (void)reaction_rows;
    (void)log_equilibrium_constants;
    (void)reaction_standard_states;
    (void)initial_phase1;
    (void)initial_phase2;
    (void)initial_phase_fraction_phase2;
    (void)has_initial_phases;
    throw ValueError("Ceres support is required for native reactive phase-equilibrium solve.");
#else
    EquilibriumOptionsNative solve_options = options;
    solve_options.jacobian_backend = "cppad";
    ReactivePhaseResidualEvaluationNative initial = evaluate_reactive_phase_equilibrium_residual_native(
        mixture,
        t,
        p,
        raw_feed,
        solve_options,
        balance_matrix_row_major,
        balance_rows,
        total_vector,
        reaction_stoichiometry_row_major,
        reaction_rows,
        log_equilibrium_constants,
        reaction_standard_states,
        {},
        false,
        initial_phase1,
        initial_phase2,
        initial_phase_fraction_phase2,
        has_initial_phases
    );
    std::vector<double> variables = initial.variables;
    ceres::Problem problem;
    auto* cost = new ReactivePhaseEquilibriumCostFunction(
        mixture,
        t,
        p,
        raw_feed,
        solve_options,
        balance_matrix_row_major,
        balance_rows,
        total_vector,
        reaction_stoichiometry_row_major,
        reaction_rows,
        log_equilibrium_constants,
        reaction_standard_states,
        static_cast<int>(initial.residual.size()),
        static_cast<int>(variables.size())
    );
    problem.AddResidualBlock(cost, nullptr, variables.data());
    ceres::Solver::Options ceres_options;
    ceres_options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    ceres_options.linear_solver_type = ceres::DENSE_QR;
    ceres_options.max_num_iterations = options.max_iterations;
    ceres_options.minimizer_progress_to_stdout = false;
    ceres_options.logging_type = ceres::SILENT;
    ceres_options.function_tolerance = std::min(1.0e-12, std::max(1.0e-18, options.tolerance * options.tolerance));
    ceres_options.gradient_tolerance = std::min(1.0e-10, std::max(1.0e-14, options.tolerance));
    ceres_options.parameter_tolerance = std::min(1.0e-10, std::max(1.0e-14, options.tolerance));
    ceres::Solver::Summary summary;
    ceres::Solve(ceres_options, &problem, &summary);
    ReactivePhaseResidualEvaluationNative final_eval = evaluate_reactive_phase_equilibrium_residual_native(
        mixture,
        t,
        p,
        raw_feed,
        solve_options,
        balance_matrix_row_major,
        balance_rows,
        total_vector,
        reaction_stoichiometry_row_major,
        reaction_rows,
        log_equilibrium_constants,
        reaction_standard_states,
        variables,
        true
    );

    EquilibriumResultNative result;
    result.backend = "reactive_phase_equilibrium";
    result.problem_kind = "reactive_phase_equilibrium";
    result.stable = false;
    result.split_detected = final_eval.phase_distance > split_distance_tolerance(options);
    EquilibriumPhaseNative phase1;
    phase1.label = "liq1";
    phase1.composition = final_eval.phase1_composition;
    phase1.density = final_eval.phase1_density;
    phase1.temperature = t;
    phase1.pressure = p;
    phase1.phase_fraction = 1.0 - final_eval.phase_fraction_phase2;
    phase1.ln_fugacity_coefficient = final_eval.phase1_ln_fugacity_coefficient;
    EquilibriumPhaseNative phase2;
    phase2.label = "liq2";
    phase2.composition = final_eval.phase2_composition;
    phase2.density = final_eval.phase2_density;
    phase2.temperature = t;
    phase2.pressure = p;
    phase2.phase_fraction = final_eval.phase_fraction_phase2;
    phase2.ln_fugacity_coefficient = final_eval.phase2_ln_fugacity_coefficient;
    result.phases = {phase1, phase2};

    result.diagnostics_string = final_eval.diagnostics_string;
    result.diagnostics_string["native_entrypoint"] = "_solve_reactive_phase_equilibrium_native";
    result.diagnostics_string["equilibrium_route"] = "reactive_phase_equilibrium";
    result.diagnostics_string["solver_backend"] = "ceres";
    result.diagnostics_string["selected_solver_backend"] = "ceres";
    result.diagnostics_string["solver_method"] = "ceres_trust_region_coupled_reactive_phase_equilibrium";
    result.diagnostics_string["ceres_trust_region_strategy"] = "levenberg_marquardt";
    result.diagnostics_string["ceres_linear_solver"] = "dense_qr";
    result.diagnostics_string["ceres_termination_type"] = ceres_termination_type_name_reactive(summary.termination_type);
    result.diagnostics_string["ceres_summary"] = summary.BriefReport();
    result.diagnostics_string["jacobian_backend"] = "cppad_implicit";
    result.diagnostics_string["derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["solved_state_sensitivity_backend"] = "cppad_implicit";
    result.diagnostics_bool = final_eval.diagnostics_bool;
    result.diagnostics_bool["ceres_accepted_solve"] = true;
    result.diagnostics_bool["jacobian_available"] = true;
    result.diagnostics_bool["derivative_available"] = true;
    result.diagnostics_bool["solved_state_sensitivity_available"] = true;
    result.diagnostics_int = final_eval.diagnostics_int;
    result.diagnostics_int["ceres_iteration_count"] = static_cast<int>(summary.iterations.size());
    result.diagnostics_int["ceres_status"] = static_cast<int>(summary.termination_type);
    result.diagnostics_double = final_eval.diagnostics_double;
    result.diagnostics_double["ceres_initial_cost"] = summary.initial_cost;
    result.diagnostics_double["ceres_final_cost"] = summary.final_cost;
    result.diagnostics_double["phase_distance"] = final_eval.phase_distance;
    result.diagnostics_vector = final_eval.diagnostics_vector;
    result.diagnostics_vector["variables"] = final_eval.variables;
    result.diagnostics_vector["residual"] = final_eval.residual;
    result.diagnostics_vector["jacobian_row_major"] = final_eval.jacobian_row_major;
    return result;
#endif
}
