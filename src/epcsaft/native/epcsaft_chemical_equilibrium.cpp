#include "epcsaft_chemical_equilibrium.h"

#include <cppad/cppad.hpp>
#include <Eigen/Dense>

#include "epcsaft_electrolyte.h"
#include "ideal_speciation_problem.h"
#include "equilibrium_nlp/second_order.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

PhaseStateCompositionSensitivityResult phase_state_ln_fugacity_composition_sensitivity_cpp(
    double t,
    double p,
    std::vector<double> x,
    int phase,
    const add_args& cppargs
);

namespace {

constexpr int STANDARD_STATE_MOLE_FRACTION_ACTIVITY = 0;
constexpr int STANDARD_STATE_IDEAL_MOLE_FRACTION = 1;
constexpr int STANDARD_STATE_CONCENTRATION = 2;

int phase_token_to_int_chemical(const std::string& phase) {
    if (phase == "liq" || phase == "liquid" || phase == "aq" || phase == "org") {
        return 0;
    }
    if (phase == "vap" || phase == "vapor" || phase == "gas") {
        return 1;
    }
    throw ValueError("phase must be 'liq' or 'vap'.");
}

double max_abs_chemical(const std::vector<double>& values) {
    double out = 0.0;
    for (double value : values) {
        out = std::max(out, std::abs(value));
    }
    return out;
}

bool standard_states_all_ideal_mole_fraction(const std::vector<int>& standard_states) {
    for (int value : standard_states) {
        if (value != STANDARD_STATE_IDEAL_MOLE_FRACTION) {
            return false;
        }
    }
    return true;
}

bool standard_states_require_phase_state(const std::vector<int>& standard_states) {
    for (int value : standard_states) {
        if (value == STANDARD_STATE_MOLE_FRACTION_ACTIVITY || value == STANDARD_STATE_CONCENTRATION) {
            return true;
        }
    }
    return false;
}

bool standard_states_include_activity(const std::vector<int>& standard_states) {
    for (int value : standard_states) {
        if (value == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
            return true;
        }
    }
    return false;
}

std::string standard_state_label(int value) {
    if (value == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return "mole_fraction_activity";
    }
    if (value == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return "ideal_mole_fraction";
    }
    if (value == STANDARD_STATE_CONCENTRATION) {
        return "concentration";
    }
    throw ValueError("reaction standard state code is outside the native speciation contract.");
}

std::string standard_state_summary(const std::vector<int>& standard_states) {
    if (standard_states.empty()) {
        return "mole_fraction";
    }
    int first = standard_states.front();
    for (int value : standard_states) {
        if (value != first) {
            return "mixed_standard_state";
        }
    }
    if (first == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return "mole_fraction";
    }
    return standard_state_label(first);
}

std::vector<double> normalize_composition_chemical(const std::vector<double>& value, double min_mole_fraction) {
    std::vector<double> out(value.size(), min_mole_fraction);
    double total = 0.0;
    for (std::size_t i = 0; i < value.size(); ++i) {
        if (!std::isfinite(value[i]) || value[i] < -min_mole_fraction) {
            throw ValueError("initial_x values must be finite and non-negative.");
        }
        out[i] = std::max(value[i], min_mole_fraction);
        total += out[i];
    }
    if (!std::isfinite(total) || total <= 0.0) {
        throw ValueError("initial_x must have a positive finite sum.");
    }
    for (double& item : out) {
        item /= total;
    }
    return out;
}

std::vector<double> moles_from_log_amounts(const Eigen::VectorXd& log_n) {
    std::vector<double> out(static_cast<std::size_t>(log_n.size()), 0.0);
    for (Eigen::Index i = 0; i < log_n.size(); ++i) {
        out[static_cast<std::size_t>(i)] = std::exp(std::max(-700.0, std::min(700.0, log_n[i])));
    }
    return out;
}

std::vector<double> composition_from_moles(const std::vector<double>& n, double min_mole_fraction) {
    double total = std::accumulate(n.begin(), n.end(), 0.0);
    if (!std::isfinite(total) || total <= 0.0) {
        throw ValueError("chemical equilibrium iterate produced invalid mole amounts.");
    }
    std::vector<double> x(n.size(), min_mole_fraction);
    double clipped_total = 0.0;
    for (std::size_t i = 0; i < n.size(); ++i) {
        x[i] = std::max(n[i] / total, min_mole_fraction);
        clipped_total += x[i];
    }
    for (double& item : x) {
        item /= clipped_total;
    }
    return x;
}

struct ChemicalEvaluation {
    std::vector<double> n;
    std::vector<double> x;
    std::vector<double> gamma;
    bool has_phase_state = false;
    PhaseStateCompositionSensitivityResult phase_state;
    bool has_activity_sensitivity = false;
    std::vector<double> ln_activity_coefficient;
    std::vector<double> ln_activity_coefficient_jacobian_row_major;
    std::vector<double> mass_residuals;
    double charge_residual = 0.0;
    std::vector<double> reaction_residuals;
    std::vector<double> residuals;
    double residual_norm = std::numeric_limits<double>::infinity();
};

struct ChemicalEvaluationCounters {
    int residual_evaluations = 0;
    int jacobian_evaluations = 0;
    int activity_evaluations = 0;
};

struct ChemicalDerivativeSelection {
    std::string backend = "";
    std::string capability_path = "";
    bool derivative_available = false;
};

ChemicalDerivativeSelection select_chemical_derivative_backend(
    const ChemicalEquilibriumOptionsNative& options,
    const std::vector<int>& reaction_standard_states
) {
    ChemicalDerivativeSelection selection;
    const std::string requested = options.jacobian_backend;
    if (requested != "auto" && requested != "analytic") {
        if (requested != "cppad") {
            throw ValueError("chemical equilibrium jacobian_backend must be 'auto', 'analytic', or 'cppad'.");
        }
    }
    if (standard_states_all_ideal_mole_fraction(reaction_standard_states)) {
        selection.backend = requested == "cppad" ? "cppad" : "analytic";
        selection.capability_path = requested == "cppad"
            ? "chemical_equilibrium:ideal_mole_fraction:cppad_log_amounts"
            : "chemical_equilibrium:ideal_mole_fraction:log_amounts";
        selection.derivative_available = true;
        return selection;
    }
    selection.backend = "cppad_implicit";
    selection.capability_path =
        "chemical_equilibrium:" + standard_state_summary(reaction_standard_states) + ":phase_state_cppad_implicit";
    selection.derivative_available = true;
    return selection;
}

bool should_evaluate_activity_coefficients(const ChemicalEquilibriumOptionsNative& options) {
    const std::string mode = options.activity_output;
    if (mode == "always") {
        return true;
    }
    if (mode == "auto" || mode == "never") {
        return false;
    }
    throw ValueError("chemical equilibrium activity_output must be 'auto', 'always', or 'never'.");
}

Eigen::MatrixXd analytic_ideal_log_amount_jacobian(
    const ChemicalEvaluation& base,
    const Eigen::MatrixXd& balances,
    const Eigen::MatrixXd& reactions,
    const std::vector<double>& charges,
    double min_mole_fraction
) {
    const Eigen::Index nvars = static_cast<Eigen::Index>(base.n.size());
    const Eigen::Index rows = balances.rows() + 1 + reactions.rows();
    if (static_cast<Eigen::Index>(base.x.size()) != nvars) {
        throw ValueError("analytic chemical-equilibrium Jacobian requires matching composition and variable sizes.");
    }
    for (double value : base.x) {
        if (!(std::isfinite(value) && value > min_mole_fraction)) {
            throw ValueError("analytic chemical-equilibrium Jacobian requires an unclipped positive composition.");
        }
    }
    Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(rows, nvars);
    for (Eigen::Index r = 0; r < balances.rows(); ++r) {
        for (Eigen::Index j = 0; j < nvars; ++j) {
            jac(r, j) = balances(r, j) * base.n[static_cast<std::size_t>(j)];
        }
    }
    const Eigen::Index charge_row = balances.rows();
    if (static_cast<Eigen::Index>(charges.size()) == nvars) {
        for (Eigen::Index j = 0; j < nvars; ++j) {
            jac(charge_row, j) = charges[static_cast<std::size_t>(j)] * base.n[static_cast<std::size_t>(j)];
        }
    }
    for (Eigen::Index r = 0; r < reactions.rows(); ++r) {
        double stoich_sum = 0.0;
        for (Eigen::Index i = 0; i < reactions.cols(); ++i) {
            stoich_sum += reactions(r, i);
        }
        const Eigen::Index row = balances.rows() + 1 + r;
        for (Eigen::Index j = 0; j < nvars; ++j) {
            jac(row, j) = reactions(r, j) - stoich_sum * base.x[static_cast<std::size_t>(j)];
        }
    }
    return jac;
}

Eigen::MatrixXd cppad_ideal_log_amount_jacobian(
    const Eigen::VectorXd& log_n,
    const Eigen::MatrixXd& balances,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    const std::vector<double>& charges
) {
    using CppADScalar = CppAD::AD<double>;
    const Eigen::Index nvars = log_n.size();
    const Eigen::Index rows = balances.rows() + 1 + reactions.rows();
    std::vector<CppADScalar> variables(static_cast<std::size_t>(nvars));
    std::vector<double> point(static_cast<std::size_t>(nvars));
    for (Eigen::Index col = 0; col < nvars; ++col) {
        const double value = log_n[col];
        if (!std::isfinite(value)) {
            throw ValueError("CppAD chemical-equilibrium residual Jacobian requires finite log amounts.");
        }
        variables[static_cast<std::size_t>(col)] = value;
        point[static_cast<std::size_t>(col)] = value;
    }
    CppAD::Independent(variables);

    std::vector<CppADScalar> amounts(static_cast<std::size_t>(nvars));
    CppADScalar total = CppADScalar(0.0);
    for (Eigen::Index col = 0; col < nvars; ++col) {
        amounts[static_cast<std::size_t>(col)] = CppAD::exp(variables[static_cast<std::size_t>(col)]);
        total += amounts[static_cast<std::size_t>(col)];
    }
    std::vector<CppADScalar> composition(static_cast<std::size_t>(nvars));
    for (Eigen::Index col = 0; col < nvars; ++col) {
        composition[static_cast<std::size_t>(col)] = amounts[static_cast<std::size_t>(col)] / total;
    }

    std::vector<CppADScalar> outputs(static_cast<std::size_t>(rows), CppADScalar(0.0));
    for (Eigen::Index row = 0; row < balances.rows(); ++row) {
        CppADScalar residual = CppADScalar(0.0);
        for (Eigen::Index col = 0; col < nvars; ++col) {
            residual += balances(row, col) * amounts[static_cast<std::size_t>(col)];
        }
        outputs[static_cast<std::size_t>(row)] = residual;
    }

    const Eigen::Index charge_row = balances.rows();
    CppADScalar charge_residual = CppADScalar(0.0);
    if (charges.size() == static_cast<std::size_t>(nvars)) {
        for (Eigen::Index col = 0; col < nvars; ++col) {
            charge_residual += charges[static_cast<std::size_t>(col)] * amounts[static_cast<std::size_t>(col)];
        }
    }
    outputs[static_cast<std::size_t>(charge_row)] = charge_residual;

    for (Eigen::Index row = 0; row < reactions.rows(); ++row) {
        CppADScalar residual = -log_k[row];
        for (Eigen::Index col = 0; col < nvars; ++col) {
            residual += reactions(row, col) * CppAD::log(composition[static_cast<std::size_t>(col)]);
        }
        outputs[static_cast<std::size_t>(balances.rows() + 1 + row)] = residual;
    }

    CppAD::ADFun<double> function(variables, outputs);
    std::vector<double> jacobian = function.Jacobian(point);
    if (jacobian.size() != static_cast<std::size_t>(rows * nvars)) {
        throw ValueError("CppAD chemical-equilibrium residual Jacobian shape did not match the ideal route.");
    }
    Eigen::MatrixXd out(rows, nvars);
    for (Eigen::Index row = 0; row < rows; ++row) {
        for (Eigen::Index col = 0; col < nvars; ++col) {
            const double value = jacobian[static_cast<std::size_t>(row * nvars + col)];
            if (!std::isfinite(value)) {
                throw ValueError("CppAD chemical-equilibrium residual Jacobian produced a non-finite value.");
            }
            out(row, col) = value;
        }
    }
    return out;
}

std::vector<double> log_mole_fraction_terms(const ChemicalEvaluation& state, double floor) {
    std::vector<double> out(state.x.size(), 0.0);
    for (std::size_t index = 0; index < out.size(); ++index) {
        out[index] = std::log(std::max(state.x[index], floor));
    }
    return out;
}

std::vector<double> log_activity_terms(
    const ChemicalEvaluation& state,
    double floor
) {
    if (state.has_activity_sensitivity) {
        if (state.ln_activity_coefficient.size() != state.x.size()) {
            throw ValueError("chemical residual activity-coefficient payload length mismatch.");
        }
    } else if (state.phase_state.ln_fugacity.size() != state.x.size()) {
        throw ValueError("chemical residual phase-state activity payload length mismatch.");
    }
    std::vector<double> out(state.x.size(), 0.0);
    for (std::size_t index = 0; index < out.size(); ++index) {
        const double ln_gamma = state.has_activity_sensitivity
            ? state.ln_activity_coefficient[index]
            : state.phase_state.ln_fugacity[index];
        out[index] = std::log(std::max(state.x[index], floor)) + ln_gamma;
    }
    return out;
}

std::vector<double> log_concentration_terms(
    const ChemicalEvaluation& state,
    const PhaseStateCompositionSensitivityResult& phase_state,
    double floor
) {
    if (!(std::isfinite(phase_state.density) && phase_state.density > 0.0)) {
        throw ValueError("concentration chemical residual requires a finite positive molar density.");
    }
    std::vector<double> out(state.x.size(), 0.0);
    for (std::size_t index = 0; index < out.size(); ++index) {
        out[index] = std::log(std::max(state.x[index] * phase_state.density, floor));
    }
    return out;
}

std::vector<double> reaction_standard_state_log_terms(
    const ChemicalEvaluation& state,
    int standard_state,
    double floor
) {
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return log_mole_fraction_terms(state, floor);
    }
    if (!state.has_phase_state) {
        throw ValueError("chemical residual evaluation requires a phase-state derivative block.");
    }
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return log_activity_terms(state, floor);
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        return log_concentration_terms(state, state.phase_state, floor);
    }
    throw ValueError("reaction standard state code is outside the native speciation contract.");
}

std::vector<double> log_mole_fraction_log_amount_jacobian_row(
    const ChemicalEvaluation& state,
    std::size_t species
) {
    const std::size_t ncomp = state.x.size();
    std::vector<double> row(ncomp, 0.0);
    for (std::size_t variable = 0; variable < ncomp; ++variable) {
        row[variable] = (species == variable ? 1.0 : 0.0) - state.x[variable];
    }
    return row;
}

std::vector<double> log_activity_log_amount_jacobian_row(
    const ChemicalEvaluation& state,
    std::size_t species
) {
    const std::size_t ncomp = state.x.size();
    if (state.has_activity_sensitivity) {
        if (state.ln_activity_coefficient_jacobian_row_major.size() != ncomp * ncomp) {
            throw ValueError("chemical activity residual Jacobian requires activity-coefficient sensitivities.");
        }
    } else if (state.phase_state.jacobian_row_major.size() != ncomp * ncomp) {
        throw ValueError("chemical activity residual Jacobian requires supported phase-state sensitivities.");
    }
    std::vector<double> row = log_mole_fraction_log_amount_jacobian_row(state, species);
    for (std::size_t variable = 0; variable < ncomp; ++variable) {
        double dlnphi = 0.0;
        for (std::size_t k = 0; k < ncomp; ++k) {
            const double dxk_dlogn = state.x[k] * ((k == variable ? 1.0 : 0.0) - state.x[variable]);
            const double dlnactivity_dx = state.has_activity_sensitivity
                ? state.ln_activity_coefficient_jacobian_row_major[species * ncomp + k]
                : state.phase_state.jacobian_row_major[species * ncomp + k];
            dlnphi += dlnactivity_dx * dxk_dlogn;
        }
        row[variable] += dlnphi;
    }
    return row;
}

std::vector<double> log_concentration_log_amount_jacobian_row(
    const ChemicalEvaluation& state,
    const PhaseStateCompositionSensitivityResult& phase_state,
    std::size_t species
) {
    const std::size_t ncomp = state.x.size();
    if (phase_state.density_composition_derivative.size() != ncomp) {
        throw ValueError("concentration chemical residual Jacobian requires supported density composition sensitivities.");
    }
    if (!(std::isfinite(phase_state.density) && phase_state.density > 0.0)) {
        throw ValueError("concentration chemical residual Jacobian requires a finite positive molar density.");
    }
    std::vector<double> row = log_mole_fraction_log_amount_jacobian_row(state, species);
    for (std::size_t variable = 0; variable < ncomp; ++variable) {
        double drho_dlogn = 0.0;
        for (std::size_t k = 0; k < ncomp; ++k) {
            const double dxk_dlogn = state.x[k] * ((k == variable ? 1.0 : 0.0) - state.x[variable]);
            drho_dlogn += phase_state.density_composition_derivative[k] * dxk_dlogn;
        }
        row[variable] += drho_dlogn / phase_state.density;
    }
    return row;
}

std::vector<double> reaction_standard_state_log_amount_jacobian_row(
    const ChemicalEvaluation& state,
    int standard_state,
    std::size_t species
) {
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return log_mole_fraction_log_amount_jacobian_row(state, species);
    }
    if (!state.has_phase_state) {
        throw ValueError("chemical residual Jacobian requires a phase-state derivative block.");
    }
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return log_activity_log_amount_jacobian_row(state, species);
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        return log_concentration_log_amount_jacobian_row(state, state.phase_state, species);
    }
    throw ValueError("reaction standard state code is outside the native speciation contract.");
}

PhaseStateCompositionSensitivityResult evaluate_phase_state_sensitivity(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ChemicalEvaluation& state,
    int phase
) {
    PhaseStateCompositionSensitivityResult phase_state =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, state.x, phase, mixture->args());
    if (!phase_state.supported) {
        const std::string message = phase_state.message.empty()
            ? "phase-state fugacity composition sensitivity was not available."
            : phase_state.message;
        throw ValueError("chemical residual " + message);
    }
    if (phase_state.ln_fugacity.size() != state.x.size()) {
        throw ValueError("chemical residual phase-state fugacity payload length mismatch.");
    }
    return phase_state;
}

struct ActivityReferenceGroups {
    std::vector<int> cations;
    std::vector<int> anions;
    std::vector<int> solvents;
};

ActivityReferenceGroups activity_reference_groups(const add_args& args, std::size_t ncomp) {
    ActivityReferenceGroups out;
    if (args.z.size() != ncomp) {
        return out;
    }
    for (std::size_t i = 0; i < ncomp; ++i) {
        const double charge = args.z[i];
        if (charge > 1.0e-12) {
            out.cations.push_back(static_cast<int>(i));
        } else if (charge < -1.0e-12) {
            out.anions.push_back(static_cast<int>(i));
        } else {
            out.solvents.push_back(static_cast<int>(i));
        }
    }
    return out;
}

bool supports_component_activity_reference(const add_args& args, std::size_t ncomp) {
    ActivityReferenceGroups groups = activity_reference_groups(args, ncomp);
    return !groups.cations.empty() && !groups.anions.empty() && !groups.solvents.empty();
}

std::vector<double> activity_reference_composition(
    const std::vector<double>& x,
    const ActivityReferenceGroups& groups,
    double eps
) {
    const std::size_t ncomp = x.size();
    std::vector<double> out(ncomp, eps);
    double solvent_sum = 0.0;
    for (int idx : groups.solvents) {
        solvent_sum += x[static_cast<std::size_t>(idx)];
    }
    if (!(std::isfinite(solvent_sum) && solvent_sum > 0.0)) {
        throw ValueError("activity chemical residual requires a positive solvent fraction.");
    }
    const double solvent_budget = std::max(
        1.0 - eps * static_cast<double>(ncomp - groups.solvents.size()),
        eps * static_cast<double>(groups.solvents.size())
    );
    for (int idx : groups.solvents) {
        out[static_cast<std::size_t>(idx)] = solvent_budget * x[static_cast<std::size_t>(idx)] / solvent_sum;
    }
    double sum = std::accumulate(out.begin(), out.end(), 0.0);
    if (!(std::isfinite(sum) && sum > 0.0)) {
        throw ValueError("activity chemical residual produced an invalid reference composition.");
    }
    for (double& value : out) {
        value /= sum;
    }
    return out;
}

std::vector<double> activity_reference_composition_jacobian_row_major(
    const std::vector<double>& x,
    const ActivityReferenceGroups& groups,
    double eps
) {
    const std::size_t ncomp = x.size();
    std::vector<double> pre(ncomp, eps);
    std::vector<double> dpre(ncomp * ncomp, 0.0);
    double solvent_sum = 0.0;
    for (int idx : groups.solvents) {
        solvent_sum += x[static_cast<std::size_t>(idx)];
    }
    if (!(std::isfinite(solvent_sum) && solvent_sum > 0.0)) {
        throw ValueError("activity chemical residual requires a positive solvent fraction.");
    }
    const double solvent_budget = std::max(
        1.0 - eps * static_cast<double>(ncomp - groups.solvents.size()),
        eps * static_cast<double>(groups.solvents.size())
    );
    for (int row_idx : groups.solvents) {
        const std::size_t row = static_cast<std::size_t>(row_idx);
        pre[row] = solvent_budget * x[row] / solvent_sum;
        for (int col_idx : groups.solvents) {
            const std::size_t col = static_cast<std::size_t>(col_idx);
            const double numerator = (row == col ? solvent_sum : 0.0) - x[row];
            dpre[row * ncomp + col] = solvent_budget * numerator / (solvent_sum * solvent_sum);
        }
    }
    const double norm = std::accumulate(pre.begin(), pre.end(), 0.0);
    if (!(std::isfinite(norm) && norm > 0.0)) {
        throw ValueError("activity chemical residual produced an invalid reference composition.");
    }
    std::vector<double> dnorm(ncomp, 0.0);
    for (std::size_t col = 0; col < ncomp; ++col) {
        for (std::size_t row = 0; row < ncomp; ++row) {
            dnorm[col] += dpre[row * ncomp + col];
        }
    }
    std::vector<double> jac(ncomp * ncomp, 0.0);
    for (std::size_t row = 0; row < ncomp; ++row) {
        for (std::size_t col = 0; col < ncomp; ++col) {
            jac[row * ncomp + col] = (dpre[row * ncomp + col] * norm - pre[row] * dnorm[col]) / (norm * norm);
        }
    }
    return jac;
}

void assign_component_activity_sensitivity(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    ChemicalEvaluation& state,
    int phase
) {
    const std::size_t ncomp = state.x.size();
    const add_args& args = mixture->args();
    if (!supports_component_activity_reference(args, ncomp)) {
        return;
    }
    if (state.phase_state.jacobian_row_major.size() != ncomp * ncomp) {
        throw ValueError("activity chemical residual requires current phase-state fugacity sensitivities.");
    }
    const ActivityReferenceGroups groups = activity_reference_groups(args, ncomp);
    const double eps = 1.0e-12;
    std::vector<double> x_ref = activity_reference_composition(state.x, groups, eps);
    std::vector<double> dxref_dx = activity_reference_composition_jacobian_row_major(state.x, groups, eps);
    PhaseStateCompositionSensitivityResult ref_state =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, x_ref, phase, args);
    if (!ref_state.supported) {
        const std::string message = ref_state.message.empty()
            ? "activity reference-state fugacity sensitivity was not available."
            : ref_state.message;
        throw ValueError("chemical residual " + message);
    }
    if (ref_state.ln_fugacity.size() != ncomp || ref_state.jacobian_row_major.size() != ncomp * ncomp) {
        throw ValueError("activity chemical residual reference-state sensitivity shape mismatch.");
    }
    state.ln_activity_coefficient.assign(ncomp, 0.0);
    state.ln_activity_coefficient_jacobian_row_major.assign(ncomp * ncomp, 0.0);
    for (std::size_t species = 0; species < ncomp; ++species) {
        state.ln_activity_coefficient[species] =
            state.phase_state.ln_fugacity[species] - ref_state.ln_fugacity[species];
        for (std::size_t variable = 0; variable < ncomp; ++variable) {
            double ref_chain = 0.0;
            for (std::size_t ref_variable = 0; ref_variable < ncomp; ++ref_variable) {
                ref_chain += ref_state.jacobian_row_major[species * ncomp + ref_variable]
                    * dxref_dx[ref_variable * ncomp + variable];
            }
            state.ln_activity_coefficient_jacobian_row_major[species * ncomp + variable] =
                state.phase_state.jacobian_row_major[species * ncomp + variable] - ref_chain;
        }
    }
    state.has_activity_sensitivity = true;
}

std::vector<double> composition_log_amount_jacobian_row_major_chemical(const std::vector<double>& x) {
    const std::size_t ncomp = x.size();
    std::vector<double> jacobian(ncomp * ncomp, 0.0);
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t variable = 0; variable < ncomp; ++variable) {
            jacobian[species * ncomp + variable] =
                x[species] * ((species == variable ? 1.0 : 0.0) - x[variable]);
        }
    }
    return jacobian;
}

std::vector<double> composition_log_amount_hessian_tensor_chemical(const std::vector<double>& x) {
    const std::size_t ncomp = x.size();
    std::vector<double> hessian(ncomp * ncomp * ncomp, 0.0);
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t first = 0; first < ncomp; ++first) {
            const double first_delta = species == first ? 1.0 : 0.0;
            const double x_first = x[first];
            for (std::size_t second = 0; second < ncomp; ++second) {
                const double second_delta = species == second ? 1.0 : 0.0;
                const double cross_delta = first == second ? 1.0 : 0.0;
                const double x_second = x[second];
                hessian[species * ncomp * ncomp + first * ncomp + second] =
                    x[species] * (second_delta - x_second) * (first_delta - x_first)
                    - x[species] * x_first * (cross_delta - x_second);
            }
        }
    }
    return hessian;
}

std::vector<double> log_mole_fraction_log_amount_hessian_tensor_chemical(const std::vector<double>& x) {
    const std::size_t ncomp = x.size();
    std::vector<double> hessian(ncomp * ncomp * ncomp, 0.0);
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t first = 0; first < ncomp; ++first) {
            const double x_first = x[first];
            for (std::size_t second = 0; second < ncomp; ++second) {
                hessian[species * ncomp * ncomp + first * ncomp + second] =
                    -x_first * ((first == second ? 1.0 : 0.0) - x[second]);
            }
        }
    }
    return hessian;
}

std::vector<double> activity_reference_composition_hessian_tensor_row_major(
    const std::vector<double>& x,
    const ActivityReferenceGroups& groups,
    double eps
) {
    using CppADScalar = CppAD::AD<double>;
    const std::size_t ncomp = x.size();
    std::vector<CppADScalar> variables(ncomp);
    std::vector<double> point(ncomp, 0.0);
    for (std::size_t index = 0; index < ncomp; ++index) {
        variables[index] = x[index];
        point[index] = x[index];
    }
    CppAD::Independent(variables);
    std::vector<CppADScalar> pre(ncomp, CppADScalar(eps));
    CppADScalar solvent_sum = CppADScalar(0.0);
    for (int idx : groups.solvents) {
        solvent_sum += variables[static_cast<std::size_t>(idx)];
    }
    const double solvent_budget = std::max(
        1.0 - eps * static_cast<double>(ncomp - groups.solvents.size()),
        eps * static_cast<double>(groups.solvents.size())
    );
    for (int row_idx : groups.solvents) {
        const std::size_t row = static_cast<std::size_t>(row_idx);
        pre[row] = solvent_budget * variables[row] / solvent_sum;
    }
    CppADScalar norm = CppADScalar(0.0);
    for (const CppADScalar& value : pre) {
        norm += value;
    }
    std::vector<CppADScalar> outputs(ncomp, CppADScalar(0.0));
    for (std::size_t row = 0; row < ncomp; ++row) {
        outputs[row] = pre[row] / norm;
    }
    CppAD::ADFun<double> function(variables, outputs);
    std::vector<double> tensor(ncomp * ncomp * ncomp, 0.0);
    for (std::size_t output = 0; output < ncomp; ++output) {
        std::vector<double> weights(ncomp, 0.0);
        weights[output] = 1.0;
        std::vector<double> hessian = function.Hessian(point, weights);
        for (std::size_t first = 0; first < ncomp; ++first) {
            for (std::size_t second = 0; second < ncomp; ++second) {
                tensor[output * ncomp * ncomp + first * ncomp + second] =
                    hessian[first * ncomp + second];
            }
        }
    }
    return tensor;
}

std::vector<double> activity_coefficient_composition_hessian_tensor_chemical(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ChemicalEvaluation& state,
    int phase
) {
    const std::size_t ncomp = state.x.size();
    if (!state.has_phase_state || state.phase_state.hessian_tensor_row_major.size() != ncomp * ncomp * ncomp) {
        throw ValueError("chemical activity Hessian requires current phase-state Hessians.");
    }
    if (!state.has_activity_sensitivity) {
        return state.phase_state.hessian_tensor_row_major;
    }
    const add_args& args = mixture->args();
    const ActivityReferenceGroups groups = activity_reference_groups(args, ncomp);
    const double eps = 1.0e-12;
    const std::vector<double> x_ref = activity_reference_composition(state.x, groups, eps);
    const std::vector<double> dxref_dx = activity_reference_composition_jacobian_row_major(state.x, groups, eps);
    const std::vector<double> d2xref_dx2 =
        activity_reference_composition_hessian_tensor_row_major(state.x, groups, eps);
    PhaseStateCompositionSensitivityResult ref_state =
        phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, x_ref, phase, args);
    if (!ref_state.supported
        || ref_state.jacobian_row_major.size() != ncomp * ncomp
        || ref_state.hessian_tensor_row_major.size() != ncomp * ncomp * ncomp) {
        throw ValueError("chemical activity Hessian requires reference phase-state Hessians.");
    }
    std::vector<double> hessian = state.phase_state.hessian_tensor_row_major;
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t first = 0; first < ncomp; ++first) {
            for (std::size_t second = 0; second < ncomp; ++second) {
                double ref_chain = 0.0;
                for (std::size_t ref_first = 0; ref_first < ncomp; ++ref_first) {
                    ref_chain += ref_state.jacobian_row_major[species * ncomp + ref_first]
                        * d2xref_dx2[ref_first * ncomp * ncomp + first * ncomp + second];
                    for (std::size_t ref_second = 0; ref_second < ncomp; ++ref_second) {
                        ref_chain += ref_state.hessian_tensor_row_major[
                            species * ncomp * ncomp + ref_first * ncomp + ref_second
                        ] * dxref_dx[ref_first * ncomp + first] * dxref_dx[ref_second * ncomp + second];
                    }
                }
                hessian[species * ncomp * ncomp + first * ncomp + second] -= ref_chain;
            }
        }
    }
    return hessian;
}

std::vector<double> log_activity_log_amount_hessian_tensor_chemical(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ChemicalEvaluation& state,
    int phase
) {
    const std::size_t ncomp = state.x.size();
    const std::vector<double> composition_jacobian =
        composition_log_amount_jacobian_row_major_chemical(state.x);
    const std::vector<double> composition_hessian =
        composition_log_amount_hessian_tensor_chemical(state.x);
    const std::vector<double> activity_hessian =
        activity_coefficient_composition_hessian_tensor_chemical(mixture, t, p, state, phase);
    const std::vector<double>& activity_jacobian = state.has_activity_sensitivity
        ? state.ln_activity_coefficient_jacobian_row_major
        : state.phase_state.jacobian_row_major;
    if (activity_jacobian.size() != ncomp * ncomp || activity_hessian.size() != ncomp * ncomp * ncomp) {
        throw ValueError("chemical activity Hessian shape does not match species count.");
    }
    std::vector<double> hessian = log_mole_fraction_log_amount_hessian_tensor_chemical(state.x);
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t first = 0; first < ncomp; ++first) {
            for (std::size_t second = 0; second < ncomp; ++second) {
                double value = hessian[species * ncomp * ncomp + first * ncomp + second];
                for (std::size_t k = 0; k < ncomp; ++k) {
                    value += activity_jacobian[species * ncomp + k]
                        * composition_hessian[k * ncomp * ncomp + first * ncomp + second];
                    for (std::size_t l = 0; l < ncomp; ++l) {
                        value += activity_hessian[species * ncomp * ncomp + k * ncomp + l]
                            * composition_jacobian[k * ncomp + first]
                            * composition_jacobian[l * ncomp + second];
                    }
                }
                hessian[species * ncomp * ncomp + first * ncomp + second] = value;
            }
        }
    }
    return hessian;
}

std::vector<double> log_concentration_log_amount_hessian_tensor_chemical(const ChemicalEvaluation& state) {
    const std::size_t ncomp = state.x.size();
    if (!state.has_phase_state
        || state.phase_state.density_composition_derivative.size() != ncomp
        || state.phase_state.density_composition_hessian_row_major.size() != ncomp * ncomp
        || !(std::isfinite(state.phase_state.density) && state.phase_state.density > 0.0)) {
        throw ValueError("chemical concentration Hessian requires density composition Hessians.");
    }
    const std::vector<double> composition_jacobian =
        composition_log_amount_jacobian_row_major_chemical(state.x);
    const std::vector<double> composition_hessian =
        composition_log_amount_hessian_tensor_chemical(state.x);
    const double density = state.phase_state.density;
    std::vector<double> dlogrho(ncomp, 0.0);
    for (std::size_t variable = 0; variable < ncomp; ++variable) {
        for (std::size_t k = 0; k < ncomp; ++k) {
            dlogrho[variable] += state.phase_state.density_composition_derivative[k]
                * composition_jacobian[k * ncomp + variable];
        }
        dlogrho[variable] /= density;
    }
    std::vector<double> d2logrho(ncomp * ncomp, 0.0);
    for (std::size_t first = 0; first < ncomp; ++first) {
        for (std::size_t second = 0; second < ncomp; ++second) {
            double d2rho = 0.0;
            for (std::size_t k = 0; k < ncomp; ++k) {
                d2rho += state.phase_state.density_composition_derivative[k]
                    * composition_hessian[k * ncomp * ncomp + first * ncomp + second];
                for (std::size_t l = 0; l < ncomp; ++l) {
                    d2rho += state.phase_state.density_composition_hessian_row_major[k * ncomp + l]
                        * composition_jacobian[k * ncomp + first]
                        * composition_jacobian[l * ncomp + second];
                }
            }
            d2logrho[first * ncomp + second] =
                d2rho / density - dlogrho[first] * dlogrho[second];
        }
    }
    std::vector<double> hessian = log_mole_fraction_log_amount_hessian_tensor_chemical(state.x);
    for (std::size_t species = 0; species < ncomp; ++species) {
        for (std::size_t first = 0; first < ncomp; ++first) {
            for (std::size_t second = 0; second < ncomp; ++second) {
                hessian[species * ncomp * ncomp + first * ncomp + second] +=
                    d2logrho[first * ncomp + second];
            }
        }
    }
    return hessian;
}

std::vector<double> reaction_standard_state_log_amount_hessian_tensor_chemical(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ChemicalEvaluation& state,
    int standard_state,
    int phase
) {
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return log_mole_fraction_log_amount_hessian_tensor_chemical(state.x);
    }
    if (!state.has_phase_state) {
        throw ValueError("chemical residual Hessian requires a phase-state derivative block.");
    }
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return log_activity_log_amount_hessian_tensor_chemical(mixture, t, p, state, phase);
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        return log_concentration_log_amount_hessian_tensor_chemical(state);
    }
    throw ValueError("reaction standard state code is outside the native speciation contract.");
}

ChemicalEvaluation evaluate_chemical(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const Eigen::VectorXd& log_n,
    const Eigen::MatrixXd& balances,
    const Eigen::VectorXd& totals,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options,
    ChemicalEvaluationCounters* counters
) {
    ChemicalEvaluation out;
    if (counters != nullptr) {
        counters->residual_evaluations += 1;
    }
    out.n = moles_from_log_amounts(log_n);
    out.x = composition_from_moles(out.n, options.min_mole_fraction);
    const bool needs_phase_state = standard_states_require_phase_state(reaction_standard_states);
    const int phase = phase_token_to_int_chemical(options.phase);
    if (needs_phase_state) {
        if (counters != nullptr) {
            counters->activity_evaluations += 1;
        }
        out.phase_state = evaluate_phase_state_sensitivity(mixture, t, p, out, phase);
        out.has_phase_state = true;
        if (standard_states_include_activity(reaction_standard_states)) {
            assign_component_activity_sensitivity(mixture, t, p, out, phase);
        }
    }
    if (should_evaluate_activity_coefficients(options)) {
        if (counters != nullptr) {
            counters->activity_evaluations += 1;
        }
        if (out.has_activity_sensitivity) {
            out.gamma.reserve(out.ln_activity_coefficient.size());
            for (double value : out.ln_activity_coefficient) {
                out.gamma.push_back(std::exp(value));
            }
        } else if (needs_phase_state) {
            out.gamma.reserve(out.phase_state.ln_fugacity.size());
            for (double value : out.phase_state.ln_fugacity) {
                out.gamma.push_back(std::exp(value));
            }
        } else {
            out.gamma = std::vector<double>(out.x.size(), 1.0);
        }
    }

    Eigen::VectorXd n_vec = Eigen::Map<const Eigen::VectorXd>(out.n.data(), static_cast<Eigen::Index>(out.n.size()));
    Eigen::VectorXd mass = balances * n_vec - totals;
    out.mass_residuals.assign(mass.data(), mass.data() + mass.size());

    const std::vector<double>& charges = mixture->args().z;
    if (charges.size() == out.n.size()) {
        for (std::size_t i = 0; i < out.n.size(); ++i) {
            out.charge_residual += charges[i] * out.n[i];
        }
    }

    Eigen::VectorXd reaction(reactions.rows());
    for (Eigen::Index r = 0; r < reactions.rows(); ++r) {
        double value = -log_k[r];
        const int standard_state = reaction_standard_states[static_cast<std::size_t>(r)];
        const std::vector<double> log_terms = reaction_standard_state_log_terms(
            out,
            standard_state,
            options.min_mole_fraction
        );
        for (Eigen::Index i = 0; i < reactions.cols(); ++i) {
            value += reactions(r, i) * log_terms[static_cast<std::size_t>(i)];
        }
        reaction[r] = value;
    }
    out.reaction_residuals.assign(reaction.data(), reaction.data() + reaction.size());

    out.residuals.reserve(out.mass_residuals.size() + 1 + out.reaction_residuals.size());
    out.residuals.insert(out.residuals.end(), out.mass_residuals.begin(), out.mass_residuals.end());
    out.residuals.push_back(out.charge_residual);
    out.residuals.insert(out.residuals.end(), out.reaction_residuals.begin(), out.reaction_residuals.end());
    out.residual_norm = max_abs_chemical(out.residuals);
    return out;
}

Eigen::MatrixXd phase_state_log_amount_jacobian(
    const ChemicalEvaluation& current,
    const Eigen::MatrixXd& balances,
    const Eigen::MatrixXd& reactions,
    const std::vector<double>& charges,
    const std::vector<int>& reaction_standard_states
) {
    const Eigen::Index nvars = static_cast<Eigen::Index>(current.n.size());
    const Eigen::Index rows = balances.rows() + 1 + reactions.rows();
    Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(rows, nvars);
    for (Eigen::Index r = 0; r < balances.rows(); ++r) {
        for (Eigen::Index j = 0; j < nvars; ++j) {
            jac(r, j) = balances(r, j) * current.n[static_cast<std::size_t>(j)];
        }
    }
    const Eigen::Index charge_row = balances.rows();
    if (static_cast<Eigen::Index>(charges.size()) == nvars) {
        for (Eigen::Index j = 0; j < nvars; ++j) {
            jac(charge_row, j) = charges[static_cast<std::size_t>(j)] * current.n[static_cast<std::size_t>(j)];
        }
    }
    for (Eigen::Index r = 0; r < reactions.rows(); ++r) {
        const Eigen::Index row = balances.rows() + 1 + r;
        const int standard_state = reaction_standard_states[static_cast<std::size_t>(r)];
        for (Eigen::Index species = 0; species < reactions.cols(); ++species) {
            const double coefficient = reactions(r, species);
            if (coefficient == 0.0) {
                continue;
            }
            const std::vector<double> species_jacobian = reaction_standard_state_log_amount_jacobian_row(
                current,
                standard_state,
                static_cast<std::size_t>(species)
            );
            for (Eigen::Index variable = 0; variable < nvars; ++variable) {
                jac(row, variable) += coefficient * species_jacobian[static_cast<std::size_t>(variable)];
            }
        }
    }
    return jac;
}

std::vector<double> chemical_residual_hessian_tensor_row_major(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const ChemicalEvaluation& current,
    const Eigen::MatrixXd& balances,
    const Eigen::MatrixXd& reactions,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options
) {
    const std::size_t nvars = current.n.size();
    const std::size_t rows = static_cast<std::size_t>(balances.rows() + 1 + reactions.rows());
    std::vector<double> tensor(rows * nvars * nvars, 0.0);
    for (Eigen::Index row = 0; row < balances.rows(); ++row) {
        for (std::size_t species = 0; species < nvars; ++species) {
            tensor[static_cast<std::size_t>(row) * nvars * nvars + species * nvars + species] =
                balances(row, static_cast<Eigen::Index>(species)) * current.n[species];
        }
    }
    const Eigen::Index charge_row = balances.rows();
    const std::vector<double>& charges = mixture->args().z;
    if (charges.size() == nvars) {
        for (std::size_t species = 0; species < nvars; ++species) {
            tensor[static_cast<std::size_t>(charge_row) * nvars * nvars + species * nvars + species] =
                charges[species] * current.n[species];
        }
    }
    const int phase = phase_token_to_int_chemical(options.phase);
    for (Eigen::Index reaction = 0; reaction < reactions.rows(); ++reaction) {
        const std::size_t residual_row = static_cast<std::size_t>(balances.rows() + 1 + reaction);
        const int standard_state = reaction_standard_states[static_cast<std::size_t>(reaction)];
        const std::vector<double> standard_state_hessian =
            reaction_standard_state_log_amount_hessian_tensor_chemical(
                mixture,
                t,
                p,
                current,
                standard_state,
                phase
            );
        if (standard_state_hessian.size() != nvars * nvars * nvars) {
            throw ValueError("chemical standard-state Hessian shape does not match species count.");
        }
        for (std::size_t species = 0; species < nvars; ++species) {
            const double coefficient = reactions(reaction, static_cast<Eigen::Index>(species));
            if (coefficient == 0.0) {
                continue;
            }
            for (std::size_t first = 0; first < nvars; ++first) {
                for (std::size_t second = 0; second < nvars; ++second) {
                    tensor[residual_row * nvars * nvars + first * nvars + second] +=
                        coefficient
                        * standard_state_hessian[species * nvars * nvars + first * nvars + second];
                }
            }
        }
    }
    return tensor;
}

Eigen::MatrixXd chemical_residual_jacobian(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const Eigen::VectorXd& log_n,
    const ChemicalEvaluation& current,
    const Eigen::MatrixXd& balances,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options,
    ChemicalEvaluationCounters* counters,
    ChemicalDerivativeSelection* selection
) {
    ChemicalDerivativeSelection selected = select_chemical_derivative_backend(options, reaction_standard_states);
    if (selection != nullptr) {
        *selection = selected;
    }
    if (selected.backend == "analytic") {
        if (counters != nullptr) {
            counters->jacobian_evaluations += 1;
        }
        return analytic_ideal_log_amount_jacobian(
            current,
            balances,
            reactions,
            mixture->args().z,
            options.min_mole_fraction
        );
    }
    if (selected.backend == "cppad") {
        if (counters != nullptr) {
            counters->jacobian_evaluations += 1;
        }
        return cppad_ideal_log_amount_jacobian(
            log_n,
            balances,
            reactions,
            log_k,
            mixture->args().z
        );
    }
    if (selected.backend == "cppad_implicit") {
        if (counters != nullptr) {
            counters->jacobian_evaluations += 1;
        }
        if (!current.has_phase_state) {
            throw ValueError("chemical-equilibrium implicit Jacobian requires a retained phase-state sensitivity.");
        }
        return phase_state_log_amount_jacobian(
            current,
            balances,
            reactions,
            mixture->args().z,
            reaction_standard_states
        );
    }
    (void)mixture;
    (void)log_n;
    (void)log_k;
    throw ValueError("chemical-equilibrium residual jacobian has no registered analytical or CppAD backend.");
}

double positive_scale_from_totals_chemical(const std::vector<double>& totals) {
    double scale = 0.0;
    for (double value : totals) {
        if (!std::isfinite(value)) {
            throw ValueError("chemical equilibrium material totals must be finite.");
        }
        scale += std::abs(value);
    }
    return std::max(1.0, scale);
}

std::vector<double> canonical_initial_amounts_chemical(
    const std::vector<double>& initial_x,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    int species_count,
    const std::vector<double>& total_vector,
    double min_mole_fraction
) {
    const std::vector<double> x = normalize_composition_chemical(initial_x, min_mole_fraction);
    double numerator = 0.0;
    double denominator = 0.0;
    for (int row = 0; row < balance_rows; ++row) {
        double projected = 0.0;
        for (int col = 0; col < species_count; ++col) {
            projected += balance_matrix_row_major[
                static_cast<std::size_t>(row) * static_cast<std::size_t>(species_count)
                + static_cast<std::size_t>(col)
            ] * x[static_cast<std::size_t>(col)];
        }
        numerator += projected * total_vector[static_cast<std::size_t>(row)];
        denominator += projected * projected;
    }
    const double total_scale = positive_scale_from_totals_chemical(total_vector);
    double scale = total_scale;
    if (std::isfinite(numerator) && std::isfinite(denominator) && denominator > 0.0) {
        scale = numerator / denominator;
    }
    if (!(std::isfinite(scale) && scale > 0.0)) {
        scale = total_scale;
    }
    const double floor = min_mole_fraction * std::max(1.0, scale);
    std::vector<double> amounts(x.size(), floor);
    for (std::size_t index = 0; index < x.size(); ++index) {
        amounts[index] = std::max(floor, x[index] * scale);
    }
    return amounts;
}

void require_size_chemical(const std::vector<double>& values, std::size_t expected, const std::string& label) {
    if (values.size() != expected) {
        throw ValueError(label + " has an invalid size.");
    }
}

void validate_nonideal_speciation_request_chemical(
    const epcsaft::native::equilibrium_nlp::IdealSpeciationRequest& request
) {
    if (request.species_count <= 0) {
        throw ValueError("Nonideal Ipopt speciation requires at least one species.");
    }
    if (request.balance_rows <= 0) {
        throw ValueError("Nonideal Ipopt speciation requires at least one material balance.");
    }
    if (request.reaction_rows <= 0) {
        throw ValueError("Nonideal Ipopt speciation requires at least one reaction.");
    }
    if (!(std::isfinite(request.min_mole_fraction) && request.min_mole_fraction > 0.0)) {
        throw ValueError("Nonideal Ipopt speciation requires a positive min_mole_fraction.");
    }
    const std::size_t species = static_cast<std::size_t>(request.species_count);
    require_size_chemical(
        request.balance_matrix_row_major,
        static_cast<std::size_t>(request.balance_rows) * species,
        "Nonideal Ipopt speciation balance matrix"
    );
    require_size_chemical(
        request.total_vector,
        static_cast<std::size_t>(request.balance_rows),
        "Nonideal Ipopt speciation total vector"
    );
    require_size_chemical(
        request.reaction_stoichiometry_row_major,
        static_cast<std::size_t>(request.reaction_rows) * species,
        "Nonideal Ipopt speciation reaction matrix"
    );
    require_size_chemical(
        request.log_equilibrium_constants,
        static_cast<std::size_t>(request.reaction_rows),
        "Nonideal Ipopt speciation log equilibrium constants"
    );
    require_size_chemical(request.initial_x, species, "Nonideal Ipopt speciation initial composition");
    if (!request.charges.empty()) {
        require_size_chemical(request.charges, species, "Nonideal Ipopt speciation charges");
    }
    for (double value : request.balance_matrix_row_major) {
        if (!std::isfinite(value)) {
            throw ValueError("Nonideal Ipopt speciation balance matrix must contain finite values.");
        }
    }
    for (double value : request.total_vector) {
        if (!std::isfinite(value)) {
            throw ValueError("Nonideal Ipopt speciation total vector must contain finite values.");
        }
    }
    for (double value : request.reaction_stoichiometry_row_major) {
        if (!std::isfinite(value)) {
            throw ValueError("Nonideal Ipopt speciation reaction matrix must contain finite values.");
        }
    }
    for (double value : request.log_equilibrium_constants) {
        if (!std::isfinite(value)) {
            throw ValueError("Nonideal Ipopt speciation log equilibrium constants must be finite.");
        }
    }
    for (double value : request.charges) {
        if (!std::isfinite(value)) {
            throw ValueError("Nonideal Ipopt speciation charges must be finite.");
        }
    }
}

bool has_nonzero_charge_chemical(
    const epcsaft::native::equilibrium_nlp::IdealSpeciationRequest& request
) {
    for (double charge : request.charges) {
        if (std::abs(charge) > 1.0e-12) {
            return true;
        }
    }
    return false;
}

bool charge_constraint_increases_rank_chemical(
    const epcsaft::native::equilibrium_nlp::IdealSpeciationRequest& request
) {
    if (!has_nonzero_charge_chemical(request)) {
        return false;
    }
    Eigen::MatrixXd balances(request.balance_rows, request.species_count);
    for (int row = 0; row < request.balance_rows; ++row) {
        for (int col = 0; col < request.species_count; ++col) {
            balances(row, col) = request.balance_matrix_row_major[
                static_cast<std::size_t>(row) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ];
        }
    }
    Eigen::MatrixXd with_charge(request.balance_rows + 1, request.species_count);
    with_charge.topRows(request.balance_rows) = balances;
    for (int col = 0; col < request.species_count; ++col) {
        with_charge(request.balance_rows, col) = request.charges[static_cast<std::size_t>(col)];
    }
    const Eigen::FullPivLU<Eigen::MatrixXd> base_rank(balances);
    const Eigen::FullPivLU<Eigen::MatrixXd> charged_rank(with_charge);
    return charged_rank.rank() > base_rank.rank();
}

std::vector<double> matrix_to_row_major_chemical(const Eigen::MatrixXd& matrix) {
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(matrix.rows() * matrix.cols()));
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
            out.push_back(matrix(row, col));
        }
    }
    return out;
}

Eigen::MatrixXd matrix_from_row_major(const std::vector<double>& values, int rows, int cols, const std::string& name);
Eigen::VectorXd vector_from_values(const std::vector<double>& values);

Eigen::VectorXd log_amount_vector_from_amounts(const std::vector<double>& amounts) {
    Eigen::VectorXd log_n(static_cast<Eigen::Index>(amounts.size()));
    for (Eigen::Index index = 0; index < log_n.size(); ++index) {
        const double amount = amounts[static_cast<std::size_t>(index)];
        if (!(std::isfinite(amount) && amount > 0.0)) {
            throw ValueError("chemical equilibrium implicit diagnostics require positive amounts.");
        }
        log_n[index] = std::log(amount);
    }
    return log_n;
}

class NonidealSpeciationResidualProblem final : public epcsaft::native::equilibrium_nlp::NlpProblem {
public:
    NonidealSpeciationResidualProblem(
        std::shared_ptr<ePCSAFTMixtureNative> mixture,
        double t,
        double p,
        epcsaft::native::equilibrium_nlp::IdealSpeciationRequest request,
        std::vector<int> reaction_standard_states,
        ChemicalEquilibriumOptionsNative options
    )
        : mixture_(std::move(mixture)),
          t_(t),
          p_(p),
          request_(std::move(request)),
          reaction_standard_states_(std::move(reaction_standard_states)),
          options_(std::move(options)),
          balances_(matrix_from_row_major(
              request_.balance_matrix_row_major,
              request_.balance_rows,
              request_.species_count,
              "nonideal residual speciation balance matrix"
          )),
          totals_(vector_from_values(request_.total_vector)),
          reactions_(matrix_from_row_major(
              request_.reaction_stoichiometry_row_major,
              request_.reaction_rows,
              request_.species_count,
              "nonideal residual speciation reaction matrix"
          )),
          log_k_(vector_from_values(request_.log_equilibrium_constants)) {
        total_scale_ = positive_scale_from_totals_chemical(request_.total_vector);
        include_charge_constraint_ = charge_constraint_increases_rank_chemical(request_);
        const std::vector<double> initial_amounts = canonical_initial_amounts_chemical(
              request_.initial_x,
              request_.balance_matrix_row_major,
              request_.balance_rows,
              request_.species_count,
              request_.total_vector,
              options_.min_mole_fraction
        );
        initial_log_amounts_ = log_amount_vector_from_amounts(initial_amounts);
    }

    std::string name() const override {
        return "nonideal_homogeneous_reactive_speciation_residual";
    }

    int variable_count() const override {
        return request_.species_count;
    }

    int constraint_count() const override {
        return request_.balance_rows + (include_charge_constraint_ ? 1 : 0);
    }

    int jacobian_nonzero_count() const override {
        int balance_nonzeros = 0;
        for (double coefficient : request_.balance_matrix_row_major) {
            if (coefficient != 0.0) {
                ++balance_nonzeros;
            }
        }
        int charge_nonzeros = 0;
        if (include_charge_constraint_) {
            for (double charge : request_.charges) {
                if (charge != 0.0) {
                    ++charge_nonzeros;
                }
            }
        }
        return balance_nonzeros + charge_nonzeros;
    }

    epcsaft::native::equilibrium_nlp::NlpBounds bounds() const override {
        epcsaft::native::equilibrium_nlp::NlpBounds out;
        const double lower = options_.min_mole_fraction * total_scale_;
        const double upper = std::max(1.0, 10.0 * total_scale_);
        out.variable_lower.assign(static_cast<std::size_t>(request_.species_count), std::log(lower));
        out.variable_upper.assign(static_cast<std::size_t>(request_.species_count), std::log(upper));
        out.constraint_lower = request_.total_vector;
        out.constraint_upper = request_.total_vector;
        if (include_charge_constraint_) {
            out.constraint_lower.push_back(0.0);
            out.constraint_upper.push_back(0.0);
        }
        return out;
    }

    std::vector<double> initial_point() const override {
        return std::vector<double>(initial_log_amounts_.data(), initial_log_amounts_.data() + initial_log_amounts_.size());
    }

    double objective(const std::vector<double>& variables) const override {
        ChemicalEvaluationCounters counters;
        ChemicalEvaluation current = evaluate(variables, &counters);
        Eigen::Map<const Eigen::VectorXd> reaction(
            current.reaction_residuals.data(),
            static_cast<Eigen::Index>(current.reaction_residuals.size())
        );
        return 0.5 * reaction.squaredNorm();
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        ChemicalEvaluationCounters counters;
        ChemicalEvaluation current = evaluate(variables, &counters);
        ChemicalDerivativeSelection selection;
        const Eigen::VectorXd log_n = vector_from_values(variables);
        const Eigen::MatrixXd jac = chemical_residual_jacobian(
            mixture_,
            t_,
            p_,
            log_n,
            current,
            balances_,
            reactions_,
            log_k_,
            reaction_standard_states_,
            options_,
            &counters,
            &selection
        );
        Eigen::VectorXd reaction_gradient = Eigen::VectorXd::Zero(jac.cols());
        const Eigen::Index reaction_offset = balances_.rows() + 1;
        for (Eigen::Index reaction = 0; reaction < static_cast<Eigen::Index>(current.reaction_residuals.size()); ++reaction) {
            reaction_gradient += current.reaction_residuals[static_cast<std::size_t>(reaction)]
                * jac.row(reaction_offset + reaction).transpose();
        }
        const Eigen::VectorXd gradient = reaction_gradient;
        return std::vector<double>(gradient.data(), gradient.data() + gradient.size());
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        std::vector<double> out(static_cast<std::size_t>(constraint_count()), 0.0);
        const std::vector<double> amounts = moles_from_log_amounts(vector_from_values(variables));
        for (int row = 0; row < request_.balance_rows; ++row) {
            for (int col = 0; col < request_.species_count; ++col) {
                out[static_cast<std::size_t>(row)] += request_.balance_matrix_row_major[
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(request_.species_count)
                    + static_cast<std::size_t>(col)
                ] * amounts[static_cast<std::size_t>(col)];
            }
        }
        if (include_charge_constraint_) {
            for (int col = 0; col < request_.species_count; ++col) {
                out[static_cast<std::size_t>(request_.balance_rows)] +=
                    request_.charges[static_cast<std::size_t>(col)] * amounts[static_cast<std::size_t>(col)];
            }
        }
        return out;
    }

    epcsaft::native::equilibrium_nlp::NlpJacobianStructure jacobian_structure() const override {
        epcsaft::native::equilibrium_nlp::NlpJacobianStructure out;
        out.rows.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        out.cols.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < request_.balance_rows; ++row) {
            for (int col = 0; col < request_.species_count; ++col) {
                const double coefficient = request_.balance_matrix_row_major[
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(request_.species_count)
                    + static_cast<std::size_t>(col)
                ];
                if (coefficient == 0.0) {
                    continue;
                }
                out.rows.push_back(row);
                out.cols.push_back(col);
            }
        }
        if (include_charge_constraint_) {
            for (int col = 0; col < request_.species_count; ++col) {
                if (request_.charges[static_cast<std::size_t>(col)] == 0.0) {
                    continue;
                }
                out.rows.push_back(request_.balance_rows);
                out.cols.push_back(col);
            }
        }
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        const std::vector<double> amounts = moles_from_log_amounts(vector_from_values(variables));
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < request_.balance_rows; ++row) {
            for (int col = 0; col < request_.species_count; ++col) {
                const double coefficient = request_.balance_matrix_row_major[
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(request_.species_count)
                    + static_cast<std::size_t>(col)
                ];
                if (coefficient == 0.0) {
                    continue;
                }
                out.push_back(coefficient * amounts[static_cast<std::size_t>(col)]);
            }
        }
        if (include_charge_constraint_) {
            for (int col = 0; col < request_.species_count; ++col) {
                const double charge = request_.charges[static_cast<std::size_t>(col)];
                if (charge == 0.0) {
                    continue;
                }
                out.push_back(charge * amounts[static_cast<std::size_t>(col)]);
            }
        }
        return out;
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
        if (constraint_multipliers.size() != static_cast<std::size_t>(constraint_count())) {
            throw ValueError("nonideal speciation Hessian multiplier vector size does not match constraints.");
        }
        ChemicalEvaluationCounters counters;
        const Eigen::VectorXd log_n = vector_from_values(variables);
        ChemicalEvaluation current = evaluate(variables, &counters);
        ChemicalDerivativeSelection selection;
        const Eigen::MatrixXd jac = chemical_residual_jacobian(
            mixture_,
            t_,
            p_,
            log_n,
            current,
            balances_,
            reactions_,
            log_k_,
            reaction_standard_states_,
            options_,
            &counters,
            &selection
        );
        const std::vector<double> residual_hessian = chemical_residual_hessian_tensor_row_major(
            mixture_,
            t_,
            p_,
            current,
            balances_,
            reactions_,
            reaction_standard_states_,
            options_
        );
        const std::size_t n = static_cast<std::size_t>(variable_count());
        const std::size_t reaction_rows = static_cast<std::size_t>(reactions_.rows());
        const std::size_t reaction_offset = static_cast<std::size_t>(balances_.rows() + 1);

        epcsaft::native::equilibrium_nlp::ResidualSecondOrderData residuals;
        residuals.residual_count = static_cast<int>(reaction_rows);
        residuals.variable_count = variable_count();
        residuals.residuals = current.reaction_residuals;
        residuals.jacobian_row_major.assign(reaction_rows * n, 0.0);
        residuals.residual_hessian_tensor_row_major.assign(reaction_rows * n * n, 0.0);
        residuals.backend = "cppad_implicit_speciation_residual";
        for (std::size_t row = 0; row < reaction_rows; ++row) {
            for (std::size_t col = 0; col < n; ++col) {
                residuals.jacobian_row_major[row * n + col] =
                    jac(static_cast<Eigen::Index>(reaction_offset + row), static_cast<Eigen::Index>(col));
            }
            std::copy(
                residual_hessian.begin()
                    + static_cast<std::ptrdiff_t>((reaction_offset + row) * n * n),
                residual_hessian.begin()
                    + static_cast<std::ptrdiff_t>((reaction_offset + row + 1) * n * n),
                residuals.residual_hessian_tensor_row_major.begin()
                    + static_cast<std::ptrdiff_t>(row * n * n)
            );
        }
        epcsaft::native::equilibrium_nlp::ObjectiveSecondOrderData objective =
            epcsaft::native::equilibrium_nlp::residual_quadratic_objective_second_order(residuals);

        epcsaft::native::equilibrium_nlp::ConstraintSecondOrderData constraints_data;
        constraints_data.constraint_count = constraint_count();
        constraints_data.variable_count = variable_count();
        constraints_data.values = constraints(variables);
        constraints_data.hessian_tensor_row_major.assign(
            static_cast<std::size_t>(constraint_count()) * n * n,
            0.0
        );
        constraints_data.has_hessian.assign(static_cast<std::size_t>(constraint_count()), true);
        constraints_data.backend = "cppad_implicit_speciation_residual";
        for (std::size_t row = 0; row < static_cast<std::size_t>(request_.balance_rows); ++row) {
            std::copy(
                residual_hessian.begin() + static_cast<std::ptrdiff_t>(row * n * n),
                residual_hessian.begin() + static_cast<std::ptrdiff_t>((row + 1) * n * n),
                constraints_data.hessian_tensor_row_major.begin() + static_cast<std::ptrdiff_t>(row * n * n)
            );
        }
        if (include_charge_constraint_) {
            const std::size_t charge_source = static_cast<std::size_t>(balances_.rows());
            const std::size_t charge_target = static_cast<std::size_t>(request_.balance_rows);
            std::copy(
                residual_hessian.begin() + static_cast<std::ptrdiff_t>(charge_source * n * n),
                residual_hessian.begin() + static_cast<std::ptrdiff_t>((charge_source + 1) * n * n),
                constraints_data.hessian_tensor_row_major.begin()
                    + static_cast<std::ptrdiff_t>(charge_target * n * n)
            );
        }
        return epcsaft::native::equilibrium_nlp::LagrangianHessianAssembler(variable_count()).values(
            objective_factor,
            objective,
            constraints_data,
            constraint_multipliers
        );
    }

    std::string hessian_backend() const override {
        return "cppad_implicit_speciation_residual";
    }

    epcsaft::native::equilibrium_nlp::NlpScaling scaling() const override {
        epcsaft::native::equilibrium_nlp::NlpScaling out;
        out.objective = 1.0;
        out.variables.assign(static_cast<std::size_t>(request_.species_count), 1.0);
        out.constraints.reserve(static_cast<std::size_t>(constraint_count()));
        for (double total : request_.total_vector) {
            out.constraints.push_back(1.0 / std::max(1.0, std::abs(total)));
        }
        if (include_charge_constraint_) {
            out.constraints.push_back(1.0);
        }
        return out;
    }

    ChemicalEvaluation evaluate(const std::vector<double>& variables, ChemicalEvaluationCounters* counters) const {
        if (variables.size() != static_cast<std::size_t>(request_.species_count)) {
            throw ValueError("nonideal residual speciation variables length must match species count.");
        }
        const Eigen::VectorXd log_n = vector_from_values(variables);
        return evaluate_chemical(
            mixture_,
            t_,
            p_,
            log_n,
            balances_,
            totals_,
            reactions_,
            log_k_,
            reaction_standard_states_,
            options_,
            counters
        );
    }

private:
    std::shared_ptr<ePCSAFTMixtureNative> mixture_;
    double t_ = 0.0;
    double p_ = 0.0;
    epcsaft::native::equilibrium_nlp::IdealSpeciationRequest request_;
    std::vector<int> reaction_standard_states_;
    ChemicalEquilibriumOptionsNative options_;
    Eigen::MatrixXd balances_;
    Eigen::VectorXd totals_;
    Eigen::MatrixXd reactions_;
    Eigen::VectorXd log_k_;
    Eigen::VectorXd initial_log_amounts_;
    double total_scale_ = 1.0;
    bool include_charge_constraint_ = false;
};

void add_chemical_implicit_sensitivity_diagnostics(
    ChemicalEquilibriumResultNative& result,
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& amounts,
    const Eigen::MatrixXd& balances,
    const Eigen::VectorXd& totals,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options,
    const std::string& derivative_backend
) {
    ChemicalEvaluationCounters counters;
    const Eigen::VectorXd log_n = log_amount_vector_from_amounts(amounts);
    ChemicalEvaluation current = evaluate_chemical(
        mixture,
        t,
        p,
        log_n,
        balances,
        totals,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        &counters
    );
    ChemicalDerivativeSelection selection;
    Eigen::MatrixXd residual_state = chemical_residual_jacobian(
        mixture,
        t,
        p,
        log_n,
        current,
        balances,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        &counters,
        &selection
    );
    Eigen::MatrixXd residual_parameter =
        Eigen::MatrixXd::Zero(residual_state.rows(), static_cast<Eigen::Index>(log_k.size()));
    const Eigen::Index reaction_offset = balances.rows() + 1;
    for (Eigen::Index reaction = 0; reaction < log_k.size(); ++reaction) {
        residual_parameter(reaction_offset + reaction, reaction) = -1.0;
    }
    Eigen::MatrixXd sensitivity = residual_state.colPivHouseholderQr().solve(-residual_parameter);
    std::vector<double> residuals = current.residuals;

    result.diagnostics_string["implicit_sensitivity_backend"] =
        derivative_backend == "cppad_implicit" ? "cppad_implicit" : derivative_backend + "_implicit";
    result.diagnostics_string["reactive_speciation_sensitivity_parameter"] = "log_equilibrium_constants";
    result.diagnostics_int["reactive_speciation_residual_rows"] = static_cast<int>(residual_state.rows());
    result.diagnostics_int["reactive_speciation_state_size"] = static_cast<int>(residual_state.cols());
    result.diagnostics_int["reactive_speciation_parameter_size"] = static_cast<int>(log_k.size());
    result.diagnostics_vector["reactive_speciation_state"] =
        std::vector<double>(log_n.data(), log_n.data() + log_n.size());
    result.diagnostics_vector["reactive_speciation_residual"] = residuals;
    result.diagnostics_vector["reactive_speciation_residual_state_jacobian_row_major"] =
        matrix_to_row_major_chemical(residual_state);
    result.diagnostics_vector["reactive_speciation_residual_parameter_jacobian_row_major"] =
        matrix_to_row_major_chemical(residual_parameter);
    result.diagnostics_vector["reactive_speciation_log_amount_sensitivity_to_log_k_row_major"] =
        matrix_to_row_major_chemical(sensitivity);
}

Eigen::MatrixXd matrix_from_row_major(const std::vector<double>& values, int rows, int cols, const std::string& name) {
    if (rows < 0 || cols < 0 || values.size() != static_cast<std::size_t>(rows * cols)) {
        throw ValueError(name + " has an invalid row-major matrix size.");
    }
    Eigen::MatrixXd out(rows, cols);
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            out(r, c) = values[static_cast<std::size_t>(r * cols + c)];
        }
    }
    return out;
}

Eigen::VectorXd vector_from_values(const std::vector<double>& values) {
    Eigen::VectorXd out(static_cast<Eigen::Index>(values.size()));
    for (Eigen::Index i = 0; i < out.size(); ++i) {
        out[i] = values[static_cast<std::size_t>(i)];
    }
    return out;
}

ChemicalEquilibriumResultNative solve_nonideal_speciation_chemical_equilibrium_ipopt(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const epcsaft::native::equilibrium_nlp::IdealSpeciationRequest& request,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options
) {
    validate_nonideal_speciation_request_chemical(request);
    if (reaction_standard_states.size() != static_cast<std::size_t>(request.reaction_rows)) {
        throw ValueError("Nonideal Ipopt speciation requires one standard state per reaction.");
    }
    const int phase = phase_token_to_int_chemical(options.phase);
    (void)phase;
    NonidealSpeciationResidualProblem problem(mixture, t, p, request, reaction_standard_states, options);

    epcsaft::native::equilibrium_nlp::IpoptSolveOptions solve_options;
    solve_options.max_iterations = options.max_iterations;
    solve_options.tolerance = options.tolerance;
    solve_options.hessian_mode = options.hessian_mode;
    solve_options.limited_memory_hessian = options.hessian_mode != "exact";
    solve_options.iteration_history_limit = options.iteration_history_limit;
    solve_options.linear_solver = options.linear_solver;
    solve_options.acceptable_tolerance = options.acceptable_tolerance;
    solve_options.constraint_violation_tolerance = options.constraint_violation_tolerance;
    solve_options.dual_infeasibility_tolerance = options.dual_infeasibility_tolerance;
    solve_options.complementarity_tolerance = options.complementarity_tolerance;
    solve_options.initial_variables = options.initial_variables;
    solve_options.initial_bound_lower_multipliers = options.initial_bound_lower_multipliers;
    solve_options.initial_bound_upper_multipliers = options.initial_bound_upper_multipliers;
    solve_options.initial_constraint_multipliers = options.initial_constraint_multipliers;
    epcsaft::native::equilibrium_nlp::IpoptSolveResult ipopt_result =
        epcsaft::native::equilibrium_nlp::solve_ipopt_nlp(problem, solve_options);
    if (ipopt_result.variables.size() != static_cast<std::size_t>(request.species_count)) {
        throw SolutionError("Ipopt returned an invalid nonideal reactive speciation variable vector.");
    }

    const std::vector<double> amounts = moles_from_log_amounts(vector_from_values(ipopt_result.variables));
    const Eigen::MatrixXd balances = matrix_from_row_major(
        request.balance_matrix_row_major,
        request.balance_rows,
        request.species_count,
        "nonideal speciation balance matrix"
    );
    const Eigen::VectorXd totals = vector_from_values(request.total_vector);
    const Eigen::MatrixXd reactions = matrix_from_row_major(
        request.reaction_stoichiometry_row_major,
        request.reaction_rows,
        request.species_count,
        "nonideal speciation reaction matrix"
    );
    const Eigen::VectorXd log_k = vector_from_values(request.log_equilibrium_constants);
    const Eigen::VectorXd log_n = log_amount_vector_from_amounts(amounts);
    ChemicalEvaluationCounters counters;
    ChemicalEvaluation current = evaluate_chemical(
        mixture,
        t,
        p,
        log_n,
        balances,
        totals,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        &counters
    );

    ChemicalEquilibriumResultNative result;
    result.success = ipopt_result.accepted && current.residual_norm <= options.tolerance;
    if (result.success) {
        result.message = "converged";
    } else if (!ipopt_result.accepted) {
        result.message = "Ipopt did not accept the nonideal reactive speciation NLP solution.";
    } else {
        result.message = "Ipopt nonideal reactive speciation residual acceptance gate failed";
    }
    result.composition = current.x;
    result.activity_coefficients = current.gamma;
    result.mass_balance_residuals = current.mass_residuals;
    result.charge_residual = current.charge_residual;
    result.reaction_residuals = current.reaction_residuals;
    result.diagnostics_double = ipopt_result.diagnostics_double;
    result.diagnostics_int = ipopt_result.diagnostics_int;
    result.diagnostics_bool = ipopt_result.diagnostics_bool;
    result.diagnostics_string = ipopt_result.diagnostics_string;
    result.continuation_variables = ipopt_result.variables;
    result.continuation_bound_lower_multipliers = ipopt_result.bound_lower_multipliers;
    result.continuation_bound_upper_multipliers = ipopt_result.bound_upper_multipliers;
    result.continuation_constraint_multipliers = ipopt_result.constraint_multipliers;
    result.iteration_history = ipopt_result.iteration_history;
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_chemical_equilibrium_native";
    result.diagnostics_string["problem_class"] = "homogeneous_nonideal_residual_speciation";
    result.diagnostics_string["activity_model"] = "eos_phase_state";
    result.diagnostics_string["variable_model"] = "log_species_amounts";
    result.diagnostics_string["activity_output"] = options.activity_output;
    result.diagnostics_string["activity_basis"] = standard_state_summary(reaction_standard_states);
    result.diagnostics_string["phase"] = options.phase;
    result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
    result.diagnostics_string["jacobian_backend"] = "cppad_implicit";
    result.diagnostics_string["derivative_backend"] = "cppad_implicit";
    result.diagnostics_string["derivative_capability_path"] =
        "chemical_equilibrium:" + standard_state_summary(reaction_standard_states)
        + ":ipopt_log_amount_residual_phase_state_cppad_implicit";
    result.diagnostics_string["selected_solver_backend"] = "native_ipopt";
    result.diagnostics_string["solver_selection_reason"] = "explicit_request";
    result.diagnostics_string["ipopt_solver_status"] = ipopt_result.solver_status;
    result.diagnostics_string["ipopt_application_status"] = ipopt_result.application_status;
    result.diagnostics_bool["derivative_available"] = true;
    result.diagnostics_bool["jacobian_available"] = true;
    result.diagnostics_bool["activity_coefficients_evaluated"] = !result.activity_coefficients.empty();
    result.diagnostics_bool["charge_constraint_in_nlp"] = charge_constraint_increases_rank_chemical(request);
    result.diagnostics_bool["reaction_residual_in_objective"] = true;
    result.diagnostics_bool["ipopt_solver_ran"] = ipopt_result.solver_ran;
    result.diagnostics_bool["ipopt_accepted"] = ipopt_result.accepted;
    result.diagnostics_int["residual_evaluation_count"] = counters.residual_evaluations;
    result.diagnostics_int["jacobian_evaluation_count"] = counters.jacobian_evaluations;
    result.diagnostics_int["activity_evaluation_count"] = counters.activity_evaluations;
    result.diagnostics_double["residual_norm"] = current.residual_norm;
    result.diagnostics_double["tolerance"] = options.tolerance;
    result.diagnostics_double["objective"] = ipopt_result.objective;
    result.diagnostics_vector["history"] = {};
    result.diagnostics_vector["phase_handoff_composition"] = current.x;
    add_chemical_implicit_sensitivity_diagnostics(
        result,
        mixture,
        t,
        p,
        amounts,
        balances,
        totals,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        "cppad_implicit"
    );
    return result;
}

} // namespace

ChemicalResidualEvaluationNative evaluate_chemical_equilibrium_residual_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& initial_x,
    const std::vector<double>& variables,
    bool has_variables,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options
) {
    const int ncomp = static_cast<int>(mixture->ncomp());
    if (initial_x.size() != static_cast<std::size_t>(ncomp)) {
        throw ValueError("initial_x length must match mixture component count.");
    }
    if (has_variables && variables.size() != static_cast<std::size_t>(ncomp)) {
        throw ValueError("chemical residual variables length must match mixture component count.");
    }
    if (balance_rows <= 0) {
        throw ValueError("chemical residual evaluation requires at least one material balance.");
    }
    if (total_vector.size() != static_cast<std::size_t>(balance_rows)) {
        throw ValueError("total_vector length must match balance row count.");
    }
    if (log_equilibrium_constants.size() != static_cast<std::size_t>(reaction_rows)) {
        throw ValueError("log equilibrium constant length must match reaction row count.");
    }
    if (reaction_standard_states.size() != static_cast<std::size_t>(reaction_rows)) {
        throw ValueError("reaction standard state length must match reaction row count.");
    }
    if (options.min_mole_fraction <= 0.0) {
        throw ValueError("chemical residual evaluation options contain invalid numerical controls.");
    }
    for (int standard_state : reaction_standard_states) {
        standard_state_label(standard_state);
    }

    Eigen::MatrixXd balances = matrix_from_row_major(balance_matrix_row_major, balance_rows, ncomp, "balance_matrix");
    Eigen::VectorXd totals = vector_from_values(total_vector);
    Eigen::MatrixXd reactions = matrix_from_row_major(
        reaction_stoichiometry_row_major,
        reaction_rows,
        ncomp,
        "reaction_stoichiometry"
    );
    Eigen::VectorXd log_k = vector_from_values(log_equilibrium_constants);
    Eigen::VectorXd log_n(static_cast<Eigen::Index>(ncomp));
    if (has_variables) {
        log_n = vector_from_values(variables);
    } else {
        std::vector<double> initial = normalize_composition_chemical(initial_x, options.min_mole_fraction);
        for (int i = 0; i < ncomp; ++i) {
            log_n[i] = std::log(std::max(initial[static_cast<std::size_t>(i)], options.min_mole_fraction));
        }
    }
    for (Eigen::Index i = 0; i < log_n.size(); ++i) {
        if (!std::isfinite(log_n[i])) {
            throw ValueError("chemical residual variables must be finite.");
        }
    }

    (void)phase_token_to_int_chemical(options.phase);
    const std::string activity_model = standard_states_require_phase_state(reaction_standard_states)
        ? "eos_phase_state"
        : "ideal";
    ChemicalDerivativeSelection derivative_selection = select_chemical_derivative_backend(
        options,
        reaction_standard_states
    );
    ChemicalEvaluationCounters counters;
    ChemicalEvaluation current = evaluate_chemical(
        mixture,
        t,
        p,
        log_n,
        balances,
        totals,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        &counters
    );
    Eigen::MatrixXd jac = chemical_residual_jacobian(
        mixture,
        t,
        p,
        log_n,
        current,
        balances,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        &counters,
        &derivative_selection
    );
    Eigen::VectorXd residual = Eigen::Map<const Eigen::VectorXd>(
        current.residuals.data(),
        static_cast<Eigen::Index>(current.residuals.size())
    );
    Eigen::VectorXd gradient = jac.transpose() * residual;

    ChemicalResidualEvaluationNative out;
    out.variables.assign(log_n.data(), log_n.data() + log_n.size());
    const double lower = std::log(options.min_mole_fraction);
    out.lower_bounds.assign(static_cast<std::size_t>(ncomp), lower);
    out.upper_bounds.assign(static_cast<std::size_t>(ncomp), 50.0);
    out.residual = current.residuals;
    out.jacobian_rows = static_cast<int>(jac.rows());
    out.jacobian_cols = static_cast<int>(jac.cols());
    out.jacobian_row_major.reserve(static_cast<std::size_t>(jac.rows() * jac.cols()));
    for (Eigen::Index r = 0; r < jac.rows(); ++r) {
        for (Eigen::Index c = 0; c < jac.cols(); ++c) {
            out.jacobian_row_major.push_back(jac(r, c));
        }
    }
    out.gradient.assign(gradient.data(), gradient.data() + gradient.size());
    out.objective = 0.5 * residual.squaredNorm();
    out.composition = current.x;
    out.activity_coefficients = current.gamma;
    out.mass_balance_residuals = current.mass_residuals;
    out.charge_residual = current.charge_residual;
    out.reaction_residuals = current.reaction_residuals;
    out.diagnostics_string["solver_language"] = "c++";
    out.diagnostics_string["native_entrypoint"] = "_evaluate_chemical_equilibrium_residual_native";
    out.diagnostics_string["problem_class"] = "homogeneous_chemical_equilibrium";
    out.diagnostics_string["activity_model"] = activity_model;
    out.diagnostics_string["activity_output"] = options.activity_output;
    out.diagnostics_string["activity_basis"] = standard_state_summary(reaction_standard_states);
    out.diagnostics_string["phase"] = options.phase;
    out.diagnostics_string["jacobian_backend"] = derivative_selection.backend;
    out.diagnostics_string["derivative_backend"] = derivative_selection.backend;
    out.diagnostics_string["derivative_capability_path"] = derivative_selection.capability_path;
    out.diagnostics_bool["derivative_available"] = derivative_selection.derivative_available;
    out.diagnostics_bool["jacobian_available"] = derivative_selection.derivative_available;
    out.diagnostics_bool["activity_coefficients_evaluated"] = !current.gamma.empty();
    out.diagnostics_int["residual_evaluation_count"] = counters.residual_evaluations;
    out.diagnostics_int["jacobian_evaluation_count"] = counters.jacobian_evaluations;
    out.diagnostics_double["residual_norm"] = current.residual_norm;
    out.diagnostics_double["objective"] = out.objective;
    out.diagnostics_vector["phase_handoff_composition"] = current.x;
    return out;
}

ChemicalEquilibriumResultNative chemical_equilibrium_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& initial_x,
    const std::vector<double>& balance_matrix_row_major,
    int balance_rows,
    const std::vector<double>& total_vector,
    const std::vector<double>& reaction_stoichiometry_row_major,
    int reaction_rows,
    const std::vector<double>& log_equilibrium_constants,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options
) {
    const int ncomp = static_cast<int>(mixture->ncomp());
    if (initial_x.size() != static_cast<std::size_t>(ncomp)) {
        throw ValueError("initial_x length must match mixture component count.");
    }
    if (balance_rows <= 0) {
        throw ValueError("chemical equilibrium requires at least one material balance.");
    }
    if (total_vector.size() != static_cast<std::size_t>(balance_rows)) {
        throw ValueError("total_vector length must match balance row count.");
    }
    if (log_equilibrium_constants.size() != static_cast<std::size_t>(reaction_rows)) {
        throw ValueError("log equilibrium constant length must match reaction row count.");
    }
    if (reaction_standard_states.size() != static_cast<std::size_t>(reaction_rows)) {
        throw ValueError("reaction standard state length must match reaction row count.");
    }
    for (int standard_state : reaction_standard_states) {
        standard_state_label(standard_state);
    }
    if (options.max_iterations < 0 || options.tolerance <= 0.0 || options.min_mole_fraction <= 0.0) {
        throw ValueError("chemical equilibrium options contain invalid numerical controls.");
    }
    if (options.solver_backend != "auto" && options.solver_backend != "ipopt") {
        throw ValueError("chemical equilibrium solver_backend must be 'auto' or 'ipopt'.");
    }
    if (options.jacobian_backend != "auto"
        && options.jacobian_backend != "analytic"
        && options.jacobian_backend != "cppad") {
        throw ValueError("chemical equilibrium jacobian_backend must be 'auto', 'analytic', or 'cppad'.");
    }
    std::vector<double> initial = normalize_composition_chemical(initial_x, options.min_mole_fraction);
    if (options.solver_backend != "ipopt") {
        throw ValueError("chemical equilibrium solve requires solver_backend='ipopt'.");
    }
    epcsaft::native::equilibrium_nlp::IdealSpeciationRequest request;
    request.species_count = ncomp;
    request.balance_rows = balance_rows;
    request.balance_matrix_row_major = balance_matrix_row_major;
    request.total_vector = total_vector;
    request.reaction_rows = reaction_rows;
    request.reaction_stoichiometry_row_major = reaction_stoichiometry_row_major;
    request.log_equilibrium_constants = log_equilibrium_constants;
    request.initial_x = initial;
    request.charges = mixture->args().z;
    request.min_mole_fraction = options.min_mole_fraction;
    if (standard_states_all_ideal_mole_fraction(reaction_standard_states)) {
        return epcsaft::native::equilibrium_nlp::solve_ideal_speciation_chemical_equilibrium_ipopt(request, options);
    }
    return solve_nonideal_speciation_chemical_equilibrium_ipopt(
        mixture,
        t,
        p,
        request,
        reaction_standard_states,
        options
    );
}
