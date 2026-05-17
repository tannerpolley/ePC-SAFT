#include "epcsaft_chemical_equilibrium.h"

#include <cppad/cppad.hpp>
#include <Eigen/Dense>

#include "epcsaft_electrolyte.h"
#include "ideal_speciation_problem.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>

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
    const PhaseStateCompositionSensitivityResult& phase_state,
    double floor
) {
    if (phase_state.ln_fugacity.size() != state.x.size()) {
        throw ValueError("chemical residual phase-state fugacity payload length mismatch.");
    }
    std::vector<double> out(state.x.size(), 0.0);
    for (std::size_t index = 0; index < out.size(); ++index) {
        out[index] = std::log(std::max(state.x[index], floor)) + phase_state.ln_fugacity[index];
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
    const PhaseStateCompositionSensitivityResult* phase_state,
    int standard_state,
    double floor
) {
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return log_mole_fraction_terms(state, floor);
    }
    if (phase_state == nullptr) {
        throw ValueError("chemical residual evaluation requires a phase-state derivative block.");
    }
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return log_activity_terms(state, *phase_state, floor);
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        return log_concentration_terms(state, *phase_state, floor);
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
    const PhaseStateCompositionSensitivityResult& phase_state,
    std::size_t species
) {
    const std::size_t ncomp = state.x.size();
    if (phase_state.jacobian_row_major.size() != ncomp * ncomp) {
        throw ValueError("chemical activity residual Jacobian requires supported phase-state fugacity sensitivities.");
    }
    std::vector<double> row = log_mole_fraction_log_amount_jacobian_row(state, species);
    for (std::size_t variable = 0; variable < ncomp; ++variable) {
        double dlnphi = 0.0;
        for (std::size_t k = 0; k < ncomp; ++k) {
            const double dxk_dlogn = state.x[k] * ((k == variable ? 1.0 : 0.0) - state.x[variable]);
            dlnphi += phase_state.jacobian_row_major[species * ncomp + k] * dxk_dlogn;
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
    const PhaseStateCompositionSensitivityResult* phase_state,
    int standard_state,
    std::size_t species
) {
    if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
        return log_mole_fraction_log_amount_jacobian_row(state, species);
    }
    if (phase_state == nullptr) {
        throw ValueError("chemical residual Jacobian requires a phase-state derivative block.");
    }
    if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return log_activity_log_amount_jacobian_row(state, *phase_state, species);
    }
    if (standard_state == STANDARD_STATE_CONCENTRATION) {
        return log_concentration_log_amount_jacobian_row(state, *phase_state, species);
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
    }
    if (should_evaluate_activity_coefficients(options)) {
        if (counters != nullptr) {
            counters->activity_evaluations += 1;
        }
        if (needs_phase_state) {
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
            needs_phase_state ? &out.phase_state : nullptr,
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
    const std::vector<int>& reaction_standard_states,
    const PhaseStateCompositionSensitivityResult& phase_state
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
                &phase_state,
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
            if (!current.has_phase_state) {
                counters->activity_evaluations += 1;
            }
        }
        const PhaseStateCompositionSensitivityResult phase_state = current.has_phase_state
            ? current.phase_state
            : evaluate_phase_state_sensitivity(
                mixture,
                t,
                p,
                current,
                phase_token_to_int_chemical(options.phase)
            );
        return phase_state_log_amount_jacobian(
            current,
            balances,
            reactions,
            mixture->args().z,
            reaction_standard_states,
            phase_state
        );
    }
    (void)mixture;
    (void)log_n;
    (void)log_k;
    throw ValueError("chemical-equilibrium residual jacobian has no registered analytical or CppAD backend.");
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
    (void)t;
    (void)p;
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
    if (!standard_states_all_ideal_mole_fraction(reaction_standard_states)) {
        throw ValueError(
            "Native Ipopt nonideal reactive speciation requires a native Ipopt Gibbs/activity NLP route builder."
        );
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
    return epcsaft::native::equilibrium_nlp::solve_ideal_speciation_chemical_equilibrium_ipopt(request, options);
}
