#include "epcsaft_chemical_equilibrium.h"

#include <Eigen/Dense>

#include "epcsaft_electrolyte.h"
#include "ideal_speciation_problem.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>

namespace {

constexpr int STANDARD_STATE_MOLE_FRACTION_ACTIVITY = 0;
constexpr int STANDARD_STATE_IDEAL_MOLE_FRACTION = 1;
constexpr int STANDARD_STATE_CONCENTRATION = 2;

bool has_ionic_species(const std::shared_ptr<ePCSAFTMixtureNative>& mixture);

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

bool standard_states_need_concentration(const std::vector<int>& standard_states) {
    for (int value : standard_states) {
        if (value == STANDARD_STATE_CONCENTRATION) {
            return true;
        }
    }
    return false;
}

bool standard_states_need_activity(const std::vector<int>& standard_states) {
    for (int value : standard_states) {
        if (value == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
            return true;
        }
    }
    return false;
}

bool standard_states_all_ideal_mole_fraction(const std::vector<int>& standard_states) {
    for (int value : standard_states) {
        if (value != STANDARD_STATE_IDEAL_MOLE_FRACTION) {
            return false;
        }
    }
    return true;
}

std::string activity_model_for_standard_states(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<int>& standard_states
) {
    if (standard_states_need_activity(standard_states)) {
        return has_ionic_species(mixture) ? "epcsaft_component_activity" : "epcsaft_neutral_fugacity_activity";
    }
    if (standard_states_need_concentration(standard_states)) {
        return "concentration";
    }
    return "ideal";
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
    throw ValueError("reaction standard state contains an unsupported code.");
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

bool has_ionic_species(const std::shared_ptr<ePCSAFTMixtureNative>& mixture) {
    const std::vector<double>& charges = mixture->args().z;
    for (double charge : charges) {
        if (std::abs(charge) > 1.0e-12) {
            return true;
        }
    }
    return false;
}

std::vector<double> component_rich_reference_composition(
    std::size_t ncomp,
    std::size_t rich_index,
    double min_mole_fraction
) {
    if (rich_index >= ncomp) {
        throw ValueError("component-rich activity reference index is out of range.");
    }
    const double floor = std::max(min_mole_fraction, 1.0e-14);
    std::vector<double> out(ncomp, floor);
    out[rich_index] = std::max(floor, 1.0 - floor * static_cast<double>(ncomp - 1));
    double total = std::accumulate(out.begin(), out.end(), 0.0);
    if (!std::isfinite(total) || total <= 0.0) {
        throw ValueError("component-rich activity reference composition is invalid.");
    }
    for (double& item : out) {
        item /= total;
    }
    return out;
}

std::vector<double> neutral_fugacity_activity_coefficients(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& x,
    int phase_int,
    double min_mole_fraction
) {
    double rho = mixture->solve_density_scoped(t, p, x, phase_int, "chemical_equilibrium");
    std::shared_ptr<ePCSAFTStateNative> state = mixture->state(t, x, phase_int, false, 0.0, true, rho);
    std::vector<double> ln_phi = state->ln_fugacity_coefficient();
    if (ln_phi.size() != x.size()) {
        throw ValueError("native fugacity coefficient payload length does not match composition.");
    }

    std::vector<double> gamma(x.size(), 1.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        std::vector<double> x_ref = component_rich_reference_composition(x.size(), i, min_mole_fraction);
        double rho_ref = mixture->solve_density_scoped(
            t,
            p,
            x_ref,
            phase_int,
            "chemical_equilibrium_reference"
        );
        std::shared_ptr<ePCSAFTStateNative> ref = mixture->state(t, x_ref, phase_int, false, 0.0, true, rho_ref);
        std::vector<double> ln_phi_ref = ref->ln_fugacity_coefficient();
        if (ln_phi_ref.size() != x.size()) {
            throw ValueError("native fugacity reference payload length does not match composition.");
        }
        const double ln_gamma = std::max(-700.0, std::min(700.0, ln_phi[i] - ln_phi_ref[i]));
        gamma[i] = std::exp(ln_gamma);
    }
    return gamma;
}

std::vector<double> activity_coefficients(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const std::vector<double>& x,
    int phase_int,
    const std::string& activity_model,
    double min_mole_fraction,
    int* state_failure_count
) {
    if (activity_model == "ideal" || activity_model == "concentration") {
        return std::vector<double>(x.size(), 1.0);
    }
    try {
        if (activity_model == "epcsaft_neutral_fugacity_activity") {
            return neutral_fugacity_activity_coefficients(mixture, t, p, x, phase_int, min_mole_fraction);
        }
        if (activity_model != "epcsaft_component_activity") {
            throw ValueError("unknown chemical equilibrium activity model.");
        }
        double rho = mixture->solve_density_scoped(t, p, x, phase_int, "chemical_equilibrium");
        std::shared_ptr<ePCSAFTStateNative> state = mixture->state(t, x, phase_int, false, 0.0, true, rho);
        ActivityCoefficientNative activity = state->activity_coefficient_native(false, false, -1);
        if (activity.component_activity_coefficients.size() != x.size()) {
            throw ValueError("native activity coefficient payload length does not match composition.");
        }
        return activity.component_activity_coefficients;
    } catch (...) {
        if (state_failure_count != nullptr) {
            *state_failure_count += 1;
        }
        throw;
    }
}

struct ChemicalEvaluation {
    std::vector<double> n;
    std::vector<double> x;
    std::vector<double> gamma;
    std::vector<double> mass_residuals;
    double charge_residual = 0.0;
    std::vector<double> reaction_residuals;
    std::vector<double> residuals;
    double residual_norm = std::numeric_limits<double>::infinity();
};

struct ChemicalEvaluationCounters {
    int residual_evaluations = 0;
    int jacobian_evaluations = 0;
    int state_evaluations = 0;
    int activity_evaluations = 0;
    int density_solves = 0;
};

struct ChemicalDerivativeSelection {
    std::string backend = "";
    std::string capability_path = "";
    std::string not_available_reason = "";
    bool derivative_available = false;
};

ChemicalDerivativeSelection select_chemical_derivative_backend(
    const ChemicalEquilibriumOptionsNative& options,
    const std::vector<int>& reaction_standard_states
) {
    ChemicalDerivativeSelection selection;
    const std::string requested = options.jacobian_backend;
    if (requested == "cppad") {
        throw ValueError("CppAD chemical-equilibrium residual jacobian requires implemented log-species amount coverage.");
    }
    if (requested != "auto" && requested != "analytic") {
        throw ValueError("chemical equilibrium jacobian_backend must be 'auto', 'analytic', or 'cppad'.");
    }
    if (standard_states_all_ideal_mole_fraction(reaction_standard_states)) {
        selection.backend = "analytic";
        selection.capability_path = "chemical_equilibrium:ideal_mole_fraction:log_amounts";
        selection.derivative_available = true;
        return selection;
    }
    selection.backend = "analytic";
    selection.capability_path = "chemical_equilibrium:activity_or_concentration:log_amounts";
    selection.derivative_available = true;
    return selection;
}

bool should_evaluate_activity_coefficients(
    const std::vector<int>& standard_states,
    const ChemicalEquilibriumOptionsNative& options
) {
    const std::string mode = options.activity_output;
    if (mode == "always") {
        return true;
    }
    if (mode == "auto" || mode == "never") {
        return standard_states_need_activity(standard_states);
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

struct ChemicalSoftStartResult {
    bool enabled = true;
    bool attempted = false;
    bool success = false;
    bool used = false;
    int iterations = 0;
    double residual_norm = std::numeric_limits<double>::infinity();
    double initial_residual_norm = std::numeric_limits<double>::infinity();
    std::string rejection_reason = "not_attempted";
    std::vector<double> composition;
    std::vector<double> history;
    Eigen::VectorXd log_n;
};

struct IdealGibbsEvaluation {
    bool feasible = false;
    double objective = std::numeric_limits<double>::infinity();
    double gradient_norm = std::numeric_limits<double>::infinity();
    Eigen::VectorXd n;
    Eigen::VectorXd gradient;
    Eigen::MatrixXd hessian;
};

bool is_finite_vector(const Eigen::VectorXd& values) {
    for (Eigen::Index i = 0; i < values.size(); ++i) {
        if (!std::isfinite(values[i])) {
            return false;
        }
    }
    return true;
}

std::vector<double> vector_from_eigen(const Eigen::VectorXd& values) {
    std::vector<double> out(static_cast<std::size_t>(values.size()), 0.0);
    for (Eigen::Index i = 0; i < values.size(); ++i) {
        out[static_cast<std::size_t>(i)] = values[i];
    }
    return out;
}

Eigen::VectorXd log_amounts_from_moles(const Eigen::VectorXd& n, double min_mole_fraction) {
    Eigen::VectorXd out(n.size());
    for (Eigen::Index i = 0; i < n.size(); ++i) {
        out[i] = std::log(std::max(n[i], min_mole_fraction));
    }
    return out;
}

Eigen::VectorXd soft_start_initial_moles(
    const std::vector<double>& initial,
    const Eigen::MatrixXd& balances,
    const Eigen::VectorXd& totals,
    double min_mole_fraction
) {
    Eigen::VectorXd n = Eigen::Map<const Eigen::VectorXd>(initial.data(), static_cast<Eigen::Index>(initial.size()));
    if (balances.rows() == 0) {
        return n;
    }

    Eigen::VectorXd balance_error = totals - balances * n;
    if (balance_error.cwiseAbs().maxCoeff() <= 1.0e-12) {
        return n;
    }

    Eigen::MatrixXd gram = balances * balances.transpose();
    const double scale = std::max(1.0, gram.cwiseAbs().maxCoeff());
    gram.diagonal().array() += 1.0e-12 * scale;
    Eigen::VectorXd correction_multiplier = gram.colPivHouseholderQr().solve(balance_error);
    Eigen::VectorXd candidate = n + balances.transpose() * correction_multiplier;
    if (!is_finite_vector(candidate)) {
        return n;
    }
    for (Eigen::Index i = 0; i < candidate.size(); ++i) {
        if (candidate[i] <= min_mole_fraction) {
            return n;
        }
    }
    Eigen::VectorXd projected_error = balances * candidate - totals;
    if (projected_error.cwiseAbs().maxCoeff() > 1.0e-8) {
        return n;
    }
    return candidate;
}

IdealGibbsEvaluation ideal_gibbs_soft_start_evaluate(
    const Eigen::VectorXd& n0,
    const Eigen::VectorXd& extents,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    double min_mole_fraction,
    double pressure_log
) {
    IdealGibbsEvaluation out;
    const Eigen::Index nreactions = reactions.rows();
    const Eigen::Index ncomp = reactions.cols();
    out.n = n0 + reactions.transpose() * extents;
    out.gradient = Eigen::VectorXd::Zero(nreactions);
    out.hessian = Eigen::MatrixXd::Zero(nreactions, nreactions);

    if (!is_finite_vector(out.n)) {
        return out;
    }
    double ntotal = 0.0;
    for (Eigen::Index i = 0; i < ncomp; ++i) {
        if (out.n[i] <= min_mole_fraction) {
            return out;
        }
        ntotal += out.n[i];
    }
    if (!std::isfinite(ntotal) || ntotal <= min_mole_fraction) {
        return out;
    }

    out.objective = -log_k.dot(extents);
    for (Eigen::Index i = 0; i < ncomp; ++i) {
        out.objective += out.n[i] * std::log(out.n[i]);
    }
    out.objective -= ntotal * std::log(ntotal);
    out.objective += ntotal * pressure_log;

    Eigen::VectorXd x = out.n / ntotal;
    Eigen::VectorXd reaction_sums = reactions.rowwise().sum();
    for (Eigen::Index r = 0; r < nreactions; ++r) {
        out.gradient[r] = -log_k[r] + reaction_sums[r] * pressure_log;
        for (Eigen::Index i = 0; i < ncomp; ++i) {
            out.gradient[r] += reactions(r, i) * std::log(std::max(x[i], min_mole_fraction));
        }
    }
    for (Eigen::Index r = 0; r < nreactions; ++r) {
        for (Eigen::Index s = 0; s < nreactions; ++s) {
            double value = -(reaction_sums[r] * reaction_sums[s]) / ntotal;
            for (Eigen::Index i = 0; i < ncomp; ++i) {
                value += reactions(r, i) * reactions(s, i) / out.n[i];
            }
            out.hessian(r, s) = value;
        }
    }
    out.gradient_norm = out.gradient.size() == 0 ? 0.0 : out.gradient.cwiseAbs().maxCoeff();
    out.feasible = std::isfinite(out.objective) && std::isfinite(out.gradient_norm) && is_finite_vector(out.gradient);
    return out;
}

ChemicalSoftStartResult ideal_gibbs_soft_start_solve(
    const std::vector<double>& initial,
    const Eigen::MatrixXd& balances,
    const Eigen::VectorXd& totals,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    const ChemicalEquilibriumOptionsNative& options,
    int phase_int,
    double p
) {
    ChemicalSoftStartResult out;
    out.composition = initial;
    out.log_n = log_amounts_from_moles(
        Eigen::Map<const Eigen::VectorXd>(initial.data(), static_cast<Eigen::Index>(initial.size())),
        options.min_mole_fraction
    );
    if (reactions.rows() == 0) {
        out.rejection_reason = "no_reactions";
        return out;
    }
    out.attempted = true;
    out.rejection_reason = "not_converged";

    Eigen::VectorXd n0 = soft_start_initial_moles(initial, balances, totals, options.min_mole_fraction);
    const double pressure_log = phase_int == 1 ? std::log(std::max(p, 1.0e-300) / 1.0e5) : 0.0;
    Eigen::VectorXd extents = Eigen::VectorXd::Zero(reactions.rows());
    IdealGibbsEvaluation current = ideal_gibbs_soft_start_evaluate(
        n0,
        extents,
        reactions,
        log_k,
        options.min_mole_fraction,
        pressure_log
    );
    if (!current.feasible) {
        out.rejection_reason = "initial_soft_start_infeasible";
        return out;
    }

    constexpr int max_soft_iterations = 40;
    constexpr double gradient_tolerance = 1.0e-10;
    for (int iteration = 0; iteration <= max_soft_iterations; ++iteration) {
        out.iterations = iteration;
        out.history.push_back(current.gradient_norm);
        if (current.gradient_norm <= gradient_tolerance) {
            out.success = true;
            out.rejection_reason = "";
            break;
        }
        if (iteration == max_soft_iterations) {
            break;
        }

        Eigen::MatrixXd hessian = current.hessian;
        const double hessian_scale = std::max(1.0, hessian.cwiseAbs().maxCoeff());
        hessian.diagonal().array() += 1.0e-12 * hessian_scale;
        Eigen::VectorXd step = hessian.colPivHouseholderQr().solve(-current.gradient);
        if (!is_finite_vector(step)) {
            out.rejection_reason = "singular_extent_step";
            break;
        }

        Eigen::VectorXd dn = reactions.transpose() * step;
        double alpha = 1.0;
        for (Eigen::Index i = 0; i < current.n.size(); ++i) {
            if (dn[i] < 0.0) {
                alpha = std::min(alpha, 0.99 * (current.n[i] - options.min_mole_fraction) / (-dn[i]));
            }
        }
        alpha = std::max(0.0, std::min(1.0, alpha));
        bool accepted = false;
        for (int trial = 0; trial < 24 && alpha > 1.0e-12; ++trial) {
            Eigen::VectorXd candidate_extents = extents + alpha * step;
            IdealGibbsEvaluation candidate = ideal_gibbs_soft_start_evaluate(
                n0,
                candidate_extents,
                reactions,
                log_k,
                options.min_mole_fraction,
                pressure_log
            );
            if (candidate.feasible && candidate.objective <= current.objective + 1.0e-14) {
                extents = candidate_extents;
                current = candidate;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            out.rejection_reason = "extent_line_search_failed";
            break;
        }
    }

    if (!out.success && current.feasible && !out.history.empty()) {
        out.success = true;
        out.rejection_reason = "";
    }
    if (out.success) {
        out.log_n = log_amounts_from_moles(current.n, options.min_mole_fraction);
        out.composition = composition_from_moles(vector_from_eigen(current.n), options.min_mole_fraction);
    }
    return out;
}

void add_soft_start_diagnostics(
    ChemicalEquilibriumResultNative& result,
    const ChemicalSoftStartResult& soft_start
) {
    result.diagnostics_bool["soft_start_enabled"] = soft_start.enabled;
    result.diagnostics_bool["soft_start_attempted"] = soft_start.attempted;
    result.diagnostics_bool["soft_start_success"] = soft_start.success;
    result.diagnostics_bool["soft_start_used"] = soft_start.used;
    result.diagnostics_int["soft_start_iterations"] = soft_start.iterations;
    result.diagnostics_double["soft_start_residual_norm"] = soft_start.residual_norm;
    result.diagnostics_double["soft_start_initial_residual_norm"] = soft_start.initial_residual_norm;
    result.diagnostics_string["soft_start_rejection_reason"] = soft_start.rejection_reason;
    result.diagnostics_vector["soft_start_composition"] = soft_start.composition;
    result.diagnostics_vector["soft_start_history"] = soft_start.history;
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
    int phase_int,
    const std::string& activity_model,
    int* state_failure_count,
    ChemicalEvaluationCounters* counters
) {
    ChemicalEvaluation out;
    if (counters != nullptr) {
        counters->residual_evaluations += 1;
    }
    out.n = moles_from_log_amounts(log_n);
    out.x = composition_from_moles(out.n, options.min_mole_fraction);
    const bool evaluate_activity = should_evaluate_activity_coefficients(reaction_standard_states, options);
    if (evaluate_activity) {
        if (counters != nullptr) {
            counters->activity_evaluations += 1;
            counters->state_evaluations += 1;
            counters->density_solves += 1;
        }
        out.gamma = activity_coefficients(
            mixture,
            t,
            p,
            out.x,
            phase_int,
            activity_model,
            options.min_mole_fraction,
            state_failure_count
        );
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

    double molar_density = 0.0;
    if (standard_states_need_concentration(reaction_standard_states)) {
        try {
            if (counters != nullptr) {
                counters->density_solves += 1;
                counters->state_evaluations += 1;
            }
            molar_density = mixture->solve_density_scoped(
                t,
                p,
                out.x,
                phase_int,
                "chemical_equilibrium_concentration_standard_state"
            );
        } catch (...) {
            if (state_failure_count != nullptr) {
                *state_failure_count += 1;
            }
            throw;
        }
    }

    Eigen::VectorXd reaction(reactions.rows());
    for (Eigen::Index r = 0; r < reactions.rows(); ++r) {
        double value = -log_k[r];
        const int standard_state = reaction_standard_states[static_cast<std::size_t>(r)];
        for (Eigen::Index i = 0; i < reactions.cols(); ++i) {
            double species_activity = out.x[static_cast<std::size_t>(i)];
            if (standard_state == STANDARD_STATE_IDEAL_MOLE_FRACTION) {
                species_activity = out.x[static_cast<std::size_t>(i)];
            } else if (standard_state == STANDARD_STATE_CONCENTRATION) {
                species_activity = out.x[static_cast<std::size_t>(i)] * molar_density;
            } else if (standard_state == STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
                if (out.gamma.size() != out.x.size()) {
                    throw ValueError("activity-coupled reaction residual requires activity coefficients.");
                }
                species_activity = out.x[static_cast<std::size_t>(i)] * out.gamma[static_cast<std::size_t>(i)];
            } else {
                throw ValueError("reaction standard state contains an unsupported code.");
            }
            value += reactions(r, i) * std::log(std::max(species_activity, options.min_mole_fraction));
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

Eigen::MatrixXd chemical_residual_jacobian(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const Eigen::VectorXd& log_n,
    const ChemicalEvaluation& current,
    const Eigen::MatrixXd& balances,
    const Eigen::VectorXd& totals,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options,
    int phase_int,
    const std::string& activity_model,
    int* state_failure_count,
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
    (void)mixture;
    (void)t;
    (void)p;
    (void)log_n;
    (void)totals;
    (void)log_k;
    (void)phase_int;
    (void)activity_model;
    (void)state_failure_count;
    throw ValueError(selected.not_available_reason.empty()
        ? "not_available: chemical-equilibrium residual jacobian is unavailable."
        : selected.not_available_reason);
}

bool is_scalar_binary_activity_problem(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const Eigen::MatrixXd& balances,
    const Eigen::VectorXd& totals,
    const Eigen::MatrixXd& reactions,
    const std::vector<int>& reaction_standard_states
) {
    if (has_ionic_species(mixture) || mixture->args().m.size() != 2 || balances.rows() != 1 || balances.cols() != 2) {
        return false;
    }
    if (totals.size() != 1 || std::abs(totals[0] - 1.0) > 1.0e-12) {
        return false;
    }
    if (reactions.rows() != 1 || reactions.cols() != 2 || reaction_standard_states.size() != 1) {
        return false;
    }
    if (reaction_standard_states[0] != STANDARD_STATE_MOLE_FRACTION_ACTIVITY) {
        return false;
    }
    return std::abs(balances(0, 0) - 1.0) <= 1.0e-12
        && std::abs(balances(0, 1) - 1.0) <= 1.0e-12
        && std::abs(reactions(0, 0) + 1.0) <= 1.0e-12
        && std::abs(reactions(0, 1) - 1.0) <= 1.0e-12;
}

Eigen::VectorXd binary_log_amounts(double x0, double min_mole_fraction) {
    Eigen::VectorXd log_n(2);
    const double clipped = std::max(min_mole_fraction, std::min(1.0 - min_mole_fraction, x0));
    log_n[0] = std::log(clipped);
    log_n[1] = std::log(std::max(min_mole_fraction, 1.0 - clipped));
    return log_n;
}

std::vector<double> matrix_to_row_major_vector(const Eigen::MatrixXd& matrix) {
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(matrix.rows() * matrix.cols()));
    for (Eigen::Index r = 0; r < matrix.rows(); ++r) {
        for (Eigen::Index c = 0; c < matrix.cols(); ++c) {
            out.push_back(matrix(r, c));
        }
    }
    return out;
}

void add_reactive_speciation_implicit_sensitivity_diagnostics(
    ChemicalEquilibriumResultNative& result,
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const Eigen::VectorXd& log_n,
    const ChemicalEvaluation& current,
    const Eigen::MatrixXd& balances,
    const Eigen::VectorXd& totals,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options,
    int phase_int,
    const std::string& activity_model
) {
    try {
        int state_failure_count = 0;
        Eigen::MatrixXd residual_state_jacobian = chemical_residual_jacobian(
            mixture,
            t,
            p,
            log_n,
            current,
            balances,
            totals,
            reactions,
            log_k,
            reaction_standard_states,
            options,
            phase_int,
            activity_model,
            &state_failure_count,
            nullptr,
            nullptr
        );
        Eigen::MatrixXd residual_parameter_jacobian =
            Eigen::MatrixXd::Zero(residual_state_jacobian.rows(), reactions.rows());
        const Eigen::Index reaction_offset = balances.rows() + 1;
        for (Eigen::Index r = 0; r < reactions.rows(); ++r) {
            residual_parameter_jacobian(reaction_offset + r, r) = -1.0;
        }
        Eigen::MatrixXd sensitivity =
            residual_state_jacobian.colPivHouseholderQr().solve(-residual_parameter_jacobian);

        result.diagnostics_string["implicit_sensitivity_backend"] = "analytic_implicit";
        result.diagnostics_string["implicit_sensitivity_status"] = "residual_jacobian_available";
        result.diagnostics_string["reactive_speciation_sensitivity_parameter"] = "log_equilibrium_constants";
        result.diagnostics_int["reactive_speciation_residual_rows"] =
            static_cast<int>(residual_state_jacobian.rows());
        result.diagnostics_int["reactive_speciation_state_size"] =
            static_cast<int>(residual_state_jacobian.cols());
        result.diagnostics_int["reactive_speciation_parameter_size"] =
            static_cast<int>(reactions.rows());
        result.diagnostics_vector["reactive_speciation_state"] = vector_from_eigen(log_n);
        result.diagnostics_vector["reactive_speciation_residual"] = current.residuals;
        result.diagnostics_vector["reactive_speciation_residual_state_jacobian_row_major"] =
            matrix_to_row_major_vector(residual_state_jacobian);
        result.diagnostics_vector["reactive_speciation_residual_parameter_jacobian_row_major"] =
            matrix_to_row_major_vector(residual_parameter_jacobian);
        result.diagnostics_vector["reactive_speciation_log_amount_sensitivity_to_log_k_row_major"] =
            matrix_to_row_major_vector(sensitivity);
    } catch (const std::exception& exc) {
        result.diagnostics_string["implicit_sensitivity_backend"] = "not_available";
        result.diagnostics_string["implicit_sensitivity_status"] = "not_available";
        result.diagnostics_string["implicit_sensitivity_reason"] = exc.what();
    }
}

bool try_scalar_binary_activity_solve_native(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    double t,
    double p,
    const Eigen::MatrixXd& balances,
    const Eigen::VectorXd& totals,
    const Eigen::MatrixXd& reactions,
    const Eigen::VectorXd& log_k,
    const std::vector<int>& reaction_standard_states,
    const ChemicalEquilibriumOptionsNative& options,
    int phase_int,
    const std::string& activity_model,
    const ChemicalDerivativeSelection& derivative_selection,
    double preferred_fraction,
    ChemicalEquilibriumResultNative& result
) {
    if (!is_scalar_binary_activity_problem(mixture, balances, totals, reactions, reaction_standard_states)) {
        return false;
    }
    int state_failure_count = 0;
    ChemicalEvaluationCounters counters;
    std::vector<double> history;
    const double low = options.min_mole_fraction;
    const double high = 1.0 - options.min_mole_fraction;

    auto evaluate_fraction = [&](double x0) {
        return evaluate_chemical(
            mixture,
            t,
            p,
            binary_log_amounts(x0, options.min_mole_fraction),
            balances,
            totals,
            reactions,
            log_k,
            reaction_standard_states,
            options,
            phase_int,
            activity_model,
            &state_failure_count,
            &counters
        );
    };

    ChemicalEvaluation best;
    double best_residual = std::numeric_limits<double>::infinity();
    double best_fraction = std::numeric_limits<double>::infinity();
    int iteration = 0;
    const int max_iterations = std::max(1, options.max_iterations);
    const int scan_points = std::max(64, max_iterations * 4);

    auto consider_candidate = [&](double fraction, const ChemicalEvaluation& evaluation, double residual) {
        const double preferred_distance = std::abs(fraction - preferred_fraction);
        const double best_distance = std::abs(best_fraction - preferred_fraction);
        if (!std::isfinite(best_residual) || preferred_distance < best_distance
            || (std::abs(preferred_distance - best_distance) <= 1.0e-14
                && std::abs(residual) < std::abs(best_residual))) {
            best = evaluation;
            best_residual = residual;
            best_fraction = fraction;
        }
    };

    double previous_x = low;
    ChemicalEvaluation previous_eval;
    double previous_residual = std::numeric_limits<double>::quiet_NaN();
    try {
        previous_eval = evaluate_fraction(previous_x);
        previous_residual = previous_eval.reaction_residuals.empty()
            ? previous_eval.residual_norm
            : previous_eval.reaction_residuals[0];
    } catch (...) {
        return false;
    }
    if (!std::isfinite(previous_residual)) {
        return false;
    }
    history.push_back(std::abs(previous_residual));
    if (std::abs(previous_residual) <= options.tolerance) {
        consider_candidate(previous_x, previous_eval, previous_residual);
    }

    for (int scan = 1; scan <= scan_points; ++scan) {
        const double current_x = low + (high - low) * static_cast<double>(scan) / static_cast<double>(scan_points);
        ChemicalEvaluation current_eval;
        double current_residual = std::numeric_limits<double>::quiet_NaN();
        try {
            current_eval = evaluate_fraction(current_x);
            current_residual = current_eval.reaction_residuals.empty()
                ? current_eval.residual_norm
                : current_eval.reaction_residuals[0];
        } catch (...) {
            previous_x = current_x;
            continue;
        }
        if (!std::isfinite(current_residual)) {
            previous_x = current_x;
            continue;
        }
        history.push_back(std::abs(current_residual));
        if (std::abs(current_residual) <= options.tolerance) {
            consider_candidate(current_x, current_eval, current_residual);
        }
        if (previous_residual * current_residual <= 0.0) {
            double left = previous_x;
            double right = current_x;
            double f_left = previous_residual;
            ChemicalEvaluation bracket_best = std::abs(previous_residual) <= std::abs(current_residual)
                ? previous_eval
                : current_eval;
            double bracket_residual = std::abs(previous_residual) <= std::abs(current_residual)
                ? previous_residual
                : current_residual;
            double bracket_fraction = std::abs(previous_residual) <= std::abs(current_residual)
                ? previous_x
                : current_x;
            for (int inner = 1; inner <= max_iterations; ++inner) {
                ++iteration;
                const double mid = 0.5 * (left + right);
                ChemicalEvaluation mid_eval = evaluate_fraction(mid);
                const double mid_residual = mid_eval.reaction_residuals.empty()
                    ? mid_eval.residual_norm
                    : mid_eval.reaction_residuals[0];
                history.push_back(std::abs(mid_residual));
                if (std::abs(mid_residual) < std::abs(bracket_residual)) {
                    bracket_best = mid_eval;
                    bracket_residual = mid_residual;
                    bracket_fraction = mid;
                }
                if (std::abs(mid_residual) <= options.tolerance) {
                    break;
                }
                if (f_left * mid_residual <= 0.0) {
                    right = mid;
                } else {
                    left = mid;
                    f_left = mid_residual;
                }
            }
            consider_candidate(bracket_fraction, bracket_best, bracket_residual);
        }
        previous_x = current_x;
        previous_eval = current_eval;
        previous_residual = current_residual;
    }
    if (!std::isfinite(best_residual)) {
        return false;
    }

    result.success = std::abs(best_residual) <= options.tolerance;
    result.message = result.success ? "converged" : "scalar binary activity solve did not converge";
    result.composition = best.x;
    result.activity_coefficients = best.gamma;
    result.mass_balance_residuals = best.mass_residuals;
    result.charge_residual = best.charge_residual;
    result.reaction_residuals = best.reaction_residuals;
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_chemical_equilibrium_native";
    result.diagnostics_string["problem_class"] = "homogeneous_chemical_equilibrium";
    result.diagnostics_string["activity_model"] = activity_model;
    result.diagnostics_string["activity_output"] = options.activity_output;
    result.diagnostics_string["activity_basis"] = standard_state_summary(reaction_standard_states);
    result.diagnostics_string["phase"] = options.phase;
    result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
    result.diagnostics_string["jacobian_backend"] = derivative_selection.backend;
    result.diagnostics_string["derivative_backend"] = derivative_selection.backend;
    result.diagnostics_string["derivative_status"] = derivative_selection.backend;
    result.diagnostics_string["derivative_capability_path"] = derivative_selection.capability_path;
    result.diagnostics_string["not_available_reason"] = derivative_selection.not_available_reason;
    result.diagnostics_string["selected_solver_backend"] = "native_scalar_binary_activity_bracket";
    result.diagnostics_string["activity_derivative_policy"] = "native_scalar_activity_root";
    result.diagnostics_bool["derivative_available"] = derivative_selection.derivative_available;
    result.diagnostics_bool["jacobian_available"] = derivative_selection.derivative_available;
    result.diagnostics_bool["jacobian_fallback_used"] = false;
    result.diagnostics_string["jacobian_fallback_reason"] = "";
    result.diagnostics_bool["activity_fixed_point"] = false;
    result.diagnostics_int["activity_fixed_point_updates"] = 0;
    result.diagnostics_bool["activity_or_fugacity_terms_in_residual"] = true;
    result.diagnostics_bool["activity_derivative_in_jacobian"] = false;
    result.diagnostics_bool["hessian_available"] = false;
    result.diagnostics_string["hessian_backend"] = "not_implemented";
    result.diagnostics_bool["hessian_fallback_used"] = false;
    result.diagnostics_string["hessian_fallback_reason"] =
        "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
    result.diagnostics_int["iterations"] = iteration;
    result.diagnostics_int["state_failure_count"] = state_failure_count;
    result.diagnostics_int["residual_evaluation_count"] = counters.residual_evaluations;
    result.diagnostics_int["jacobian_evaluation_count"] = 0;
    result.diagnostics_int["state_evaluation_count"] = counters.state_evaluations;
    result.diagnostics_int["activity_evaluation_count"] = counters.activity_evaluations;
    result.diagnostics_int["density_solve_count"] = counters.density_solves;
    result.diagnostics_bool["activity_coefficients_evaluated"] = !best.gamma.empty();
    result.diagnostics_double["residual_norm"] = best.residual_norm;
    result.diagnostics_double["tolerance"] = options.tolerance;
    result.diagnostics_vector["history"] = history;
    result.diagnostics_vector["phase_handoff_composition"] = best.x;
    add_reactive_speciation_implicit_sensitivity_diagnostics(
        result,
        mixture,
        t,
        p,
        binary_log_amounts(best_fraction, options.min_mole_fraction),
        best,
        balances,
        totals,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        phase_int,
        activity_model
    );
    ChemicalSoftStartResult no_soft_start;
    no_soft_start.enabled = false;
    add_soft_start_diagnostics(result, no_soft_start);
    return true;
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

    const int phase_int = phase_token_to_int_chemical(options.phase);
    const std::string activity_model = activity_model_for_standard_states(mixture, reaction_standard_states);
    ChemicalDerivativeSelection derivative_selection = select_chemical_derivative_backend(
        options,
        reaction_standard_states
    );
    int state_failure_count = 0;
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
        phase_int,
        activity_model,
        &state_failure_count,
        &counters
    );
    Eigen::MatrixXd jac = chemical_residual_jacobian(
        mixture,
        t,
        p,
        log_n,
        current,
        balances,
        totals,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        phase_int,
        activity_model,
        &state_failure_count,
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
    out.diagnostics_string["derivative_status"] = derivative_selection.backend;
    out.diagnostics_string["derivative_capability_path"] = derivative_selection.capability_path;
    out.diagnostics_string["not_available_reason"] = derivative_selection.not_available_reason;
    out.diagnostics_string["hessian_backend"] = "gauss_newton";
    out.diagnostics_string["hessian_kind"] = "approximate_least_squares_gauss_newton";
    out.diagnostics_string["hessian_structure"] = "dense_lower_triangular";
    out.diagnostics_bool["derivative_available"] = derivative_selection.derivative_available;
    out.diagnostics_bool["jacobian_available"] = derivative_selection.derivative_available;
    out.diagnostics_bool["hessian_available"] = true;
    out.diagnostics_bool["exact_hessian_available"] = false;
    out.diagnostics_bool["hessian_callback_available"] = true;
    out.diagnostics_bool["hessian_includes_second_residual_derivatives"] = false;
    out.diagnostics_bool["sparse_hessian_available"] = false;
    out.diagnostics_bool["activity_coefficients_evaluated"] = !current.gamma.empty();
    out.diagnostics_int["state_failure_count"] = state_failure_count;
    out.diagnostics_int["residual_evaluation_count"] = counters.residual_evaluations;
    out.diagnostics_int["jacobian_evaluation_count"] = counters.jacobian_evaluations;
    out.diagnostics_int["state_evaluation_count"] = counters.state_evaluations;
    out.diagnostics_int["activity_evaluation_count"] = counters.activity_evaluations;
    out.diagnostics_int["density_solve_count"] = counters.density_solves;
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
    if (options.max_iterations < 0 || options.tolerance <= 0.0 || options.damping <= 0.0 || options.damping > 1.0
        || options.min_mole_fraction <= 0.0) {
        throw ValueError("chemical equilibrium options contain invalid numerical controls.");
    }
    if (options.solver_backend != "auto" && options.solver_backend != "ipopt") {
        throw ValueError("chemical equilibrium solver_backend must be 'auto' or 'ipopt'.");
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
    std::vector<double> initial = normalize_composition_chemical(initial_x, options.min_mole_fraction);
    Eigen::VectorXd log_n(static_cast<Eigen::Index>(ncomp));
    for (int i = 0; i < ncomp; ++i) {
        log_n[i] = std::log(std::max(initial[static_cast<std::size_t>(i)], options.min_mole_fraction));
    }

    const int phase_int = phase_token_to_int_chemical(options.phase);
    const std::string activity_model = activity_model_for_standard_states(mixture, reaction_standard_states);
    ChemicalDerivativeSelection derivative_selection = select_chemical_derivative_backend(
        options,
        reaction_standard_states
    );
    if (options.solver_backend == "ipopt") {
        if (!standard_states_all_ideal_mole_fraction(reaction_standard_states)) {
            throw ValueError(
                "Native Ipopt reactive speciation currently supports ideal_mole_fraction standard states; "
                "activity and concentration routes require the EOS derivative NLP blocks."
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
        request.initial_x = initial_x;
        request.charges = mixture->args().z;
        request.min_mole_fraction = options.min_mole_fraction;
        return epcsaft::native::equilibrium_nlp::solve_ideal_speciation_chemical_equilibrium_ipopt(
            request,
            options
        );
    }
    ChemicalEquilibriumResultNative scalar_result;
    if (try_scalar_binary_activity_solve_native(
        mixture,
        t,
        p,
        balances,
        totals,
        reactions,
        log_k,
        reaction_standard_states,
        options,
        phase_int,
        activity_model,
        derivative_selection,
        initial[0],
        scalar_result
    )) {
        return scalar_result;
    }
    int state_failure_count = 0;
    ChemicalEvaluationCounters counters;
    ChemicalSoftStartResult soft_start = ideal_gibbs_soft_start_solve(
        initial,
        balances,
        totals,
        reactions,
        log_k,
        options,
        phase_int,
        p
    );
    ChemicalEvaluation initial_evaluation = evaluate_chemical(
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
        phase_int,
        activity_model,
        &state_failure_count,
        &counters
    );
    soft_start.initial_residual_norm = initial_evaluation.residual_norm;
    if (soft_start.composition.empty()) {
        soft_start.composition = initial_evaluation.x;
    }
    if (!std::isfinite(soft_start.residual_norm)) {
        soft_start.residual_norm = initial_evaluation.residual_norm;
    }
    if (soft_start.success) {
        try {
            ChemicalEvaluation soft_evaluation = evaluate_chemical(
                mixture,
                t,
                p,
                soft_start.log_n,
                balances,
                totals,
                reactions,
                log_k,
                reaction_standard_states,
                options,
                phase_int,
                activity_model,
                &state_failure_count,
                &counters
            );
            soft_start.residual_norm = soft_evaluation.residual_norm;
            soft_start.composition = soft_evaluation.x;
            if (std::isfinite(soft_evaluation.residual_norm)
                && soft_evaluation.residual_norm <= initial_evaluation.residual_norm) {
                log_n = soft_start.log_n;
                soft_start.used = true;
                soft_start.rejection_reason = "";
            } else {
                soft_start.used = false;
                soft_start.rejection_reason = "full_residual_not_improved";
            }
        } catch (...) {
            soft_start.used = false;
            soft_start.success = false;
            soft_start.rejection_reason = "full_residual_evaluation_failed";
        }
    }
    std::vector<double> history;
    ChemicalEvaluation best;

    for (int iteration = 0; iteration <= options.max_iterations; ++iteration) {
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
            phase_int,
            activity_model,
            &state_failure_count,
            &counters
        );
        best = current;
        history.push_back(current.residual_norm);
        if (current.residual_norm <= options.tolerance) {
            ChemicalEquilibriumResultNative result;
            result.success = true;
            result.message = "converged";
            result.composition = current.x;
            result.activity_coefficients = current.gamma;
            result.mass_balance_residuals = current.mass_residuals;
            result.charge_residual = current.charge_residual;
            result.reaction_residuals = current.reaction_residuals;
            result.diagnostics_string["solver_language"] = "c++";
            result.diagnostics_string["native_entrypoint"] = "_solve_chemical_equilibrium_native";
            result.diagnostics_string["problem_class"] = "homogeneous_chemical_equilibrium";
            result.diagnostics_string["activity_model"] = activity_model;
            result.diagnostics_string["activity_output"] = options.activity_output;
            result.diagnostics_string["activity_basis"] = standard_state_summary(reaction_standard_states);
            result.diagnostics_string["phase"] = options.phase;
            result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
            result.diagnostics_string["jacobian_backend"] = derivative_selection.backend;
            result.diagnostics_string["derivative_backend"] = derivative_selection.backend;
            result.diagnostics_string["derivative_status"] = derivative_selection.backend;
            result.diagnostics_string["derivative_capability_path"] = derivative_selection.capability_path;
            result.diagnostics_string["not_available_reason"] = derivative_selection.not_available_reason;
            result.diagnostics_bool["derivative_available"] = derivative_selection.derivative_available;
            result.diagnostics_bool["jacobian_available"] = derivative_selection.derivative_available;
            result.diagnostics_bool["jacobian_fallback_used"] = false;
            result.diagnostics_string["jacobian_fallback_reason"] = "";
            const bool coupled_standard_states = !standard_states_all_ideal_mole_fraction(reaction_standard_states);
            result.diagnostics_bool["activity_fixed_point"] = false;
            result.diagnostics_int["activity_fixed_point_updates"] = 0;
            result.diagnostics_bool["activity_or_fugacity_terms_in_residual"] = coupled_standard_states;
            result.diagnostics_bool["activity_derivative_in_jacobian"] = coupled_standard_states;
            result.diagnostics_string["activity_derivative_policy"] = coupled_standard_states
                ? "native_log_amount_residual_jacobian"
                : "not_required_for_ideal_mole_fraction";
            result.diagnostics_bool["hessian_available"] = false;
            result.diagnostics_string["hessian_backend"] = "not_implemented";
            result.diagnostics_bool["hessian_fallback_used"] = false;
            result.diagnostics_string["hessian_fallback_reason"] =
                "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
            result.diagnostics_int["iterations"] = iteration;
            result.diagnostics_int["state_failure_count"] = state_failure_count;
            result.diagnostics_int["residual_evaluation_count"] = counters.residual_evaluations;
            result.diagnostics_int["jacobian_evaluation_count"] = counters.jacobian_evaluations;
            result.diagnostics_int["state_evaluation_count"] = counters.state_evaluations;
            result.diagnostics_int["activity_evaluation_count"] = counters.activity_evaluations;
            result.diagnostics_int["density_solve_count"] = counters.density_solves;
            result.diagnostics_bool["activity_coefficients_evaluated"] = !current.gamma.empty();
            result.diagnostics_double["residual_norm"] = current.residual_norm;
            result.diagnostics_double["tolerance"] = options.tolerance;
            result.diagnostics_vector["history"] = history;
            result.diagnostics_vector["phase_handoff_composition"] = current.x;
            add_reactive_speciation_implicit_sensitivity_diagnostics(
                result,
                mixture,
                t,
                p,
                log_n,
                current,
                balances,
                totals,
                reactions,
                log_k,
                reaction_standard_states,
                options,
                phase_int,
                activity_model
            );
            add_soft_start_diagnostics(result, soft_start);
            return result;
        }
        if (iteration == options.max_iterations) {
            break;
        }
        Eigen::MatrixXd jac = chemical_residual_jacobian(
            mixture,
            t,
            p,
            log_n,
            current,
            balances,
            totals,
            reactions,
            log_k,
            reaction_standard_states,
            options,
            phase_int,
            activity_model,
            &state_failure_count,
            &counters,
            &derivative_selection
        );
        Eigen::VectorXd residual = Eigen::Map<const Eigen::VectorXd>(
            current.residuals.data(),
            static_cast<Eigen::Index>(current.residuals.size())
        );
        Eigen::VectorXd step = jac.colPivHouseholderQr().solve(-residual);
        double alpha = options.damping;
        bool accepted = false;
        for (int trial = 0; trial < 12; ++trial) {
            Eigen::VectorXd candidate_log_n = log_n + alpha * step;
            ChemicalEvaluation candidate = evaluate_chemical(
                mixture,
                t,
                p,
                candidate_log_n,
                balances,
                totals,
                reactions,
                log_k,
                reaction_standard_states,
                options,
                phase_int,
                activity_model,
                &state_failure_count,
                &counters
            );
            if (candidate.residual_norm <= current.residual_norm || alpha <= 1.0e-4) {
                log_n = candidate_log_n;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!accepted) {
            log_n += options.damping * step;
        }
    }

    ChemicalEquilibriumResultNative result;
    result.success = false;
    result.message = "chemical equilibrium did not converge";
    result.composition = best.x;
    result.activity_coefficients = best.gamma;
    result.mass_balance_residuals = best.mass_residuals;
    result.charge_residual = best.charge_residual;
    result.reaction_residuals = best.reaction_residuals;
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_chemical_equilibrium_native";
    result.diagnostics_string["problem_class"] = "homogeneous_chemical_equilibrium";
    result.diagnostics_string["activity_model"] = activity_model;
    result.diagnostics_string["activity_output"] = options.activity_output;
    result.diagnostics_string["activity_basis"] = standard_state_summary(reaction_standard_states);
    result.diagnostics_string["phase"] = options.phase;
    result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
    result.diagnostics_string["jacobian_backend"] = derivative_selection.backend;
    result.diagnostics_string["derivative_backend"] = derivative_selection.backend;
    result.diagnostics_string["derivative_status"] = derivative_selection.backend;
    result.diagnostics_string["derivative_capability_path"] = derivative_selection.capability_path;
    result.diagnostics_string["not_available_reason"] = derivative_selection.not_available_reason;
    result.diagnostics_bool["derivative_available"] = derivative_selection.derivative_available;
    result.diagnostics_bool["jacobian_available"] = derivative_selection.derivative_available;
    result.diagnostics_bool["jacobian_fallback_used"] = false;
    result.diagnostics_string["jacobian_fallback_reason"] = "";
    result.diagnostics_bool["hessian_available"] = false;
    result.diagnostics_string["hessian_backend"] = "not_implemented";
    result.diagnostics_bool["hessian_fallback_used"] = false;
    result.diagnostics_string["hessian_fallback_reason"] =
        "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
    result.diagnostics_int["iterations"] = options.max_iterations;
    result.diagnostics_int["state_failure_count"] = state_failure_count;
    result.diagnostics_int["residual_evaluation_count"] = counters.residual_evaluations;
    result.diagnostics_int["jacobian_evaluation_count"] = counters.jacobian_evaluations;
    result.diagnostics_int["state_evaluation_count"] = counters.state_evaluations;
    result.diagnostics_int["activity_evaluation_count"] = counters.activity_evaluations;
    result.diagnostics_int["density_solve_count"] = counters.density_solves;
    result.diagnostics_bool["activity_coefficients_evaluated"] = !best.gamma.empty();
    result.diagnostics_double["residual_norm"] = best.residual_norm;
    result.diagnostics_double["tolerance"] = options.tolerance;
    result.diagnostics_vector["history"] = history;
    result.diagnostics_vector["phase_handoff_composition"] = best.x;
    add_soft_start_diagnostics(result, soft_start);
    return result;
}
