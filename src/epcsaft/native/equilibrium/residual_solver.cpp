#include "residual_solver.h"

#include "equilibrium_helpers.h"

#ifdef EPCSAFT_HAS_CERES
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
#endif

#include <algorithm>
#include <cmath>
#include <sstream>

namespace epcsaft::native::equilibrium {

namespace {

#ifdef EPCSAFT_HAS_CERES
std::string ceres_termination_type_name(ceres::TerminationType type) {
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

class NativeResidualCeresCostFunction final : public ceres::CostFunction {
public:
    explicit NativeResidualCeresCostFunction(const NativeResidualProblem& problem)
        : problem_(problem) {
        set_num_residuals(problem_.residual_count());
        mutable_parameter_block_sizes()->push_back(problem_.variable_count());
    }

    bool Evaluate(double const* const* parameters, double* residuals, double** jacobians) const override {
        try {
            std::vector<double> variables(parameters[0], parameters[0] + parameter_block_sizes()[0]);
            const bool need_jacobian = jacobians != nullptr && jacobians[0] != nullptr;
            NativeResidualEvaluation eval = problem_.evaluate(variables, need_jacobian);
            if (!eval.success || eval.residual.size() != static_cast<std::size_t>(num_residuals())) {
                return false;
            }
            std::copy(eval.residual.begin(), eval.residual.end(), residuals);
            if (need_jacobian) {
                if (eval.jacobian_row_major.size()
                    != static_cast<std::size_t>(num_residuals() * parameter_block_sizes()[0])) {
                    return false;
                }
                std::copy(eval.jacobian_row_major.begin(), eval.jacobian_row_major.end(), jacobians[0]);
            }
            return true;
        } catch (...) {
            return false;
        }
    }

private:
    const NativeResidualProblem& problem_;
};
#endif

double residual_cost(const std::vector<double>& residual) {
    double cost = 0.0;
    for (double value : residual) {
        cost += 0.5 * value * value;
    }
    return cost;
}

bool residuals_finite(const std::vector<double>& residual) {
    return std::all_of(residual.begin(), residual.end(), [](double value) {
        return std::isfinite(value);
    });
}

}  // namespace

NativeResidualSolveResult solve_native_residual_problem(
    const NativeResidualProblem& problem,
    const std::vector<double>& initial_variables,
    const NativeResidualSolveOptions& options
) {
    if (initial_variables.size() != static_cast<std::size_t>(problem.variable_count())) {
        throw SolutionError("Native residual solver received an invalid variable vector size.");
    }
    NativeResidualSolveResult out;
    out.variables = initial_variables;
#ifndef EPCSAFT_HAS_CERES
    (void)problem;
    (void)options;
    out.solver_backend = "none";
    out.solver_method = "none";
    out.termination_type = "ceres_not_available";
    out.rejection_reason = "Ceres support is required for native residual solving.";
    throw SolutionError(out.rejection_reason);
#else
    try {
        NativeResidualEvaluation initial = problem.evaluate(out.variables, false);
        out.initial_cost = residual_cost(initial.residual);
    } catch (const std::exception& exc) {
        out.rejection_reason = exc.what();
        out.final_cost = 1.0e300;
        out.residual_norm_linf = 1.0e300;
        out.residual_norm_l2 = 1.0e300;
        return out;
    }

    ceres::Problem ceres_problem;
    auto* cost = new NativeResidualCeresCostFunction(problem);
    ceres_problem.AddResidualBlock(cost, nullptr, out.variables.data());
    for (int index = 0; index < static_cast<int>(out.variables.size()); ++index) {
        ceres_problem.SetParameterLowerBound(out.variables.data(), index, options.lower_bound);
        ceres_problem.SetParameterUpperBound(out.variables.data(), index, options.upper_bound);
    }

    ceres::Solver::Options ceres_options;
    ceres_options.linear_solver_type = ceres::DENSE_QR;
    ceres_options.max_num_iterations = options.max_iterations;
    ceres_options.minimizer_progress_to_stdout = false;
    ceres_options.logging_type = ceres::SILENT;
    ceres_options.function_tolerance = options.function_tolerance;
    ceres_options.gradient_tolerance = options.gradient_tolerance;
    ceres_options.parameter_tolerance = options.parameter_tolerance;

    ceres::Solver::Summary summary;
    ceres::Solve(ceres_options, &ceres_problem, &summary);
    out.status = static_cast<int>(summary.termination_type);
    out.termination_type = ceres_termination_type_name(summary.termination_type);
    out.summary = summary.BriefReport();
    out.iterations = static_cast<int>(summary.iterations.size());
    out.residual_evaluations = static_cast<int>(summary.num_residual_evaluations);
    out.jacobian_evaluations = static_cast<int>(summary.num_jacobian_evaluations);
    out.solver_usable = summary.IsSolutionUsable();
    out.termination_accepted = summary.termination_type == ceres::CONVERGENCE;

    try {
        out.final_evaluation = problem.evaluate(out.variables, true);
        out.final_cost = residual_cost(out.final_evaluation.residual);
        out.residual_norm_linf = max_abs(out.final_evaluation.residual);
        out.residual_norm_l2 = l2_norm(out.final_evaluation.residual);
        out.residuals_accepted =
            out.final_evaluation.success
            && residuals_finite(out.final_evaluation.residual)
            && out.residual_norm_linf <= options.residual_tolerance;
    } catch (const std::exception& exc) {
        out.rejection_reason = exc.what();
        out.final_cost = 1.0e300;
        out.residual_norm_linf = 1.0e300;
        out.residual_norm_l2 = 1.0e300;
    }
    out.accepted = out.solver_usable && out.termination_accepted && out.residuals_accepted;
    if (!out.accepted && out.rejection_reason.empty()) {
        std::ostringstream reason;
        reason << "native residual solve was not accepted"
               << ": termination=" << out.termination_type
               << ", solver_usable=" << (out.solver_usable ? "true" : "false")
               << ", residual_norm_linf=" << out.residual_norm_linf
               << ", residual_tolerance=" << options.residual_tolerance;
        out.rejection_reason = reason.str();
    }
    return out;
#endif
}

}  // namespace epcsaft::native::equilibrium
