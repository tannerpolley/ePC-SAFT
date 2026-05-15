#pragma once

#include <map>
#include <string>
#include <vector>

namespace epcsaft::native::equilibrium {

struct NativeResidualSolveOptions {
    int max_iterations = 80;
    double residual_tolerance = 1.0e-8;
    double function_tolerance = 1.0e-12;
    double gradient_tolerance = 1.0e-10;
    double parameter_tolerance = 1.0e-10;
    double lower_bound = -100.0;
    double upper_bound = 100.0;
};

struct NativeResidualEvaluation {
    bool success = false;
    std::vector<double> residual;
    std::vector<double> jacobian_row_major;
    int rows = 0;
    int cols = 0;
    std::vector<std::string> residual_names;
    std::vector<std::string> variable_names;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
};

class NativeResidualProblem {
public:
    virtual ~NativeResidualProblem() = default;
    virtual int variable_count() const = 0;
    virtual int residual_count() const = 0;
    virtual NativeResidualEvaluation evaluate(const std::vector<double>& variables, bool need_jacobian) const = 0;
};

struct NativeResidualSolveResult {
    bool solver_usable = false;
    bool termination_accepted = false;
    bool residuals_accepted = false;
    bool accepted = false;
    std::string solver_backend = "ceres";
    std::string solver_method = "ceres_trust_region_residual_solve";
    std::string termination_type = "not_started";
    std::string rejection_reason;
    std::string summary;
    int status = 0;
    int iterations = 0;
    int residual_evaluations = 0;
    int jacobian_evaluations = 0;
    double initial_cost = 0.0;
    double final_cost = 0.0;
    double residual_norm_linf = 0.0;
    double residual_norm_l2 = 0.0;
    std::vector<double> variables;
    NativeResidualEvaluation final_evaluation;
};

NativeResidualSolveResult solve_native_residual_problem(
    const NativeResidualProblem& problem,
    const std::vector<double>& initial_variables,
    const NativeResidualSolveOptions& options
);

}  // namespace epcsaft::native::equilibrium
