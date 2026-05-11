#pragma once

#include <map>
#include <string>
#include <vector>

enum class NativeRegressionStatus {
    CONVERGED = 0,
    MAX_ITERATIONS = 1,
    LINE_SEARCH_FAILED = 2,
    SINGULAR_JACOBIAN = 3,
    ALL_ROWS_FAILED = 4,
    NONFINITE_OBJECTIVE = 5,
    BOUNDS_INCONSISTENT = 6,
    INVALID_INPUT = 7,
    BACKEND_UNAVAILABLE = 8,
};

struct NativeRegressionParameterSpec {
    std::string name;
    std::string path;
    std::string kind;
    double initial = 0.0;
    double lower = 0.0;
    double upper = 0.0;
    double scale = 1.0;
    bool fixed = false;
    std::map<std::string, std::string> metadata;
};

struct NativeRegressionResidualSchemaEntry {
    std::string name;
    std::string row_id;
    std::string family;
    std::string target;
    double scale = 1.0;
    int residual_index = -1;
    bool required = true;
};

struct NativeRegressionRowDiagnostic {
    std::string row_id;
    bool success = false;
    std::string status = "invalid_input";
    std::string message;
    int residual_start = 0;
    int residual_count = 0;
    bool penalty_applied = false;
    std::string solve_backend;
    std::string derivative_backend;
};

struct NativeRegressionResidualRecord {
    std::string row_id;
    std::string row_kind;
    std::string family;
    std::string target;
    std::string name;
    double predicted = 0.0;
    double observed = 0.0;
    double scale = 1.0;
    std::map<std::string, double> sensitivities;
    bool success = true;
    bool recoverable_failure = false;
    std::string failure_message;
};

struct NativeRegressionResidualEvaluation {
    std::vector<double> residuals;
    std::vector<NativeRegressionResidualSchemaEntry> residual_schema;
    std::vector<NativeRegressionRowDiagnostic> row_diagnostics;
    double cost = 0.0;
    double residual_norm = 0.0;
    int success_count = 0;
    int failure_count = 0;
    bool fixed_shape_residuals = true;
};

struct NativeRegressionFitOptions {
    int max_iterations = 50;
    double gradient_tolerance = 1.0e-10;
    double function_tolerance = 1.0e-10;
    double parameter_tolerance = 1.0e-10;
    double penalty_residual = 1.0e6;
    std::string derivative_backend = "analytic";
    std::string optimizer_backend = "auto";
};

struct NativeRegressionFitResult {
    bool success = false;
    std::string status = "invalid_input";
    std::string message;
    std::string optimizer_backend;
    std::string derivative_backend;
    std::vector<double> parameters;
    std::vector<std::string> parameter_names;
    std::vector<double> lower_bounds;
    std::vector<double> upper_bounds;
    std::vector<bool> active_bounds;
    double initial_cost = 0.0;
    double final_cost = 0.0;
    double residual_norm = 0.0;
    double gradient_norm = 0.0;
    int iterations = 0;
    int function_evaluations = 0;
    NativeRegressionResidualEvaluation objective_result;
};

struct NativeRegressionProblemContract {
    std::vector<NativeRegressionParameterSpec> parameters;
    std::vector<NativeRegressionResidualSchemaEntry> residual_schema;
    std::vector<std::string> supported_target_families;
    std::vector<std::string> supported_parameter_kinds;
    bool fixed_shape_residuals = true;
    bool production_finite_difference_allowed = false;
};

std::string native_regression_status_name(NativeRegressionStatus status);
std::vector<std::string> native_regression_status_names();
NativeRegressionProblemContract native_regression_contract_schema();
NativeRegressionResidualEvaluation evaluate_native_regression_residual_records(
    const std::vector<NativeRegressionResidualRecord>& records,
    double penalty_residual = 1.0e6
);
NativeRegressionFitResult solve_native_regression_residual_records(
    const std::vector<NativeRegressionResidualRecord>& records,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    const NativeRegressionFitOptions& options
);
