#include "regression_types.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef EPCSAFT_HAS_CERES
#include <ceres/ceres.h>
#endif

namespace {

std::vector<double> solve_dense_linear_system(std::vector<std::vector<double>> matrix, std::vector<double> rhs) {
    const std::size_t n = rhs.size();
    for (std::size_t pivot = 0; pivot < n; ++pivot) {
        std::size_t best = pivot;
        double best_abs = std::abs(matrix[pivot][pivot]);
        for (std::size_t row = pivot + 1; row < n; ++row) {
            const double value = std::abs(matrix[row][pivot]);
            if (value > best_abs) {
                best = row;
                best_abs = value;
            }
        }
        if (best_abs <= 1.0e-18 || !std::isfinite(best_abs)) {
            throw std::runtime_error("singular normal equation matrix");
        }
        if (best != pivot) {
            std::swap(matrix[pivot], matrix[best]);
            std::swap(rhs[pivot], rhs[best]);
        }
        const double diagonal = matrix[pivot][pivot];
        for (std::size_t col = pivot; col < n; ++col) {
            matrix[pivot][col] /= diagonal;
        }
        rhs[pivot] /= diagonal;
        for (std::size_t row = 0; row < n; ++row) {
            if (row == pivot) {
                continue;
            }
            const double factor = matrix[row][pivot];
            for (std::size_t col = pivot; col < n; ++col) {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }
    return rhs;
}

double gradient_norm_from_records(
    const std::vector<NativeRegressionResidualRecord>& records,
    const std::vector<std::string>& parameter_names,
    double penalty_residual
) {
    std::vector<double> gradient(parameter_names.size(), 0.0);
    for (const auto& record : records) {
        const bool finite_values = std::isfinite(record.predicted) && std::isfinite(record.observed);
        const bool success = record.success && finite_values && !record.recoverable_failure;
        const double residual = success ? record.scale * (record.predicted - record.observed) : penalty_residual;
        for (std::size_t j = 0; j < parameter_names.size(); ++j) {
            const auto found = record.sensitivities.find(parameter_names[j]);
            if (found != record.sensitivities.end()) {
                gradient[j] += record.scale * found->second * residual;
            }
        }
    }
    double norm2 = 0.0;
    for (const double value : gradient) {
        norm2 += value * value;
    }
    return std::sqrt(norm2);
}

} // namespace

std::string native_regression_status_name(NativeRegressionStatus status) {
    switch (status) {
        case NativeRegressionStatus::CONVERGED:
            return "converged";
        case NativeRegressionStatus::MAX_ITERATIONS:
            return "max_iterations";
        case NativeRegressionStatus::LINE_SEARCH_FAILED:
            return "line_search_failed";
        case NativeRegressionStatus::SINGULAR_JACOBIAN:
            return "singular_jacobian";
        case NativeRegressionStatus::ALL_ROWS_FAILED:
            return "all_rows_failed";
        case NativeRegressionStatus::NONFINITE_OBJECTIVE:
            return "nonfinite_objective";
        case NativeRegressionStatus::BOUNDS_INCONSISTENT:
            return "bounds_inconsistent";
        case NativeRegressionStatus::INVALID_INPUT:
            return "invalid_input";
    }
    throw std::invalid_argument("Unknown native regression status.");
}

std::vector<std::string> native_regression_status_names() {
    return {
        native_regression_status_name(NativeRegressionStatus::CONVERGED),
        native_regression_status_name(NativeRegressionStatus::MAX_ITERATIONS),
        native_regression_status_name(NativeRegressionStatus::LINE_SEARCH_FAILED),
        native_regression_status_name(NativeRegressionStatus::SINGULAR_JACOBIAN),
        native_regression_status_name(NativeRegressionStatus::ALL_ROWS_FAILED),
        native_regression_status_name(NativeRegressionStatus::NONFINITE_OBJECTIVE),
        native_regression_status_name(NativeRegressionStatus::BOUNDS_INCONSISTENT),
        native_regression_status_name(NativeRegressionStatus::INVALID_INPUT),
    };
}

NativeRegressionProblemContract native_regression_contract_schema() {
    NativeRegressionProblemContract contract;
    contract.supported_target_families = {
        "pressure",
        "speciation",
        "activity",
        "fugacity",
        "density",
        "relative_permittivity",
    };
    contract.supported_parameter_kinds = {
        "pure_component",
        "binary_interaction",
        "reaction_equilibrium_constant",
        "born_radius",
        "dielectric_parameter",
    };
    contract.fixed_shape_residuals = true;
    contract.production_finite_difference_allowed = false;
    return contract;
}

NativeRegressionResidualEvaluation evaluate_native_regression_residual_records(
    const std::vector<NativeRegressionResidualRecord>& records,
    double penalty_residual
) {
    if (!(penalty_residual > 0.0) || !std::isfinite(penalty_residual)) {
        throw std::invalid_argument("Native regression penalty_residual must be finite and positive.");
    }

    NativeRegressionResidualEvaluation out;
    out.residuals.reserve(records.size());
    out.residual_schema.reserve(records.size());
    std::map<std::string, NativeRegressionRowDiagnostic> diagnostics_by_row;

    for (std::size_t i = 0; i < records.size(); ++i) {
        const NativeRegressionResidualRecord& record = records[i];
        if (!(record.scale > 0.0) || !std::isfinite(record.scale)) {
            throw std::invalid_argument("Native regression residual scale must be finite and positive.");
        }
        if (record.row_id.empty() || record.family.empty() || record.target.empty()) {
            throw std::invalid_argument("Native regression residual records require row_id, family, and target.");
        }

        NativeRegressionResidualSchemaEntry schema;
        schema.name = record.name.empty()
            ? record.row_id + ":" + record.family + ":" + record.target
            : record.name;
        schema.row_id = record.row_id;
        schema.family = record.family;
        schema.target = record.target;
        schema.scale = record.scale;
        schema.residual_index = static_cast<int>(i);
        schema.required = true;
        out.residual_schema.push_back(schema);

        double residual = penalty_residual;
        const bool finite_values = std::isfinite(record.predicted) && std::isfinite(record.observed);
        const bool success = record.success && finite_values && !record.recoverable_failure;
        if (success) {
            residual = record.scale * (record.predicted - record.observed);
            ++out.success_count;
        } else {
            ++out.failure_count;
        }
        out.residuals.push_back(residual);
        out.cost += 0.5 * residual * residual;

        NativeRegressionRowDiagnostic& diagnostic = diagnostics_by_row[record.row_id];
        if (diagnostic.row_id.empty()) {
            diagnostic.row_id = record.row_id;
            diagnostic.success = true;
            diagnostic.status = native_regression_status_name(NativeRegressionStatus::CONVERGED);
            diagnostic.residual_start = static_cast<int>(i);
            diagnostic.residual_count = 0;
        }
        diagnostic.residual_count += 1;
        diagnostic.success = diagnostic.success && success;
        if (!success) {
            diagnostic.status = finite_values
                ? native_regression_status_name(NativeRegressionStatus::LINE_SEARCH_FAILED)
                : native_regression_status_name(NativeRegressionStatus::NONFINITE_OBJECTIVE);
            diagnostic.penalty_applied = true;
            diagnostic.message = record.failure_message.empty() ? "recoverable row failure" : record.failure_message;
        }
    }

    out.residual_norm = std::sqrt(std::max(0.0, 2.0 * out.cost));
    out.row_diagnostics.reserve(diagnostics_by_row.size());
    for (auto& item : diagnostics_by_row) {
        out.row_diagnostics.push_back(item.second);
    }
    out.fixed_shape_residuals = true;
    return out;
}

NativeRegressionFitResult solve_native_regression_residual_records(
    const std::vector<NativeRegressionResidualRecord>& records,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    const NativeRegressionFitOptions& options
) {
    NativeRegressionFitResult out;
    out.optimizer_backend = options.optimizer_backend == "auto" ? "analytic_linear_native" : options.optimizer_backend;
    out.derivative_backend = options.derivative_backend;
    out.objective_result = evaluate_native_regression_residual_records(records, options.penalty_residual);
    out.initial_cost = out.objective_result.cost;
    out.final_cost = out.initial_cost;
    out.residual_norm = out.objective_result.residual_norm;

    if (options.derivative_backend != "analytic"
        && options.derivative_backend != "cppad"
        && options.derivative_backend != "implicit") {
        out.status = native_regression_status_name(NativeRegressionStatus::INVALID_INPUT);
        out.message = "native production regression derivatives must be analytic, implicit, or CppAD; finite_difference is debug-only.";
        return out;
    }
    if (parameters.empty()) {
        out.status = native_regression_status_name(NativeRegressionStatus::INVALID_INPUT);
        out.message = "native regression requires at least one parameter.";
        return out;
    }
    if (records.empty()) {
        out.status = native_regression_status_name(NativeRegressionStatus::ALL_ROWS_FAILED);
        out.message = "native regression requires at least one residual record.";
        return out;
    }
    if (out.objective_result.failure_count == static_cast<int>(records.size())) {
        out.status = native_regression_status_name(NativeRegressionStatus::ALL_ROWS_FAILED);
        out.message = "all native regression rows failed before optimization.";
        return out;
    }

    out.parameters.reserve(parameters.size());
    out.parameter_names.reserve(parameters.size());
    out.lower_bounds.reserve(parameters.size());
    out.upper_bounds.reserve(parameters.size());
    out.active_bounds.reserve(parameters.size());
    for (const auto& parameter : parameters) {
        if (!(parameter.scale > 0.0) || !std::isfinite(parameter.scale)) {
            out.status = native_regression_status_name(NativeRegressionStatus::INVALID_INPUT);
            out.message = "native regression parameter scale must be finite and positive.";
            return out;
        }
        if (!std::isfinite(parameter.lower) || !std::isfinite(parameter.upper) || parameter.upper < parameter.lower) {
            out.status = native_regression_status_name(NativeRegressionStatus::BOUNDS_INCONSISTENT);
            out.message = "native regression parameter bounds must be finite and ordered.";
            return out;
        }
        double value = parameter.initial;
        if (!std::isfinite(value)) {
            out.status = native_regression_status_name(NativeRegressionStatus::INVALID_INPUT);
            out.message = "native regression initial parameter values must be finite.";
            return out;
        }
        value = std::max(parameter.lower, std::min(parameter.upper, value));
        out.parameters.push_back(value);
        out.parameter_names.push_back(parameter.name);
        out.lower_bounds.push_back(parameter.lower);
        out.upper_bounds.push_back(parameter.upper);
        out.active_bounds.push_back(value == parameter.lower || value == parameter.upper);
    }

    out.gradient_norm = gradient_norm_from_records(records, out.parameter_names, options.penalty_residual);
    bool has_sensitivity = false;
    for (const auto& record : records) {
        for (const auto& parameter_name : out.parameter_names) {
            const auto found = record.sensitivities.find(parameter_name);
            if (found != record.sensitivities.end() && std::isfinite(found->second) && found->second != 0.0) {
                has_sensitivity = true;
                break;
            }
        }
        if (has_sensitivity) {
            break;
        }
    }

    if (has_sensitivity) {
        const std::size_t n = out.parameter_names.size();
        std::vector<std::vector<double>> jtj(n, std::vector<double>(n, 0.0));
        std::vector<double> jtr(n, 0.0);
        for (const auto& record : records) {
            const bool finite_values = std::isfinite(record.predicted) && std::isfinite(record.observed);
            const bool success = record.success && finite_values && !record.recoverable_failure;
            if (!success) {
                continue;
            }
            const double residual = record.scale * (record.predicted - record.observed);
            std::vector<double> jacobian_row(n, 0.0);
            for (std::size_t j = 0; j < n; ++j) {
                const auto found = record.sensitivities.find(out.parameter_names[j]);
                if (found != record.sensitivities.end()) {
                    if (!std::isfinite(found->second)) {
                        out.status = native_regression_status_name(NativeRegressionStatus::INVALID_INPUT);
                        out.message = "native regression sensitivities must be finite.";
                        return out;
                    }
                    jacobian_row[j] = record.scale * found->second;
                }
            }
            for (std::size_t j = 0; j < n; ++j) {
                jtr[j] += jacobian_row[j] * residual;
                for (std::size_t k = 0; k < n; ++k) {
                    jtj[j][k] += jacobian_row[j] * jacobian_row[k];
                }
            }
        }
        for (std::size_t j = 0; j < n; ++j) {
            jtj[j][j] += 1.0e-12;
            jtr[j] = -jtr[j];
        }
        std::vector<double> step;
        try {
            step = solve_dense_linear_system(jtj, jtr);
        } catch (const std::runtime_error&) {
            out.status = native_regression_status_name(NativeRegressionStatus::SINGULAR_JACOBIAN);
            out.message = "native regression analytic sensitivity matrix is singular.";
            return out;
        }
        std::vector<double> delta(n, 0.0);
        for (std::size_t j = 0; j < n; ++j) {
            const double trial = std::max(out.lower_bounds[j], std::min(out.upper_bounds[j], out.parameters[j] + step[j]));
            delta[j] = trial - out.parameters[j];
            out.parameters[j] = trial;
            out.active_bounds[j] = trial == out.lower_bounds[j] || trial == out.upper_bounds[j];
        }
        std::vector<NativeRegressionResidualRecord> adjusted_records = records;
        for (auto& record : adjusted_records) {
            for (std::size_t j = 0; j < n; ++j) {
                const auto found = record.sensitivities.find(out.parameter_names[j]);
                if (found != record.sensitivities.end()) {
                    record.predicted += found->second * delta[j];
                }
            }
        }
        const NativeRegressionResidualEvaluation adjusted =
            evaluate_native_regression_residual_records(adjusted_records, options.penalty_residual);
        if (!std::isfinite(adjusted.cost)) {
            out.status = native_regression_status_name(NativeRegressionStatus::NONFINITE_OBJECTIVE);
            out.message = "native regression produced a nonfinite objective.";
            return out;
        }
        if (adjusted.cost > out.initial_cost + options.function_tolerance) {
            out.status = native_regression_status_name(NativeRegressionStatus::LINE_SEARCH_FAILED);
            out.message = "native regression analytic sensitivity step did not improve the objective.";
            return out;
        }
        out.objective_result = adjusted;
        out.final_cost = adjusted.cost;
        out.residual_norm = adjusted.residual_norm;
        out.gradient_norm = gradient_norm_from_records(adjusted_records, out.parameter_names, options.penalty_residual);
        out.iterations = 1;
        out.function_evaluations = 2;
    } else {
        out.iterations = 0;
        out.function_evaluations = 1;
    }

    out.success = true;
    out.status = native_regression_status_name(NativeRegressionStatus::CONVERGED);
    if (has_sensitivity) {
        out.message = "native analytic sensitivity regression converged without production finite differences.";
    } else {
        out.message = "analytic fixed-shape residual contract evaluated without production finite differences.";
    }
    return out;
}
