#include "thermo_regression.h"

#include "epcsaft_electrolyte.h"
#include "implicit_sensitivity.h"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef EPCSAFT_HAS_CERES
#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>
#endif

namespace {

double scaled_residual(double predicted, double observed, double scale) {
    return (predicted - observed) * scale;
}

NativeRegressionResidualSchemaEntry schema_entry(
    const NativeThermoRegressionRow& row,
    const NativeThermoRegressionTarget& target,
    int index
) {
    NativeRegressionResidualSchemaEntry entry;
    entry.row_id = row.row_id;
    entry.family = target.family;
    entry.target = target.target;
    entry.scale = target.scale;
    entry.residual_index = index;
    entry.name = row.row_id + ":" + target.family + ":" + target.target;
    entry.required = true;
    return entry;
}

void append_penalty_row(
    NativeRegressionResidualEvaluation& out,
    const NativeThermoRegressionRow& row,
    const std::string& message,
    double penalty_residual
) {
    NativeRegressionRowDiagnostic diagnostic;
    diagnostic.row_id = row.row_id;
    diagnostic.success = false;
    diagnostic.status = native_regression_status_name(NativeRegressionStatus::LINE_SEARCH_FAILED);
    diagnostic.message = message;
    diagnostic.residual_start = static_cast<int>(out.residuals.size());
    diagnostic.residual_count = static_cast<int>(row.targets.size());
    diagnostic.penalty_applied = true;
    diagnostic.solve_backend = row.row_mode;
    diagnostic.derivative_backend = "backend_unavailable";
    for (const auto& target : row.targets) {
        const int index = static_cast<int>(out.residuals.size());
        out.residual_schema.push_back(schema_entry(row, target, index));
        out.residuals.push_back(penalty_residual);
    }
    out.row_diagnostics.push_back(diagnostic);
    out.failure_count += 1;
}

double speciation_prediction(const ChemicalEquilibriumResultNative& result, const NativeThermoRegressionTarget& target) {
    if (target.family == "speciation") {
        if (target.index < 0 || static_cast<std::size_t>(target.index) >= result.composition.size()) {
            throw std::runtime_error("speciation target index is out of range.");
        }
        return result.composition[static_cast<std::size_t>(target.index)];
    }
    if (target.family == "reaction") {
        if (target.index < 0 || static_cast<std::size_t>(target.index) >= result.reaction_residuals.size()) {
            throw std::runtime_error("reaction target index is out of range.");
        }
        return result.reaction_residuals[static_cast<std::size_t>(target.index)];
    }
    if (target.family == "activity") {
        if (target.index < 0 || static_cast<std::size_t>(target.index) >= result.activity_coefficients.size()) {
            throw std::runtime_error("activity target index is out of range.");
        }
        return result.activity_coefficients[static_cast<std::size_t>(target.index)];
    }
    throw std::runtime_error("unsupported reactive speciation target family: " + target.family);
}

double bubble_prediction(const EquilibriumResultNative& result, const NativeThermoRegressionTarget& target) {
    if (target.family == "pressure") {
        auto it = result.diagnostics_double.find("best_P");
        if (it != result.diagnostics_double.end()) {
            return it->second;
        }
        if (!result.phases.empty()) {
            return result.phases.front().pressure;
        }
        throw std::runtime_error("bubble pressure result did not include a pressure.");
    }
    if (target.family == "vapor_composition") {
        if (target.index < 0) {
            throw std::runtime_error("vapor composition target index is out of range.");
        }
        auto it = result.diagnostics_vector.find("best_y_vap");
        if (it == result.diagnostics_vector.end()
            || static_cast<std::size_t>(target.index) >= it->second.size()) {
            throw std::runtime_error("bubble result did not include requested vapor composition.");
        }
        return it->second[static_cast<std::size_t>(target.index)];
    }
    throw std::runtime_error("unsupported reactive electrolyte bubble target family: " + target.family);
}

int metadata_int(const NativeRegressionParameterSpec& parameter, const std::string& key, int fallback = -1) {
    auto found = parameter.metadata.find(key);
    if (found == parameter.metadata.end() || found->second.empty()) {
        return fallback;
    }
    return std::stoi(found->second);
}

bool parameter_matches_row(const NativeRegressionParameterSpec& parameter, const NativeThermoRegressionRow& row) {
    auto found = parameter.metadata.find("row_id");
    return found == parameter.metadata.end() || found->second.empty() || found->second == row.row_id;
}

void set_vector_parameter(std::vector<double>& values, int index, double value, const std::string& label) {
    if (index < 0 || static_cast<std::size_t>(index) >= values.size()) {
        throw std::runtime_error(label + " parameter index is out of range.");
    }
    values[static_cast<std::size_t>(index)] = value;
}

void set_symmetric_matrix_parameter(
    std::vector<double>& values,
    std::size_t n,
    int i,
    int j,
    double value,
    const std::string& label
) {
    if (i < 0 || j < 0 || static_cast<std::size_t>(i) >= n || static_cast<std::size_t>(j) >= n) {
        throw std::runtime_error(label + " parameter indices are out of range.");
    }
    if (values.size() != n * n) {
        throw std::runtime_error(label + " parameter requires a dense n-by-n matrix.");
    }
    values[static_cast<std::size_t>(i) * n + static_cast<std::size_t>(j)] = value;
    values[static_cast<std::size_t>(j) * n + static_cast<std::size_t>(i)] = value;
}

void apply_native_thermo_parameters(
    add_args& args,
    std::vector<NativeThermoRegressionRow>& rows,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    const std::vector<double>& theta
) {
    if (parameters.size() != theta.size()) {
        throw std::runtime_error("native thermodynamic regression parameter metadata must match theta length.");
    }
    const std::size_t n = args.m.size();
    for (std::size_t j = 0; j < parameters.size(); ++j) {
        const auto& parameter = parameters[j];
        const double value = theta[j];
        if (!std::isfinite(value)) {
            throw std::runtime_error("native thermodynamic regression parameter values must be finite.");
        }
        if (parameter.kind == "reaction_equilibrium_constant"
            || parameter.kind == "log_equilibrium_constant") {
            bool applied = false;
            const int reaction_index = metadata_int(parameter, "reaction_index", 0);
            for (auto& row : rows) {
                if (!parameter_matches_row(parameter, row)) {
                    continue;
                }
                set_vector_parameter(row.log_equilibrium_constants, reaction_index, value, "log_equilibrium_constant");
                applied = true;
            }
            if (!applied) {
                throw std::runtime_error("reaction_equilibrium_constant parameter did not match any regression row.");
            }
        } else if (parameter.kind == "born_radius" || parameter.kind == "born_diameter") {
            set_vector_parameter(args.d_born, metadata_int(parameter, "component_index"), value, parameter.kind);
        } else if (parameter.kind == "binary_interaction" || parameter.kind == "k_ij") {
            set_symmetric_matrix_parameter(
                args.k_ij,
                n,
                metadata_int(parameter, "component_index"),
                metadata_int(parameter, "other_component_index"),
                value,
                "k_ij"
            );
        } else if (parameter.kind == "segment_number" || parameter.kind == "m") {
            set_vector_parameter(args.m, metadata_int(parameter, "component_index"), value, "m");
        } else if (parameter.kind == "segment_diameter" || parameter.kind == "s") {
            set_vector_parameter(args.s, metadata_int(parameter, "component_index"), value, "s");
        } else if (parameter.kind == "dispersion_energy" || parameter.kind == "e") {
            set_vector_parameter(args.e, metadata_int(parameter, "component_index"), value, "e");
        } else {
            throw std::runtime_error("unsupported native thermodynamic regression parameter kind: " + parameter.kind);
        }
    }
}

NativeRegressionFitResult unavailable_native_thermo_fit_result(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<std::string>& species,
    const std::vector<NativeThermoRegressionRow>& rows,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    const NativeRegressionFitOptions& options,
    const std::string& message
) {
    NativeRegressionFitResult out;
    out.success = false;
    out.status = native_regression_status_name(NativeRegressionStatus::BACKEND_UNAVAILABLE);
    out.message = message;
    out.optimizer_backend = "backend_unavailable";
    out.derivative_backend = options.derivative_backend;
    out.objective_result = evaluate_native_thermo_regression_rows(mixture, species, rows, options.penalty_residual);
    out.initial_cost = out.objective_result.cost;
    out.final_cost = out.initial_cost;
    out.residual_norm = out.objective_result.residual_norm;
    for (const auto& parameter : parameters) {
        out.parameters.push_back(parameter.initial);
        out.parameter_names.push_back(parameter.name);
        out.lower_bounds.push_back(parameter.lower);
        out.upper_bounds.push_back(parameter.upper);
        out.active_bounds.push_back(false);
    }
    return out;
}

bool thermo_derivative_supported(
    const std::vector<NativeThermoRegressionRow>& rows,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    std::string* reason
) {
    for (const auto& parameter : parameters) {
        if (parameter.kind != "reaction_equilibrium_constant"
            && parameter.kind != "log_equilibrium_constant") {
            if (reason != nullptr) {
                *reason = "native Ceres thermodynamic derivatives currently support reaction log-equilibrium constants only.";
            }
            return false;
        }
    }
    for (const auto& row : rows) {
        if (row.row_mode != "reactive_speciation") {
            if (reason != nullptr) {
                *reason = "native Ceres thermodynamic derivatives currently support reactive_speciation rows only.";
            }
            return false;
        }
        for (int standard_state : row.reaction_standard_states) {
            if (standard_state != 1) {
                if (reason != nullptr) {
                    *reason = "implicit reactive-speciation derivatives require ideal mole-fraction reaction standard states.";
                }
                return false;
            }
        }
        for (const auto& target : row.targets) {
            if (target.family != "speciation") {
                if (reason != nullptr) {
                    *reason = "native Ceres thermodynamic derivatives currently support speciation targets only.";
                }
                return false;
            }
        }
    }
    return true;
}

std::vector<double> initial_theta(const std::vector<NativeRegressionParameterSpec>& parameters) {
    std::vector<double> theta;
    theta.reserve(parameters.size());
    for (const auto& parameter : parameters) {
        theta.push_back(std::max(parameter.lower, std::min(parameter.upper, parameter.initial)));
    }
    return theta;
}

std::shared_ptr<ePCSAFTMixtureNative> mixture_for_theta(
    const std::shared_ptr<ePCSAFTMixtureNative>& base_mixture,
    std::vector<NativeThermoRegressionRow>& rows,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    const std::vector<double>& theta
) {
    add_args args = base_mixture->args();
    apply_native_thermo_parameters(args, rows, parameters, theta);
    return std::make_shared<ePCSAFTMixtureNative>(args);
}

void fill_implicit_speciation_jacobian(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const NativeThermoRegressionRow& row,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    double* jacobian,
    int residual_offset,
    int parameter_count
) {
    std::vector<double> log_x;
    log_x.reserve(row.initial_x.size());
    for (double x : row.initial_x) {
        log_x.push_back(std::log(std::max(x, row.speciation_options.min_mole_fraction)));
    }
    ChemicalResidualEvaluationNative residual_eval = evaluate_chemical_equilibrium_residual_native(
        mixture,
        row.t,
        row.p,
        row.initial_x,
        log_x,
        true,
        row.balance_matrix_row_major,
        row.balance_rows,
        row.total_vector,
        row.reaction_stoichiometry_row_major,
        row.reaction_rows,
        row.log_equilibrium_constants,
        row.reaction_standard_states,
        row.speciation_options
    );
    std::vector<double> r_theta(static_cast<std::size_t>(residual_eval.jacobian_rows * parameter_count), 0.0);
    for (int p = 0; p < parameter_count; ++p) {
        const auto& parameter = parameters[static_cast<std::size_t>(p)];
        if (!parameter_matches_row(parameter, row)) {
            continue;
        }
        const int reaction_index = metadata_int(parameter, "reaction_index", 0);
        const int residual_row = row.balance_rows + 1 + reaction_index;
        if (residual_row < 0 || residual_row >= residual_eval.jacobian_rows) {
            throw std::runtime_error("reaction_equilibrium_constant derivative row is out of range.");
        }
        r_theta[static_cast<std::size_t>(residual_row * parameter_count + p)] = -1.0;
    }
    NativeImplicitSensitivityResult sensitivity = solve_native_implicit_sensitivity(
        residual_eval.jacobian_row_major,
        residual_eval.jacobian_rows,
        residual_eval.jacobian_cols,
        r_theta,
        parameter_count
    );
    if (!sensitivity.success) {
        throw std::runtime_error("implicit reactive-speciation sensitivity failed: " + sensitivity.message);
    }
    for (std::size_t target_index = 0; target_index < row.targets.size(); ++target_index) {
        const auto& target = row.targets[target_index];
        const int species_index = target.index;
        if (species_index < 0 || static_cast<std::size_t>(species_index) >= residual_eval.composition.size()) {
            throw std::runtime_error("speciation target index is out of range for sensitivity calculation.");
        }
        for (int p = 0; p < parameter_count; ++p) {
            double mean_log_sensitivity = 0.0;
            for (std::size_t i = 0; i < residual_eval.composition.size(); ++i) {
                mean_log_sensitivity += residual_eval.composition[i]
                    * sensitivity.sensitivities_row_major[i * static_cast<std::size_t>(parameter_count) + static_cast<std::size_t>(p)];
            }
            const double log_sensitivity =
                sensitivity.sensitivities_row_major[static_cast<std::size_t>(species_index) * static_cast<std::size_t>(parameter_count) + static_cast<std::size_t>(p)];
            const double dx_dtheta = residual_eval.composition[static_cast<std::size_t>(species_index)]
                * (log_sensitivity - mean_log_sensitivity);
            jacobian[(residual_offset + static_cast<int>(target_index)) * parameter_count + p] = target.scale * dx_dtheta;
        }
    }
}

#ifdef EPCSAFT_HAS_CERES
class NativeThermoCeresCostFunction final : public ceres::CostFunction {
public:
    NativeThermoCeresCostFunction(
        std::shared_ptr<ePCSAFTMixtureNative> mixture,
        std::vector<std::string> species,
        std::vector<NativeThermoRegressionRow> rows,
        std::vector<NativeRegressionParameterSpec> parameters,
        double penalty_residual,
        int residual_count
    )
        : mixture_(std::move(mixture)),
          species_(std::move(species)),
          rows_(std::move(rows)),
          parameters_(std::move(parameters)),
          penalty_residual_(penalty_residual) {
        set_num_residuals(residual_count);
        mutable_parameter_block_sizes()->push_back(static_cast<int>(parameters_.size()));
    }

    bool Evaluate(double const* const* theta_blocks, double* residuals, double** jacobians) const override {
        std::vector<double> theta(theta_blocks[0], theta_blocks[0] + parameters_.size());
        std::vector<NativeThermoRegressionRow> rows = rows_;
        auto mixture = mixture_for_theta(mixture_, rows, parameters_, theta);
        NativeRegressionResidualEvaluation eval =
            evaluate_native_thermo_regression_rows(mixture, species_, rows, penalty_residual_);
        if (eval.residuals.size() != static_cast<std::size_t>(num_residuals())) {
            return false;
        }
        for (int i = 0; i < num_residuals(); ++i) {
            residuals[i] = eval.residuals[static_cast<std::size_t>(i)];
        }
        if (jacobians != nullptr && jacobians[0] != nullptr) {
            const int parameter_count = static_cast<int>(parameters_.size());
            std::fill(jacobians[0], jacobians[0] + num_residuals() * parameter_count, 0.0);
            int residual_offset = 0;
            for (const auto& row : rows) {
                fill_implicit_speciation_jacobian(mixture, row, parameters_, jacobians[0], residual_offset, parameter_count);
                residual_offset += static_cast<int>(row.targets.size());
            }
        }
        return true;
    }

private:
    std::shared_ptr<ePCSAFTMixtureNative> mixture_;
    std::vector<std::string> species_;
    std::vector<NativeThermoRegressionRow> rows_;
    std::vector<NativeRegressionParameterSpec> parameters_;
    double penalty_residual_;
};
#endif

}  // namespace

NativeRegressionResidualEvaluation evaluate_native_thermo_regression_rows(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<std::string>& species,
    const std::vector<NativeThermoRegressionRow>& rows,
    double penalty_residual
) {
    NativeRegressionResidualEvaluation out;
    out.fixed_shape_residuals = true;
    for (const auto& row : rows) {
        NativeRegressionRowDiagnostic diagnostic;
        diagnostic.row_id = row.row_id;
        diagnostic.residual_start = static_cast<int>(out.residuals.size());
        diagnostic.residual_count = static_cast<int>(row.targets.size());
        diagnostic.penalty_applied = false;
        diagnostic.solve_backend = row.row_mode;
        diagnostic.derivative_backend = "backend_unavailable";
        try {
            if (row.row_mode == "reactive_speciation") {
                ChemicalEquilibriumResultNative result = chemical_equilibrium_native(
                    mixture,
                    row.t,
                    row.p,
                    row.initial_x,
                    row.balance_matrix_row_major,
                    row.balance_rows,
                    row.total_vector,
                    row.reaction_stoichiometry_row_major,
                    row.reaction_rows,
                    row.log_equilibrium_constants,
                    row.reaction_standard_states,
                    row.speciation_options
                );
                diagnostic.success = result.success;
                diagnostic.status = result.success
                    ? native_regression_status_name(NativeRegressionStatus::CONVERGED)
                    : native_regression_status_name(NativeRegressionStatus::LINE_SEARCH_FAILED);
                diagnostic.message = result.message;
                diagnostic.solve_backend = "native_chemical_equilibrium";
                auto backend = result.diagnostics_string.find("jacobian_backend");
                diagnostic.derivative_backend = backend == result.diagnostics_string.end()
                    ? "backend_unavailable"
                    : backend->second;
                if (!result.success) {
                    append_penalty_row(out, row, result.message, penalty_residual);
                    continue;
                }
                for (const auto& target : row.targets) {
                    const double predicted = speciation_prediction(result, target);
                    const int index = static_cast<int>(out.residuals.size());
                    out.residual_schema.push_back(schema_entry(row, target, index));
                    out.residuals.push_back(scaled_residual(predicted, target.observed, target.scale));
                }
                out.row_diagnostics.push_back(diagnostic);
                out.success_count += 1;
            } else if (row.row_mode == "reactive_electrolyte_bubble") {
                const std::vector<std::string>& vapor_species =
                    row.vapor_species.empty() ? row.targets_species : row.vapor_species;
                EquilibriumResultNative result = electrolyte_bubble_pressure_native(
                    mixture,
                    row.t,
                    row.x_liq.empty() ? row.initial_x : row.x_liq,
                    row.bubble_options,
                    species,
                    vapor_species
                );
                const bool success = result.diagnostics_bool.count("success")
                    ? result.diagnostics_bool.at("success")
                    : result.split_detected;
                diagnostic.success = success;
                diagnostic.status = success
                    ? native_regression_status_name(NativeRegressionStatus::CONVERGED)
                    : native_regression_status_name(NativeRegressionStatus::LINE_SEARCH_FAILED);
                diagnostic.message = result.diagnostics_string.count("message")
                    ? result.diagnostics_string.at("message")
                    : "";
                diagnostic.solve_backend = "native_electrolyte_bubble";
                diagnostic.derivative_backend = "not_differentiated";
                if (!success) {
                    append_penalty_row(out, row, diagnostic.message, penalty_residual);
                    continue;
                }
                for (const auto& target : row.targets) {
                    const double predicted = bubble_prediction(result, target);
                    const int index = static_cast<int>(out.residuals.size());
                    out.residual_schema.push_back(schema_entry(row, target, index));
                    out.residuals.push_back(scaled_residual(predicted, target.observed, target.scale));
                }
                out.row_diagnostics.push_back(diagnostic);
                out.success_count += 1;
            } else {
                append_penalty_row(out, row, "unsupported native thermodynamic regression row mode", penalty_residual);
            }
        } catch (const std::exception& exc) {
            append_penalty_row(out, row, exc.what(), penalty_residual);
        }
    }
    double cost = 0.0;
    for (double residual : out.residuals) {
        cost += 0.5 * residual * residual;
    }
    out.cost = cost;
    out.residual_norm = std::sqrt(2.0 * cost);
    return out;
}

NativeRegressionFitResult fit_native_thermo_regression(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const std::vector<std::string>& species,
    const std::vector<NativeThermoRegressionRow>& rows,
    const std::vector<NativeRegressionParameterSpec>& parameters,
    const NativeRegressionFitOptions& options
) {
    NativeRegressionFitResult out;
    out.optimizer_backend = "ceres";
    out.derivative_backend = "analytic_implicit";
    for (const auto& parameter : parameters) {
        out.parameter_names.push_back(parameter.name);
        out.lower_bounds.push_back(parameter.lower);
        out.upper_bounds.push_back(parameter.upper);
    }
    if (parameters.empty()) {
        out.status = native_regression_status_name(NativeRegressionStatus::INVALID_INPUT);
        out.message = "native thermodynamic regression requires at least one parameter.";
        return out;
    }
    if (rows.empty()) {
        out.status = native_regression_status_name(NativeRegressionStatus::ALL_ROWS_FAILED);
        out.message = "native thermodynamic regression requires at least one row.";
        return out;
    }
    for (const auto& parameter : parameters) {
        if (!std::isfinite(parameter.initial) || !std::isfinite(parameter.lower) || !std::isfinite(parameter.upper)
            || parameter.upper < parameter.lower) {
            out.status = native_regression_status_name(NativeRegressionStatus::BOUNDS_INCONSISTENT);
            out.message = "native thermodynamic regression parameter bounds must be finite and ordered.";
            return out;
        }
    }
    std::string unsupported_reason;
    if (!thermo_derivative_supported(rows, parameters, &unsupported_reason)) {
        return unavailable_native_thermo_fit_result(
            mixture,
            species,
            rows,
            parameters,
            options,
            unsupported_reason
        );
    }
    std::vector<double> theta = initial_theta(parameters);
    std::vector<NativeThermoRegressionRow> initial_rows = rows;
    auto initial_mixture = mixture_for_theta(mixture, initial_rows, parameters, theta);
    out.objective_result = evaluate_native_thermo_regression_rows(
        initial_mixture,
        species,
        initial_rows,
        options.penalty_residual
    );
    out.initial_cost = out.objective_result.cost;
    out.final_cost = out.initial_cost;
    out.residual_norm = out.objective_result.residual_norm;
    if (out.objective_result.residuals.empty()) {
        out.status = native_regression_status_name(NativeRegressionStatus::ALL_ROWS_FAILED);
        out.message = "native thermodynamic regression generated no residuals.";
        return out;
    }
    if (out.objective_result.failure_count == static_cast<int>(rows.size())) {
        out.status = native_regression_status_name(NativeRegressionStatus::ALL_ROWS_FAILED);
        out.message = "all native thermodynamic regression rows failed before optimization.";
        return out;
    }

#ifndef EPCSAFT_HAS_CERES
    return unavailable_native_thermo_fit_result(
        mixture,
        species,
        rows,
        parameters,
        options,
        "native thermodynamic Ceres backend is not compiled; rebuild with EPCSAFT_ENABLE_CERES=ON."
    );
#else
    ceres::Problem problem;
    auto* cost = new NativeThermoCeresCostFunction(
        mixture,
        species,
        rows,
        parameters,
        options.penalty_residual,
        static_cast<int>(out.objective_result.residuals.size())
    );
    problem.AddResidualBlock(cost, nullptr, theta.data());
    for (std::size_t j = 0; j < parameters.size(); ++j) {
        problem.SetParameterLowerBound(theta.data(), static_cast<int>(j), parameters[j].lower);
        problem.SetParameterUpperBound(theta.data(), static_cast<int>(j), parameters[j].upper);
    }

    ceres::Solver::Options solver_options;
    solver_options.max_num_iterations = options.max_iterations;
    solver_options.function_tolerance = options.function_tolerance;
    solver_options.gradient_tolerance = options.gradient_tolerance;
    solver_options.parameter_tolerance = options.parameter_tolerance;
    solver_options.minimizer_progress_to_stdout = false;
    solver_options.logging_type = ceres::SILENT;
    ceres::Solver::Summary summary;
    ceres::Solve(solver_options, &problem, &summary);

    std::vector<NativeThermoRegressionRow> final_rows = rows;
    auto final_mixture = mixture_for_theta(mixture, final_rows, parameters, theta);
    NativeRegressionResidualEvaluation final_eval = evaluate_native_thermo_regression_rows(
        final_mixture,
        species,
        final_rows,
        options.penalty_residual
    );
    out.parameters = theta;
    out.active_bounds.reserve(theta.size());
    for (std::size_t j = 0; j < theta.size(); ++j) {
        out.active_bounds.push_back(theta[j] == parameters[j].lower || theta[j] == parameters[j].upper);
    }
    out.objective_result = final_eval;
    out.final_cost = final_eval.cost;
    out.residual_norm = final_eval.residual_norm;
    out.gradient_norm = 0.0;
    out.iterations = static_cast<int>(summary.iterations.size());
    out.function_evaluations = summary.num_residual_evaluations + summary.num_jacobian_evaluations;
    out.success = summary.IsSolutionUsable() && std::isfinite(out.final_cost);
    out.status = out.success
        ? native_regression_status_name(NativeRegressionStatus::CONVERGED)
        : native_regression_status_name(NativeRegressionStatus::LINE_SEARCH_FAILED);
    if (summary.termination_type == ceres::NO_CONVERGENCE) {
        out.status = native_regression_status_name(NativeRegressionStatus::MAX_ITERATIONS);
    }
    out.message = summary.BriefReport()
        + "; native_hot_loop=true; python_objective_used=false; finite_difference_used=false";
    return out;
#endif
}
