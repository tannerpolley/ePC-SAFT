#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "epcsaft_chemical_equilibrium.h"
#include "epcsaft_equilibrium.h"
#include "ad_derivative_checks.h"
#include "implicit_sensitivity.h"
#include "regression_types.h"
#include "thermo_regression.h"

namespace py = pybind11;

namespace {

py::object native_solution_error_type;

ChemicalEquilibriumOptionsNative chemical_options_from_request(const py::dict& request);
ElectrolyteBubbleOptionsNative electrolyte_bubble_options_from_request(const py::dict& request);

std::vector<double> array_to_double_vector(const py::array& array) {
    py::array_t<double, py::array::forcecast> casted(array);
    py::buffer_info info = casted.request();
    const auto* data = static_cast<const double*>(info.ptr);
    std::size_t size = 1;
    for (py::ssize_t dim : info.shape) {
        size *= static_cast<std::size_t>(dim);
    }
    return std::vector<double>(data, data + size);
}

std::vector<int> array_to_int_vector(const py::array& array) {
    py::array_t<int, py::array::forcecast> casted(array);
    py::buffer_info info = casted.request();
    const auto* data = static_cast<const int*>(info.ptr);
    std::size_t size = 1;
    for (py::ssize_t dim : info.shape) {
        size *= static_cast<std::size_t>(dim);
    }
    return std::vector<int>(data, data + size);
}

std::vector<PureNeutralRegressionDensityRecord> density_records_from_arrays(
    const py::array& density_t,
    const py::array& density_p,
    const py::array& density_rho_exp,
    const py::array& density_phase
) {
    auto t = array_to_double_vector(density_t);
    auto p = array_to_double_vector(density_p);
    auto rho = array_to_double_vector(density_rho_exp);
    auto phase = array_to_int_vector(density_phase);
    if (t.size() != p.size() || t.size() != rho.size() || t.size() != phase.size()) {
        throw std::invalid_argument("density record arrays must have matching lengths");
    }
    std::vector<PureNeutralRegressionDensityRecord> records;
    records.reserve(t.size());
    for (std::size_t i = 0; i < t.size(); ++i) {
        PureNeutralRegressionDensityRecord record;
        record.t = t[i];
        record.p = p[i];
        record.rho_exp = rho[i];
        record.phase = phase[i];
        records.push_back(record);
    }
    return records;
}

std::vector<PureNeutralRegressionVLERecord> vle_records_from_arrays(
    const py::array& vle_t,
    const py::array& vle_p
) {
    auto t = array_to_double_vector(vle_t);
    auto p = array_to_double_vector(vle_p);
    if (t.size() != p.size()) {
        throw std::invalid_argument("pure VLE record arrays must have matching lengths");
    }
    std::vector<PureNeutralRegressionVLERecord> records;
    records.reserve(t.size());
    for (std::size_t i = 0; i < t.size(); ++i) {
        PureNeutralRegressionVLERecord record;
        record.t = t[i];
        record.p = p[i];
        records.push_back(record);
    }
    return records;
}

py::dict regression_result_to_dict(const PureNeutralRegressionResult& result) {
    py::dict out;
    out["x"] = result.x;
    out["cost"] = result.cost;
    out["residual_norm"] = result.residual_norm;
    out["density_metric"] = result.density_metric;
    out["pure_vle_metric"] = result.pure_vle_metric;
    out["initial_cost"] = result.initial_cost;
    out["initial_density_metric"] = result.initial_density_metric;
    out["initial_pure_vle_metric"] = result.initial_pure_vle_metric;
    out["success"] = result.success;
    out["status"] = result.status;
    out["nfev"] = result.nfev;
    out["iterations"] = result.iterations;
    out["starts_tried"] = result.starts_tried;
    out["objective_evaluations"] = result.objective_evaluations;
    out["gradient_evaluations"] = result.gradient_evaluations;
    out["residual_evaluations"] = result.residual_evaluations;
    out["density_solves"] = result.density_solves;
    out["fused_state_evaluations"] = result.fused_state_evaluations;
    out["callback_wall_time_s"] = result.callback_wall_time_s;
    out["solve_wall_time_s"] = result.solve_wall_time_s;
    out["message"] = result.message;
    out["backend"] = result.backend;
    out["jacobian_available"] = result.jacobian_available;
    out["jacobian_backend"] = result.jacobian_backend;
    out["jacobian_fallback_used"] = result.jacobian_fallback_used;
    out["jacobian_fallback_reason"] = result.jacobian_fallback_reason;
    out["hessian_available"] = result.hessian_available;
    out["hessian_backend"] = result.hessian_backend;
    out["hessian_fallback_used"] = result.hessian_fallback_used;
    out["hessian_fallback_reason"] = result.hessian_fallback_reason;
    return out;
}

py::dict regression_debug_to_dict(const PureNeutralRegressionDebugResult& result) {
    py::dict out;
    out["objective"] = result.objective;
    out["gradient"] = result.gradient;
    out["residuals"] = result.residuals;
    out["jacobian_row_major"] = result.jacobian_row_major;
    out["jacobian_shape"] = py::make_tuple(result.jacobian_rows, result.jacobian_cols);
    out["jacobian_available"] = result.jacobian_available;
    out["jacobian_backend"] = result.jacobian_backend;
    out["jacobian_fallback_used"] = result.jacobian_fallback_used;
    out["jacobian_fallback_reason"] = result.jacobian_fallback_reason;
    out["hessian_row_major"] = result.hessian_row_major;
    out["hessian_shape"] = py::make_tuple(result.hessian_rows, result.hessian_cols);
    out["hessian_available"] = result.hessian_available;
    out["hessian_backend"] = result.hessian_backend;
    out["hessian_fallback_used"] = result.hessian_fallback_used;
    out["hessian_fallback_reason"] = result.hessian_fallback_reason;
    out["density_raw_residuals"] = result.density_raw_residuals;
    out["pure_vle_raw_residuals"] = result.pure_vle_raw_residuals;
    out["residual_evaluations"] = result.residual_evaluations;
    out["density_solves"] = result.density_solves;
    out["fused_state_evaluations"] = result.fused_state_evaluations;
    out["callback_wall_time_s"] = result.callback_wall_time_s;
    return out;
}

py::dict generic_regression_result_to_dict(const GenericRegressionResult& result) {
    py::dict out;
    out["x"] = result.x;
    out["cost"] = result.cost;
    out["residual_norm"] = result.residual_norm;
    out["initial_cost"] = result.initial_cost;
    out["initial_residual_norm"] = result.initial_residual_norm;
    py::dict metrics;
    for (const auto& item : result.metrics_by_term) {
        metrics[py::str(item.first)] = item.second;
    }
    out["metrics_by_term"] = metrics;
    out["success"] = result.success;
    out["status"] = result.status;
    out["nfev"] = result.nfev;
    out["iterations"] = result.iterations;
    out["starts_tried"] = result.starts_tried;
    out["message"] = result.message;
    out["backend"] = result.backend;
    out["jacobian_available"] = result.jacobian_available;
    out["jacobian_backend"] = result.jacobian_backend;
    out["jacobian_fallback_used"] = result.jacobian_fallback_used;
    out["jacobian_fallback_reason"] = result.jacobian_fallback_reason;
    out["unsupported_derivative_fallback_count"] = result.unsupported_derivative_fallback_count;
    out["hessian_available"] = result.hessian_available;
    out["hessian_backend"] = result.hessian_backend;
    out["hessian_fallback_used"] = result.hessian_fallback_used;
    out["hessian_fallback_reason"] = result.hessian_fallback_reason;
    return out;
}

py::dict generic_regression_debug_to_dict(const GenericRegressionDebugResult& result) {
    py::dict out;
    out["cost"] = result.cost;
    out["residual_norm"] = result.residual_norm;
    out["residuals"] = result.residuals;
    py::dict metrics;
    for (const auto& item : result.metrics_by_term) {
        metrics[py::str(item.first)] = item.second;
    }
    out["metrics_by_term"] = metrics;
    out["jacobian_row_major"] = result.jacobian_row_major;
    out["jacobian_shape"] = py::make_tuple(result.jacobian_rows, result.jacobian_cols);
    out["jacobian_available"] = result.jacobian_available;
    out["jacobian_backend"] = result.jacobian_backend;
    out["jacobian_fallback_used"] = result.jacobian_fallback_used;
    out["jacobian_fallback_reason"] = result.jacobian_fallback_reason;
    out["unsupported_derivative_fallback_count"] = result.unsupported_derivative_fallback_count;
    out["hessian_row_major"] = result.hessian_row_major;
    out["hessian_shape"] = py::make_tuple(result.hessian_rows, result.hessian_cols);
    out["hessian_available"] = result.hessian_available;
    out["hessian_backend"] = result.hessian_backend;
    out["hessian_fallback_used"] = result.hessian_fallback_used;
    out["hessian_fallback_reason"] = result.hessian_fallback_reason;
    return out;
}

py::dict native_regression_parameter_spec_to_dict(const NativeRegressionParameterSpec& spec) {
    py::dict out;
    out["name"] = spec.name;
    out["path"] = spec.path;
    out["kind"] = spec.kind;
    out["initial"] = spec.initial;
    out["lower"] = spec.lower;
    out["upper"] = spec.upper;
    out["scale"] = spec.scale;
    out["fixed"] = spec.fixed;
    py::dict metadata;
    for (const auto& item : spec.metadata) {
        metadata[py::str(item.first)] = item.second;
    }
    out["metadata"] = metadata;
    return out;
}

py::dict native_autodiff_derivative_check_to_dict(
    const epcsaft::autodiff::NativeAutodiffDerivativeCheckResult& result
) {
    py::dict out;
    out["cppad_compiled"] = result.cppad_compiled;
    out["cppad_used"] = result.cppad_used;
    out["unsupported_derivative_used"] = result.unsupported_derivative_used;
    out["status"] = result.status;
    out["derivative_backend"] = result.derivative_backend;
    out["checked_residuals"] = result.checked_residuals;
    py::dict derivatives;
    for (const auto& item : result.derivative_by_residual) {
        derivatives[py::str(item.first)] = item.second;
    }
    out["derivative_by_residual"] = derivatives;
    out["max_abs_error"] = result.max_abs_error;
    return out;
}

py::dict native_implicit_sensitivity_to_dict(const NativeImplicitSensitivityResult& result) {
    py::dict out;
    out["success"] = result.success;
    out["status"] = result.status;
    out["message"] = result.message;
    out["sensitivities_row_major"] = result.sensitivities_row_major;
    out["shape"] = py::make_tuple(result.rows, result.cols);
    out["residual_jacobian_condition_proxy"] = result.residual_jacobian_condition_proxy;
    out["unsupported_derivative_used"] = false;
    out["sensitivity_backend"] = "analytic_implicit";
    return out;
}

py::dict native_regression_residual_schema_entry_to_dict(const NativeRegressionResidualSchemaEntry& entry) {
    py::dict out;
    out["name"] = entry.name;
    out["row_id"] = entry.row_id;
    out["family"] = entry.family;
    out["target"] = entry.target;
    out["scale"] = entry.scale;
    out["residual_index"] = entry.residual_index;
    out["required"] = entry.required;
    return out;
}

py::dict native_regression_row_diagnostic_to_dict(const NativeRegressionRowDiagnostic& diagnostic) {
    py::dict out;
    out["row_id"] = diagnostic.row_id;
    out["success"] = diagnostic.success;
    out["status"] = diagnostic.status;
    out["message"] = diagnostic.message;
    out["residual_start"] = diagnostic.residual_start;
    out["residual_count"] = diagnostic.residual_count;
    out["penalty_applied"] = diagnostic.penalty_applied;
    out["solve_backend"] = diagnostic.solve_backend;
    out["derivative_backend"] = diagnostic.derivative_backend;
    return out;
}

py::dict native_regression_contract_to_dict(const NativeRegressionProblemContract& contract) {
    py::dict out;
    py::list parameters;
    for (const auto& spec : contract.parameters) {
        parameters.append(native_regression_parameter_spec_to_dict(spec));
    }
    py::list residual_schema;
    for (const auto& entry : contract.residual_schema) {
        residual_schema.append(native_regression_residual_schema_entry_to_dict(entry));
    }
    out["statuses"] = native_regression_status_names();
    out["parameters"] = parameters;
    out["residual_schema"] = residual_schema;
    out["supported_target_families"] = contract.supported_target_families;
    out["supported_parameter_kinds"] = contract.supported_parameter_kinds;
    out["fixed_shape_residuals"] = contract.fixed_shape_residuals;
    out["production_unsupported_derivative_allowed"] = contract.production_unsupported_derivative_allowed;
    out["row_diagnostic_fields"] = std::vector<std::string>{
        "row_id",
        "success",
        "status",
        "message",
        "residual_start",
        "residual_count",
        "penalty_applied",
        "solve_backend",
        "derivative_backend",
    };
    return out;
}

NativeRegressionResidualRecord native_regression_residual_record_from_dict(const py::dict& input) {
    NativeRegressionResidualRecord record;
    record.row_id = input["row_id"].cast<std::string>();
    record.family = input["family"].cast<std::string>();
    record.target = input["target"].cast<std::string>();
    if (input.contains("row_kind") && !input["row_kind"].is_none()) {
        record.row_kind = input["row_kind"].cast<std::string>();
    }
    if (input.contains("name") && !input["name"].is_none()) {
        record.name = input["name"].cast<std::string>();
    }
    record.predicted = input["predicted"].cast<double>();
    record.observed = input["observed"].cast<double>();
    if (input.contains("scale") && !input["scale"].is_none()) {
        record.scale = input["scale"].cast<double>();
    }
    if (input.contains("sensitivities") && !input["sensitivities"].is_none()) {
        py::dict sensitivities = input["sensitivities"].cast<py::dict>();
        for (const auto& item : sensitivities) {
            record.sensitivities[item.first.cast<std::string>()] = item.second.cast<double>();
        }
    }
    if (input.contains("success") && !input["success"].is_none()) {
        record.success = input["success"].cast<bool>();
    }
    if (input.contains("recoverable_failure") && !input["recoverable_failure"].is_none()) {
        record.recoverable_failure = input["recoverable_failure"].cast<bool>();
    }
    if (input.contains("failure_message") && !input["failure_message"].is_none()) {
        record.failure_message = input["failure_message"].cast<std::string>();
    }
    return record;
}

NativeRegressionParameterSpec native_regression_parameter_spec_from_dict(const py::dict& input) {
    NativeRegressionParameterSpec spec;
    spec.name = input["name"].cast<std::string>();
    if (input.contains("path") && !input["path"].is_none()) {
        spec.path = input["path"].cast<std::string>();
    }
    if (input.contains("kind") && !input["kind"].is_none()) {
        spec.kind = input["kind"].cast<std::string>();
    }
    spec.initial = input["initial"].cast<double>();
    spec.lower = input["lower"].cast<double>();
    spec.upper = input["upper"].cast<double>();
    if (input.contains("scale") && !input["scale"].is_none()) {
        spec.scale = input["scale"].cast<double>();
    }
    if (input.contains("fixed") && !input["fixed"].is_none()) {
        spec.fixed = input["fixed"].cast<bool>();
    }
    if (input.contains("metadata") && !input["metadata"].is_none()) {
        py::dict metadata = input["metadata"].cast<py::dict>();
        for (const auto& item : metadata) {
            spec.metadata[item.first.cast<std::string>()] = item.second.cast<std::string>();
        }
    }
    return spec;
}

std::vector<NativeRegressionResidualRecord> native_regression_residual_records_from_list(const py::list& records) {
    std::vector<NativeRegressionResidualRecord> out;
    out.reserve(py::len(records));
    for (const py::handle item : records) {
        out.push_back(native_regression_residual_record_from_dict(item.cast<py::dict>()));
    }
    return out;
}

std::vector<NativeRegressionParameterSpec> native_regression_parameter_specs_from_list(const py::list& parameters) {
    std::vector<NativeRegressionParameterSpec> out;
    out.reserve(py::len(parameters));
    for (const py::handle item : parameters) {
        out.push_back(native_regression_parameter_spec_from_dict(item.cast<py::dict>()));
    }
    return out;
}

NativeRegressionFitOptions native_regression_fit_options_from_dict(const py::dict& input) {
    NativeRegressionFitOptions options;
    if (input.contains("max_iterations") && !input["max_iterations"].is_none()) {
        options.max_iterations = input["max_iterations"].cast<int>();
    }
    if (input.contains("gradient_tolerance") && !input["gradient_tolerance"].is_none()) {
        options.gradient_tolerance = input["gradient_tolerance"].cast<double>();
    }
    if (input.contains("function_tolerance") && !input["function_tolerance"].is_none()) {
        options.function_tolerance = input["function_tolerance"].cast<double>();
    }
    if (input.contains("parameter_tolerance") && !input["parameter_tolerance"].is_none()) {
        options.parameter_tolerance = input["parameter_tolerance"].cast<double>();
    }
    if (input.contains("penalty_residual") && !input["penalty_residual"].is_none()) {
        options.penalty_residual = input["penalty_residual"].cast<double>();
    }
    if (input.contains("derivative_backend") && !input["derivative_backend"].is_none()) {
        options.derivative_backend = input["derivative_backend"].cast<std::string>();
    }
    if (input.contains("optimizer_backend") && !input["optimizer_backend"].is_none()) {
        options.optimizer_backend = input["optimizer_backend"].cast<std::string>();
    }
    return options;
}

py::dict native_regression_residual_evaluation_to_dict(const NativeRegressionResidualEvaluation& result) {
    py::dict out;
    out["residuals"] = result.residuals;
    py::list residual_schema;
    for (const auto& entry : result.residual_schema) {
        residual_schema.append(native_regression_residual_schema_entry_to_dict(entry));
    }
    py::list row_diagnostics;
    for (const auto& diagnostic : result.row_diagnostics) {
        row_diagnostics.append(native_regression_row_diagnostic_to_dict(diagnostic));
    }
    out["residual_schema"] = residual_schema;
    out["row_diagnostics"] = row_diagnostics;
    out["cost"] = result.cost;
    out["residual_norm"] = result.residual_norm;
    out["success_count"] = result.success_count;
    out["failure_count"] = result.failure_count;
    out["fixed_shape_residuals"] = result.fixed_shape_residuals;
    return out;
}

py::dict native_regression_fit_result_to_dict(const NativeRegressionFitResult& result) {
    py::dict out;
    out["success"] = result.success;
    out["status"] = result.status;
    out["message"] = result.message;
    out["optimizer_backend"] = result.optimizer_backend;
    out["derivative_backend"] = result.derivative_backend;
    out["parameters"] = result.parameters;
    out["parameter_names"] = result.parameter_names;
    out["lower_bounds"] = result.lower_bounds;
    out["upper_bounds"] = result.upper_bounds;
    out["active_bounds"] = result.active_bounds;
    out["initial_cost"] = result.initial_cost;
    out["final_cost"] = result.final_cost;
    out["residual_norm"] = result.residual_norm;
    out["gradient_norm"] = result.gradient_norm;
    out["iterations"] = result.iterations;
    out["function_evaluations"] = result.function_evaluations;
    out["objective_result"] = native_regression_residual_evaluation_to_dict(result.objective_result);
    return out;
}

py::dict evaluate_native_regression_residual_records_binding(const py::list& records, double penalty_residual) {
    NativeRegressionResidualEvaluation result =
        evaluate_native_regression_residual_records(native_regression_residual_records_from_list(records), penalty_residual);
    return native_regression_residual_evaluation_to_dict(result);
}

py::dict solve_native_regression_residual_records_binding(
    const py::list& records,
    const py::list& parameters,
    const py::dict& options
) {
    NativeRegressionFitResult result = solve_native_regression_residual_records(
        native_regression_residual_records_from_list(records),
        native_regression_parameter_specs_from_list(parameters),
        native_regression_fit_options_from_dict(options)
    );
    return native_regression_fit_result_to_dict(result);
}

NativeThermoRegressionTarget native_thermo_regression_target_from_dict(const py::dict& input) {
    NativeThermoRegressionTarget target;
    target.family = input["family"].cast<std::string>();
    if (input.contains("target") && !input["target"].is_none()) {
        target.target = input["target"].cast<std::string>();
    }
    if (input.contains("index") && !input["index"].is_none()) {
        target.index = input["index"].cast<int>();
    }
    target.observed = input["observed"].cast<double>();
    if (input.contains("scale") && !input["scale"].is_none()) {
        target.scale = input["scale"].cast<double>();
    }
    return target;
}

NativeThermoRegressionRow native_thermo_regression_row_from_dict(const py::dict& input) {
    NativeThermoRegressionRow row;
    row.row_id = input["row_id"].cast<std::string>();
    row.row_mode = input["row_mode"].cast<std::string>();
    row.t = input["T"].cast<double>();
    if (input.contains("P") && !input["P"].is_none()) {
        row.p = input["P"].cast<double>();
    }
    if (input.contains("initial_x") && !input["initial_x"].is_none()) {
        row.initial_x = input["initial_x"].cast<std::vector<double>>();
    }
    if (input.contains("x_liq") && !input["x_liq"].is_none()) {
        row.x_liq = input["x_liq"].cast<std::vector<double>>();
    }
    if (input.contains("balance_matrix") && !input["balance_matrix"].is_none()) {
        row.balance_matrix_row_major = input["balance_matrix"].cast<std::vector<double>>();
    }
    if (input.contains("balance_rows") && !input["balance_rows"].is_none()) {
        row.balance_rows = input["balance_rows"].cast<int>();
    }
    if (input.contains("total_vector") && !input["total_vector"].is_none()) {
        row.total_vector = input["total_vector"].cast<std::vector<double>>();
    }
    if (input.contains("reaction_stoichiometry") && !input["reaction_stoichiometry"].is_none()) {
        row.reaction_stoichiometry_row_major = input["reaction_stoichiometry"].cast<std::vector<double>>();
    }
    if (input.contains("reaction_rows") && !input["reaction_rows"].is_none()) {
        row.reaction_rows = input["reaction_rows"].cast<int>();
    }
    if (input.contains("log_equilibrium_constants") && !input["log_equilibrium_constants"].is_none()) {
        row.log_equilibrium_constants = input["log_equilibrium_constants"].cast<std::vector<double>>();
    }
    if (input.contains("reaction_standard_states") && !input["reaction_standard_states"].is_none()) {
        row.reaction_standard_states = input["reaction_standard_states"].cast<std::vector<int>>();
    } else {
        row.reaction_standard_states = std::vector<int>(static_cast<std::size_t>(row.reaction_rows), 0);
    }
    if (input.contains("vapor_species") && !input["vapor_species"].is_none()) {
        row.vapor_species = input["vapor_species"].cast<std::vector<std::string>>();
    }
    if (input.contains("target_species") && !input["target_species"].is_none()) {
        row.targets_species = input["target_species"].cast<std::vector<std::string>>();
    }
    row.speciation_options = chemical_options_from_request(input);
    row.bubble_options = electrolyte_bubble_options_from_request(input);
    py::list targets = input["targets"].cast<py::list>();
    row.targets.reserve(py::len(targets));
    for (const py::handle item : targets) {
        row.targets.push_back(native_thermo_regression_target_from_dict(item.cast<py::dict>()));
    }
    return row;
}

py::dict evaluate_native_thermo_regression_rows_binding(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const py::dict& request
) {
    std::vector<std::string> species = request["species"].cast<std::vector<std::string>>();
    double penalty_residual = 1.0e6;
    if (request.contains("penalty_residual") && !request["penalty_residual"].is_none()) {
        penalty_residual = request["penalty_residual"].cast<double>();
    }
    py::list rows_input = request["rows"].cast<py::list>();
    std::vector<NativeThermoRegressionRow> rows;
    rows.reserve(py::len(rows_input));
    for (const py::handle item : rows_input) {
        rows.push_back(native_thermo_regression_row_from_dict(item.cast<py::dict>()));
    }
    NativeRegressionResidualEvaluation result;
    {
        py::gil_scoped_release release;
        result = evaluate_native_thermo_regression_rows(mixture, species, rows, penalty_residual);
    }
    return native_regression_residual_evaluation_to_dict(result);
}

py::dict fit_native_thermo_regression_binding(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const py::dict& request
) {
    std::vector<std::string> species = request["species"].cast<std::vector<std::string>>();
    py::list rows_input = request["rows"].cast<py::list>();
    std::vector<NativeThermoRegressionRow> rows;
    rows.reserve(py::len(rows_input));
    for (const py::handle item : rows_input) {
        rows.push_back(native_thermo_regression_row_from_dict(item.cast<py::dict>()));
    }
    py::list parameters_input = request["parameters"].cast<py::list>();
    std::vector<NativeRegressionParameterSpec> parameters;
    parameters.reserve(py::len(parameters_input));
    for (const py::handle item : parameters_input) {
        parameters.push_back(native_regression_parameter_spec_from_dict(item.cast<py::dict>()));
    }
    NativeRegressionFitOptions options;
    if (request.contains("options") && !request["options"].is_none()) {
        options = native_regression_fit_options_from_dict(request["options"].cast<py::dict>());
    }
    NativeRegressionFitResult result;
    {
        py::gil_scoped_release release;
        result = fit_native_thermo_regression(mixture, species, rows, parameters, options);
    }
    return native_regression_fit_result_to_dict(result);
}

GenericRegressionRecord generic_record_from_dict(const py::dict& input) {
    GenericRegressionRecord record;
    if (input.contains("term_name") && !input["term_name"].is_none()) {
        record.term_name = input["term_name"].cast<std::string>();
    }
    record.term = input["term"].cast<int>();
    record.t = input["T"].cast<double>();
    record.p = input["P"].cast<double>();
    record.phase = input.contains("phase") ? input["phase"].cast<int>() : 0;
    if (input.contains("x") && !input["x"].is_none()) {
        record.x = input["x"].cast<std::vector<double>>();
    }
    if (input.contains("y") && !input["y"].is_none()) {
        record.y = input["y"].cast<std::vector<double>>();
    }
    if (input.contains("target") && !input["target"].is_none()) {
        record.target = input["target"].cast<double>();
    }
    if (input.contains("target_index") && !input["target_index"].is_none()) {
        record.target_index = input["target_index"].cast<int>();
    }
    if (input.contains("target_index_2") && !input["target_index_2"].is_none()) {
        record.target_index_2 = input["target_index_2"].cast<int>();
    }
    if (input.contains("density_kind") && !input["density_kind"].is_none()) {
        record.density_kind = input["density_kind"].cast<int>();
    }
    if (input.contains("activity_basis") && !input["activity_basis"].is_none()) {
        record.activity_basis = input["activity_basis"].cast<int>();
    }
    if (input.contains("solvent_index") && !input["solvent_index"].is_none()) {
        record.solvent_index = input["solvent_index"].cast<int>();
    }
    if (input.contains("scale") && !input["scale"].is_none()) {
        record.scale = input["scale"].cast<double>();
    }
    return record;
}

std::vector<GenericRegressionRecord> generic_records_from_list(const py::list& records) {
    std::vector<GenericRegressionRecord> out;
    out.reserve(static_cast<std::size_t>(py::len(records)));
    for (py::handle item : records) {
        out.push_back(generic_record_from_dict(item.cast<py::dict>()));
    }
    return out;
}

std::vector<add_args> native_args_from_list(const py::list& args_by_record) {
    std::vector<add_args> out;
    out.reserve(static_cast<std::size_t>(py::len(args_by_record)));
    for (py::handle item : args_by_record) {
        out.push_back(item.cast<add_args>());
    }
    return out;
}

py::dict native_diagnostics_to_dict(
    const std::map<std::string, double>& doubles,
    const std::map<std::string, int>& ints,
    const std::map<std::string, bool>& bools,
    const std::map<std::string, std::string>& strings,
    const std::map<std::string, std::vector<double>>& vectors
) {
    py::dict out;
    for (const auto& item : doubles) {
        out[py::str(item.first)] = item.second;
    }
    for (const auto& item : ints) {
        out[py::str(item.first)] = item.second;
    }
    for (const auto& item : bools) {
        out[py::str(item.first)] = item.second;
    }
    for (const auto& item : strings) {
        out[py::str(item.first)] = item.second;
    }
    for (const auto& item : vectors) {
        out[py::str(item.first)] = item.second;
    }
    if (strings.find("stability_analysis") != strings.end() && strings.at("stability_analysis") == "not_run") {
        out["stability_stable"] = py::none();
    }
    return out;
}

double json_safe_native_double(double value) {
    return std::isfinite(value) ? value : 1.0e300;
}

py::dict native_density_candidate_to_dict(const DensityCandidateDiagnostics& candidate) {
    py::dict out;
    out["rho_sort"] = json_safe_native_double(candidate.rho_sort);
    out["rho"] = json_safe_native_double(candidate.rho);
    out["gres"] = json_safe_native_double(candidate.gres);
    out["rel_resid"] = json_safe_native_double(candidate.rel_resid);
    out["abs_p_error"] = json_safe_native_double(candidate.abs_p_error);
    out["dpdrho"] = json_safe_native_double(candidate.dpdrho);
    out["valid"] = candidate.valid;
    return out;
}

py::dict native_density_diagnostics_to_dict(const DensitySolveDiagnostics& diagnostics) {
    py::dict out;
    out["phase_label"] = diagnostics.phase_label;
    out["phase_kind"] = diagnostics.phase_kind;
    out["T"] = json_safe_native_double(diagnostics.t);
    out["P"] = json_safe_native_double(diagnostics.p);
    out["composition"] = diagnostics.composition;
    out["scan_point_count"] = diagnostics.scan_point_count;
    out["finite_point_count"] = diagnostics.finite_point_count;
    out["coarse_bracket_count"] = diagnostics.coarse_bracket_count;
    out["refined_bracket_count"] = diagnostics.refined_bracket_count;
    out["candidate_root_count"] = diagnostics.candidate_root_count;
    out["best_near_root_pressure_error"] = json_safe_native_double(diagnostics.best_near_root.abs_p_error);
    out["best_near_root"] = native_density_candidate_to_dict(diagnostics.best_near_root);
    out["dpdrho"] = json_safe_native_double(diagnostics.best_near_root.dpdrho);
    out["gres"] = json_safe_native_double(diagnostics.best_near_root.gres);
    out["rejection_reason"] = diagnostics.rejection_reason;
    out["density_fallback_used"] = diagnostics.fallback_used;
    out["density_fallback_rejected_reason"] = diagnostics.fallback_rejected_reason;
    out["density_warm_start_source"] = diagnostics.warm_start_source;
    out["density_validity_gate"] = diagnostics.validity_gate;
    py::list roots;
    for (const auto& candidate : diagnostics.candidate_roots) {
        roots.append(native_density_candidate_to_dict(candidate));
    }
    out["density_candidate_roots"] = roots;
    return out;
}

py::dict native_density_failure_payload(const DensitySolveDiagnostics& diagnostics) {
    py::dict out;
    py::list contexts;
    contexts.append(native_density_diagnostics_to_dict(diagnostics));
    out["density_failure_count"] = diagnostics.validity_gate == "failed" ? 1 : 0;
    out["density_failure_contexts"] = contexts;
    out["density_scan_summary"] = native_density_diagnostics_to_dict(diagnostics);
    out["density_candidate_roots"] = contexts[0].cast<py::dict>()["density_candidate_roots"];
    out["density_best_near_root"] = contexts[0].cast<py::dict>()["best_near_root"];
    out["density_fallback_used"] = diagnostics.fallback_used;
    out["density_fallback_rejected_reason"] = diagnostics.fallback_rejected_reason;
    out["density_warm_start_source"] = diagnostics.warm_start_source;
    out["density_validity_gate"] = diagnostics.validity_gate;
    return out;
}

py::dict native_phase_to_dict(const EquilibriumPhaseNative& phase) {
    py::dict out;
    out["label"] = phase.label;
    out["composition"] = phase.composition;
    out["density"] = phase.density;
    out["temperature"] = phase.temperature;
    out["pressure"] = phase.pressure;
    out["phase_fraction"] = phase.phase_fraction;
    out["ln_fugacity_coefficient"] = phase.ln_fugacity_coefficient;
    out["diagnostics"] = native_diagnostics_to_dict(
        phase.diagnostics_double,
        {},
        {},
        phase.diagnostics_string,
        {}
    );
    return out;
}

py::dict native_trial_to_dict(const StabilityTrialNative& trial) {
    py::dict out;
    out["parent_phase"] = trial.parent_phase;
    out["trial_phase"] = trial.trial_phase;
    out["seed_name"] = trial.seed_name;
    out["composition"] = trial.composition;
    out["tpd"] = trial.tpd;
    out["iterations"] = trial.iterations;
    out["converged"] = trial.converged;
    out["unstable"] = trial.unstable;
    out["diagnostics"] = native_diagnostics_to_dict(
        trial.diagnostics_double,
        {},
        {},
        trial.diagnostics_string,
        {}
    );
    return out;
}

py::list native_attempts_to_list(const std::vector<EquilibriumAttemptDiagnosticsNative>& attempts) {
    auto json_safe_double = [](double value) {
        return std::isfinite(value) ? value : 1.0e300;
    };
    py::list out;
    for (const auto& attempt : attempts) {
        py::dict item;
        item["seed_name"] = attempt.seed_name;
        item["rejection_reason"] = attempt.rejection_reason;
        item["beta_org"] = json_safe_double(attempt.beta_org);
        item["phase_distance"] = json_safe_double(attempt.phase_distance);
        item["solver_residual_norm"] = json_safe_double(attempt.solver_residual_norm);
        item["material_balance_error"] = json_safe_double(attempt.material_balance_error);
        item["charge_balance_error"] = json_safe_double(attempt.charge_balance_error);
        item["gibbs_delta"] = json_safe_double(attempt.gibbs_delta);
        item["iterations"] = attempt.iterations;
        out.append(item);
    }
    return out;
}

py::dict native_stability_to_dict(const StabilityResultNative& result) {
    py::dict out;
    out["result_type"] = "stability";
    out["backend"] = result.backend;
    out["problem_kind"] = result.problem_kind;
    out["stable"] = result.stable;
    out["min_tpd"] = result.min_tpd;
    out["parent_phase"] = result.parent_phase;
    out["trial_phase"] = result.trial_phase;
    out["trial_composition"] = result.trial_composition;
    py::list trials;
    for (const auto& trial : result.trials) {
        trials.append(native_trial_to_dict(trial));
    }
    out["trials"] = trials;
    out["diagnostics"] = native_diagnostics_to_dict(
        result.diagnostics_double,
        result.diagnostics_int,
        result.diagnostics_bool,
        result.diagnostics_string,
        result.diagnostics_vector
    );
    if (result.backend == "electrolyte_tpd") {
        py::dict diagnostics = out["diagnostics"].cast<py::dict>();
        py::dict charge_balance;
        charge_balance["feed"] = diagnostics.contains("phase_charge_balance_feed") ? diagnostics["phase_charge_balance_feed"] : py::float_(0.0);
        charge_balance["trial"] = diagnostics.contains("phase_charge_balance_trial") ? diagnostics["phase_charge_balance_trial"] : py::float_(0.0);
        diagnostics["phase_charge_balance"] = charge_balance;
        out["diagnostics"] = diagnostics;
    }
    return out;
}

py::dict native_equilibrium_to_dict(const EquilibriumResultNative& result) {
    py::dict out;
    out["result_type"] = "equilibrium";
    out["backend"] = result.backend;
    out["problem_kind"] = result.problem_kind;
    py::list phases;
    for (const auto& phase : result.phases) {
        phases.append(native_phase_to_dict(phase));
    }
    out["phases"] = phases;
    out["stable"] = result.stable;
    out["split_detected"] = result.split_detected;
    out["diagnostics"] = native_diagnostics_to_dict(
        result.diagnostics_double,
        result.diagnostics_int,
        result.diagnostics_bool,
        result.diagnostics_string,
        result.diagnostics_vector
    );
    py::dict diagnostics = out["diagnostics"].cast<py::dict>();
    diagnostics["seed_attempts"] = native_attempts_to_list(result.attempt_diagnostics);
    if (!result.density_diagnostics.empty()) {
        py::list contexts;
        for (const auto& density : result.density_diagnostics) {
            contexts.append(native_density_diagnostics_to_dict(density));
        }
        diagnostics["density_failure_contexts"] = contexts;
        diagnostics["density_failure_count"] = static_cast<int>(result.density_diagnostics.size());
        diagnostics["density_scan_summary"] = native_density_diagnostics_to_dict(result.density_diagnostics.back());
        diagnostics["density_candidate_roots"] = native_density_diagnostics_to_dict(result.density_diagnostics.back())["density_candidate_roots"];
        diagnostics["density_best_near_root"] = native_density_diagnostics_to_dict(result.density_diagnostics.back())["best_near_root"];
        diagnostics["density_fallback_used"] = result.density_diagnostics.back().fallback_used;
        diagnostics["density_fallback_rejected_reason"] = result.density_diagnostics.back().fallback_rejected_reason;
        diagnostics["density_warm_start_source"] = result.density_diagnostics.back().warm_start_source;
        diagnostics["density_validity_gate"] = result.density_diagnostics.back().validity_gate;
    }
    out["diagnostics"] = diagnostics;
    return out;
}

py::dict native_chemical_equilibrium_to_dict(const ChemicalEquilibriumResultNative& result) {
    py::dict out;
    out["success"] = result.success;
    out["message"] = result.message;
    out["composition"] = result.composition;
    out["activity_coefficients"] = result.activity_coefficients;
    out["mass_balance_residuals"] = result.mass_balance_residuals;
    out["charge_residual"] = result.charge_residual;
    out["reaction_residuals"] = result.reaction_residuals;
    py::dict diagnostics = native_diagnostics_to_dict(
        result.diagnostics_double,
        result.diagnostics_int,
        result.diagnostics_bool,
        result.diagnostics_string,
        result.diagnostics_vector
    );
    py::dict handoff;
    auto handoff_it = result.diagnostics_vector.find("phase_handoff_composition");
    handoff["composition"] = handoff_it == result.diagnostics_vector.end() ? result.composition : handoff_it->second;
    handoff["activity_coefficients"] = result.activity_coefficients;
    handoff["activity_basis"] = result.diagnostics_string.count("activity_basis")
        ? result.diagnostics_string.at("activity_basis")
        : "mole_fraction";
    diagnostics["phase_equilibrium_handoff"] = handoff;
    out["diagnostics"] = diagnostics;
    return out;
}

py::dict native_chemical_residual_evaluation_to_dict(const ChemicalResidualEvaluationNative& result) {
    py::dict out;
    out["variable_model"] = result.variable_model;
    out["variables"] = result.variables;
    out["lower_bounds"] = result.lower_bounds;
    out["upper_bounds"] = result.upper_bounds;
    out["residual"] = result.residual;
    out["objective"] = result.objective;
    out["gradient"] = result.gradient;
    out["jacobian_row_major"] = result.jacobian_row_major;
    out["jacobian_shape"] = py::make_tuple(result.jacobian_rows, result.jacobian_cols);
    out["jacobian_backend"] = result.diagnostics_string.count("jacobian_backend")
        ? result.diagnostics_string.at("jacobian_backend")
        : "unsupported_derivative";
    out["hessian_backend"] = result.diagnostics_string.count("hessian_backend")
        ? result.diagnostics_string.at("hessian_backend")
        : "gauss_newton";
    out["composition"] = result.composition;
    out["activity_coefficients"] = result.activity_coefficients;
    out["mass_balance_residuals"] = result.mass_balance_residuals;
    out["charge_residual"] = result.charge_residual;
    out["reaction_residuals"] = result.reaction_residuals;
    out["diagnostics"] = native_diagnostics_to_dict(
        result.diagnostics_double,
        result.diagnostics_int,
        result.diagnostics_bool,
        result.diagnostics_string,
        result.diagnostics_vector
    );
    return out;
}

py::dict native_electrolyte_lle_residual_evaluation_to_dict(const ElectrolyteLLEResidualEvaluationNative& result) {
    py::dict out;
    out["variable_model"] = result.variable_model;
    out["variables"] = result.variables;
    out["lower_bounds"] = result.lower_bounds;
    out["upper_bounds"] = result.upper_bounds;
    out["residual"] = result.residual;
    out["objective"] = result.objective;
    out["gradient"] = result.gradient;
    out["jacobian_row_major"] = result.jacobian_row_major;
    out["jacobian_shape"] = py::make_tuple(result.jacobian_rows, result.jacobian_cols);
    out["jacobian_backend"] = result.diagnostics_string.count("jacobian_backend")
        ? result.diagnostics_string.at("jacobian_backend")
        : "unsupported_derivative";
    out["hessian_backend"] = result.diagnostics_string.count("hessian_backend")
        ? result.diagnostics_string.at("hessian_backend")
        : "gauss_newton";
    out["aq_composition"] = result.aq_composition;
    out["org_composition"] = result.org_composition;
    out["aq_ln_fugacity_coefficient"] = result.aq_ln_fugacity_coefficient;
    out["org_ln_fugacity_coefficient"] = result.org_ln_fugacity_coefficient;
    out["aq_density"] = result.aq_density;
    out["org_density"] = result.org_density;
    out["phase_fraction_org"] = result.phase_fraction_org;
    out["material_balance_error"] = result.material_balance_error;
    out["charge_balance_error"] = result.charge_balance_error;
    out["phase_distance"] = result.phase_distance;
    out["gibbs_delta"] = result.gibbs_delta;
    out["diagnostics"] = native_diagnostics_to_dict(
        result.diagnostics_double,
        result.diagnostics_int,
        result.diagnostics_bool,
        result.diagnostics_string,
        result.diagnostics_vector
    );
    return out;
}

[[noreturn]] void raise_native_solution_error_with_diagnostics(
    const std::string& message,
    const EquilibriumResultNative& result
) {
    py::dict diagnostics = native_diagnostics_to_dict(
        result.diagnostics_double,
        result.diagnostics_int,
        result.diagnostics_bool,
        result.diagnostics_string,
        result.diagnostics_vector
    );
    diagnostics["seed_attempts"] = native_attempts_to_list(result.attempt_diagnostics);
    if (!result.density_diagnostics.empty()) {
        py::list contexts;
        for (const auto& density : result.density_diagnostics) {
            contexts.append(native_density_diagnostics_to_dict(density));
        }
        diagnostics["density_failure_contexts"] = contexts;
        diagnostics["density_failure_count"] = static_cast<int>(result.density_diagnostics.size());
        diagnostics["density_scan_summary"] = native_density_diagnostics_to_dict(result.density_diagnostics.back());
        diagnostics["density_candidate_roots"] = native_density_diagnostics_to_dict(result.density_diagnostics.back())["density_candidate_roots"];
        diagnostics["density_best_near_root"] = native_density_diagnostics_to_dict(result.density_diagnostics.back())["best_near_root"];
        diagnostics["density_fallback_used"] = result.density_diagnostics.back().fallback_used;
        diagnostics["density_fallback_rejected_reason"] = result.density_diagnostics.back().fallback_rejected_reason;
        diagnostics["density_warm_start_source"] = result.density_diagnostics.back().warm_start_source;
        diagnostics["density_validity_gate"] = result.density_diagnostics.back().validity_gate;
    }
    py::tuple args(2);
    args[0] = py::str(message);
    args[1] = diagnostics;
    PyErr_SetObject(native_solution_error_type.ptr(), args.ptr());
    throw py::error_already_set();
}

EquilibriumOptionsNative options_from_request(const py::dict& request) {
    EquilibriumOptionsNative options;
    if (!request.contains("options") || request["options"].is_none()) {
        return options;
    }
    py::dict input = request["options"].cast<py::dict>();
    if (input.contains("max_iterations")) {
        options.max_iterations = input["max_iterations"].cast<int>();
    }
    if (input.contains("tolerance")) {
        options.tolerance = input["tolerance"].cast<double>();
    }
    if (input.contains("damping")) {
        options.damping = input["damping"].cast<double>();
    }
    if (input.contains("min_composition")) {
        options.min_composition = input["min_composition"].cast<double>();
    }
    if (input.contains("include_phase_diagnostics")) {
        options.include_phase_diagnostics = input["include_phase_diagnostics"].cast<bool>();
    }
    if (input.contains("stability_precheck")) {
        options.stability_precheck = input["stability_precheck"].cast<bool>();
    }
    if (input.contains("density_diagnostics")) {
        options.density_diagnostics = input["density_diagnostics"].cast<std::string>();
    }
    if (input.contains("experimental_coupled_density_lle")) {
        options.experimental_coupled_density_lle = input["experimental_coupled_density_lle"].cast<bool>();
    }
    if (input.contains("jacobian_backend")) {
        options.jacobian_backend = input["jacobian_backend"].cast<std::string>();
    }
    if (input.contains("timeout_seconds") && !input["timeout_seconds"].is_none()) {
        options.timeout_seconds = input["timeout_seconds"].cast<double>();
    }
    if (input.contains("max_seed_attempts") && !input["max_seed_attempts"].is_none()) {
        options.max_seed_attempts = input["max_seed_attempts"].cast<int>();
    }
    if (input.contains("max_density_failures") && !input["max_density_failures"].is_none()) {
        options.max_density_failures = input["max_density_failures"].cast<int>();
    }
    if (input.contains("max_total_objective_evaluations") && !input["max_total_objective_evaluations"].is_none()) {
        options.max_total_objective_evaluations = input["max_total_objective_evaluations"].cast<int>();
    }
    if (input.contains("return_best_effort")) {
        options.return_best_effort = input["return_best_effort"].cast<bool>();
    }
    return options;
}

ElectrolyteBubbleOptionsNative electrolyte_bubble_options_from_request(const py::dict& request) {
    ElectrolyteBubbleOptionsNative options;
    if (!request.contains("options") || request["options"].is_none()) {
        return options;
    }
    py::dict input = request["options"].cast<py::dict>();
    if (input.contains("initial_pressure")) {
        options.initial_pressure = input["initial_pressure"].cast<double>();
    }
    if (input.contains("min_pressure")) {
        options.min_pressure = input["min_pressure"].cast<double>();
    }
    if (input.contains("max_pressure")) {
        options.max_pressure = input["max_pressure"].cast<double>();
    }
    if (input.contains("max_iterations")) {
        options.max_iterations = input["max_iterations"].cast<int>();
    }
    if (input.contains("max_vapor_iterations")) {
        options.max_vapor_iterations = input["max_vapor_iterations"].cast<int>();
    }
    if (input.contains("max_bracket_expansions")) {
        options.max_bracket_expansions = input["max_bracket_expansions"].cast<int>();
    }
    if (input.contains("tolerance")) {
        options.tolerance = input["tolerance"].cast<double>();
    }
    if (input.contains("vapor_tolerance")) {
        options.vapor_tolerance = input["vapor_tolerance"].cast<double>();
    }
    if (input.contains("pressure_factor")) {
        options.pressure_factor = input["pressure_factor"].cast<double>();
    }
    if (input.contains("min_composition")) {
        options.min_composition = input["min_composition"].cast<double>();
    }
    if (input.contains("charge_tolerance")) {
        options.charge_tolerance = input["charge_tolerance"].cast<double>();
    }
    if (input.contains("return_best_effort")) {
        options.return_best_effort = input["return_best_effort"].cast<bool>();
    }
    if (input.contains("initial_y_vap") && !input["initial_y_vap"].is_none()) {
        options.initial_y_vap = input["initial_y_vap"].cast<std::vector<double>>();
    }
    return options;
}

ChemicalEquilibriumOptionsNative chemical_options_from_request(const py::dict& request) {
    ChemicalEquilibriumOptionsNative options;
    if (!request.contains("options") || request["options"].is_none()) {
        return options;
    }
    py::dict input = request["options"].cast<py::dict>();
    if (input.contains("max_iterations")) {
        options.max_iterations = input["max_iterations"].cast<int>();
    }
    if (input.contains("tolerance")) {
        options.tolerance = input["tolerance"].cast<double>();
    }
    if (input.contains("damping")) {
        options.damping = input["damping"].cast<double>();
    }
    if (input.contains("min_mole_fraction")) {
        options.min_mole_fraction = input["min_mole_fraction"].cast<double>();
    }
    if (input.contains("unsupported_derivative_step")) {
        options.unsupported_derivative_step = input["unsupported_derivative_step"].cast<double>();
    }
    if (input.contains("jacobian_backend")) {
        options.jacobian_backend = input["jacobian_backend"].cast<std::string>();
    }
    if (input.contains("phase")) {
        options.phase = input["phase"].cast<std::string>();
    }
    if (input.contains("activity_output")) {
        options.activity_output = input["activity_output"].cast<std::string>();
    }
    return options;
}

py::dict solve_electrolyte_bubble_native_binding(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const py::dict& request
) {
    double t = request["T"].cast<double>();
    std::vector<double> x_liq = request["x_liq"].cast<std::vector<double>>();
    std::vector<std::string> species = request["species"].cast<std::vector<std::string>>();
    std::vector<std::string> vapor_species = request["vapor_species"].cast<std::vector<std::string>>();
    ElectrolyteBubbleOptionsNative options = electrolyte_bubble_options_from_request(request);
    EquilibriumResultNative result;
    {
        py::gil_scoped_release release;
        result = electrolyte_bubble_pressure_native(mixture, t, x_liq, options, species, vapor_species);
    }
    return native_equilibrium_to_dict(result);
}

py::dict solve_chemical_equilibrium_native_binding(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const py::dict& request
) {
    double t = request["T"].cast<double>();
    double p = request["P"].cast<double>();
    std::vector<double> initial_x = request["initial_x"].cast<std::vector<double>>();
    std::vector<double> balance_matrix = request["balance_matrix"].cast<std::vector<double>>();
    int balance_rows = request["balance_rows"].cast<int>();
    std::vector<double> total_vector = request["total_vector"].cast<std::vector<double>>();
    std::vector<double> reaction_stoichiometry = request["reaction_stoichiometry"].cast<std::vector<double>>();
    int reaction_rows = request["reaction_rows"].cast<int>();
    std::vector<double> log_equilibrium_constants = request["log_equilibrium_constants"].cast<std::vector<double>>();
    std::vector<int> reaction_standard_states;
    if (request.contains("reaction_standard_states") && !request["reaction_standard_states"].is_none()) {
        reaction_standard_states = request["reaction_standard_states"].cast<std::vector<int>>();
    } else {
        reaction_standard_states = std::vector<int>(static_cast<std::size_t>(reaction_rows), 0);
    }
    ChemicalEquilibriumOptionsNative options = chemical_options_from_request(request);
    ChemicalEquilibriumResultNative result;
    {
        py::gil_scoped_release release;
        result = chemical_equilibrium_native(
            mixture,
            t,
            p,
            initial_x,
            balance_matrix,
            balance_rows,
            total_vector,
            reaction_stoichiometry,
            reaction_rows,
            log_equilibrium_constants,
            reaction_standard_states,
            options
        );
    }
    return native_chemical_equilibrium_to_dict(result);
}

py::dict evaluate_chemical_equilibrium_residual_native_binding(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const py::dict& request
) {
    double t = request["T"].cast<double>();
    double p = request["P"].cast<double>();
    std::vector<double> initial_x = request["initial_x"].cast<std::vector<double>>();
    std::vector<double> variables;
    bool has_variables = false;
    if (request.contains("variables") && !request["variables"].is_none()) {
        variables = request["variables"].cast<std::vector<double>>();
        has_variables = true;
    }
    std::vector<double> balance_matrix = request["balance_matrix"].cast<std::vector<double>>();
    int balance_rows = request["balance_rows"].cast<int>();
    std::vector<double> total_vector = request["total_vector"].cast<std::vector<double>>();
    std::vector<double> reaction_stoichiometry = request["reaction_stoichiometry"].cast<std::vector<double>>();
    int reaction_rows = request["reaction_rows"].cast<int>();
    std::vector<double> log_equilibrium_constants = request["log_equilibrium_constants"].cast<std::vector<double>>();
    std::vector<int> reaction_standard_states;
    if (request.contains("reaction_standard_states") && !request["reaction_standard_states"].is_none()) {
        reaction_standard_states = request["reaction_standard_states"].cast<std::vector<int>>();
    } else {
        reaction_standard_states = std::vector<int>(static_cast<std::size_t>(reaction_rows), 0);
    }
    ChemicalEquilibriumOptionsNative options = chemical_options_from_request(request);
    ChemicalResidualEvaluationNative result;
    {
        py::gil_scoped_release release;
        result = evaluate_chemical_equilibrium_residual_native(
            mixture,
            t,
            p,
            initial_x,
            variables,
            has_variables,
            balance_matrix,
            balance_rows,
            total_vector,
            reaction_stoichiometry,
            reaction_rows,
            log_equilibrium_constants,
            reaction_standard_states,
            options
        );
    }
    return native_chemical_residual_evaluation_to_dict(result);
}

std::vector<std::string> string_vector_from_request(const py::dict& request, const char* key, std::vector<std::string> fallback) {
    if (!request.contains(key) || request[key].is_none()) {
        return fallback;
    }
    return request[key].cast<std::vector<std::string>>();
}

py::dict evaluate_electrolyte_lle_residual_native_binding(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const py::dict& request
) {
    double t = request["T"].cast<double>();
    double p = request["P"].cast<double>();
    std::vector<double> feed = request["z"].cast<std::vector<double>>();
    std::vector<std::string> species = string_vector_from_request(request, "species", {});
    EquilibriumOptionsNative options = options_from_request(request);
    std::vector<double> variables;
    bool has_variables = false;
    if (request.contains("variables") && !request["variables"].is_none()) {
        variables = request["variables"].cast<std::vector<double>>();
        has_variables = true;
    }
    std::vector<double> aq;
    std::vector<double> org;
    double beta = 0.5;
    bool has_initial = false;
    if (request.contains("initial_phases") && !request["initial_phases"].is_none()) {
        py::dict initial = request["initial_phases"].cast<py::dict>();
        aq = initial["aq"].cast<std::vector<double>>();
        org = initial["org"].cast<std::vector<double>>();
        beta = initial["phase_fraction"].cast<double>();
        has_initial = true;
    }
    ElectrolyteLLEResidualEvaluationNative result;
    {
        py::gil_scoped_release release;
        result = evaluate_electrolyte_lle_residual_native(
            mixture,
            t,
            p,
            feed,
            options,
            species,
            variables,
            has_variables,
            aq,
            org,
            beta,
            has_initial
        );
    }
    return native_electrolyte_lle_residual_evaluation_to_dict(result);
}

py::dict solve_equilibrium_native_binding(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const py::dict& request
) {
    std::string kind = request["kind"].cast<std::string>();
    double t = request["T"].cast<double>();
    double p = request["P"].cast<double>();
    std::vector<double> feed = request["z"].cast<std::vector<double>>();
    EquilibriumOptionsNative options = options_from_request(request);
    if (kind == "tp_flash" || kind == "neutral_vle") {
        EquilibriumResultNative result;
        {
            py::gil_scoped_release release;
            result = tp_flash_native(mixture, t, p, feed, options);
        }
        return native_equilibrium_to_dict(result);
    }
    if (kind == "lle_flash" || kind == "neutral_lle") {
        std::vector<double> liq1;
        std::vector<double> liq2;
        double beta = 0.5;
        bool has_initial = false;
        if (request.contains("initial_phases") && !request["initial_phases"].is_none()) {
            py::dict initial = request["initial_phases"].cast<py::dict>();
            liq1 = initial["liq1"].cast<std::vector<double>>();
            liq2 = initial["liq2"].cast<std::vector<double>>();
            beta = initial["phase_fraction"].cast<double>();
            has_initial = true;
        }
        EquilibriumResultNative result;
        {
            py::gil_scoped_release release;
            result = lle_flash_native(mixture, t, p, feed, options, liq1, liq2, beta, has_initial);
        }
        return native_equilibrium_to_dict(result);
    }
    if (kind == "stability" || kind == "neutral_tpd") {
        std::vector<std::string> parent_phases = string_vector_from_request(request, "parent_phases", {"liq", "vap"});
        std::vector<std::string> trial_phases = string_vector_from_request(request, "trial_phases", {"liq", "vap"});
        StabilityResultNative result;
        {
            py::gil_scoped_release release;
            result = neutral_stability_native(mixture, t, p, feed, options, parent_phases, trial_phases);
        }
        return native_stability_to_dict(result);
    }
    if (kind == "electrolyte_stability" || kind == "electrolyte_tpd") {
        std::vector<std::string> species;
        if (request.contains("species") && !request["species"].is_none()) {
            species = request["species"].cast<std::vector<std::string>>();
        }
        StabilityResultNative result;
        {
            py::gil_scoped_release release;
            result = electrolyte_stability_native(mixture, t, p, feed, options, species);
        }
        return native_stability_to_dict(result);
    }
    if (kind == "electrolyte_lle" || kind == "electrolyte_lle_flash") {
        std::vector<std::string> species;
        if (request.contains("species") && !request["species"].is_none()) {
            species = request["species"].cast<std::vector<std::string>>();
        }
        std::vector<double> aq;
        std::vector<double> org;
        double beta_org = 0.5;
        bool has_initial = false;
        if (request.contains("initial_phases") && !request["initial_phases"].is_none()) {
            py::dict initial = request["initial_phases"].cast<py::dict>();
            aq = initial["aq"].cast<std::vector<double>>();
            org = initial["org"].cast<std::vector<double>>();
            beta_org = initial["phase_fraction"].cast<double>();
            has_initial = true;
        }
        EquilibriumResultNative result;
        {
            py::gil_scoped_release release;
            result = electrolyte_lle_native(mixture, t, p, feed, options, species, aq, org, beta_org, has_initial);
        }
        const std::string acceptance_gate = result.diagnostics_string["acceptance_gate"];
        if (!result.split_detected
            && !options.return_best_effort
            && (acceptance_gate == "predictive_solve_failed" || acceptance_gate == "predictive_budget_exhausted")) {
            raise_native_solution_error_with_diagnostics("electrolyte LLE flash did not converge", result);
        }
        return native_equilibrium_to_dict(result);
    }
    throw ValueError("Native equilibrium kind is not implemented: " + kind);
}

py::dict fit_pure_neutral_native_least_squares_binding(
    const add_args& args,
    const py::array& density_t,
    const py::array& density_p,
    const py::array& density_rho_exp,
    const py::array& density_phase,
    double density_scale,
    const py::array& vle_t,
    const py::array& vle_p,
    double pure_vle_scale,
    const py::array& x0,
    const py::array& lower,
    const py::array& upper,
    int multistart
) {
    auto density_records = density_records_from_arrays(density_t, density_p, density_rho_exp, density_phase);
    auto pure_vle_records = vle_records_from_arrays(vle_t, vle_p);
    auto cpp_x0 = array_to_double_vector(x0);
    auto cpp_lower = array_to_double_vector(lower);
    auto cpp_upper = array_to_double_vector(upper);
    PureNeutralRegressionResult result;
    {
        py::gil_scoped_release release;
        result = fit_pure_neutral_least_squares_cpp(
            args,
            density_records,
            density_scale,
            pure_vle_records,
            pure_vle_scale,
            cpp_x0,
            cpp_lower,
            cpp_upper,
            multistart
        );
    }
    return regression_result_to_dict(result);
}

py::dict evaluate_pure_neutral_objective_debug_binding(
    const add_args& args,
    const py::array& density_t,
    const py::array& density_p,
    const py::array& density_rho_exp,
    const py::array& density_phase,
    double density_scale,
    const py::array& vle_t,
    const py::array& vle_p,
    double pure_vle_scale,
    const py::array& x
) {
    auto density_records = density_records_from_arrays(density_t, density_p, density_rho_exp, density_phase);
    auto pure_vle_records = vle_records_from_arrays(vle_t, vle_p);
    auto cpp_x = array_to_double_vector(x);
    PureNeutralRegressionDebugResult result;
    {
        py::gil_scoped_release release;
        result = evaluate_pure_neutral_objective_debug_cpp(
            args,
            density_records,
            density_scale,
            pure_vle_records,
            pure_vle_scale,
            cpp_x
        );
    }
    return regression_debug_to_dict(result);
}

py::dict fit_generic_native_least_squares_binding(
    const py::list& args_by_record,
    const py::list& records,
    const py::array& target_kinds,
    const py::array& target_indices,
    const py::array& target_indices_2,
    const py::array& x0,
    const py::array& lower,
    const py::array& upper,
    int multistart,
    int max_nfev
) {
    auto cpp_args = native_args_from_list(args_by_record);
    auto cpp_records = generic_records_from_list(records);
    auto cpp_target_kinds = array_to_int_vector(target_kinds);
    auto cpp_target_indices = array_to_int_vector(target_indices);
    auto cpp_target_indices_2 = array_to_int_vector(target_indices_2);
    auto cpp_x0 = array_to_double_vector(x0);
    auto cpp_lower = array_to_double_vector(lower);
    auto cpp_upper = array_to_double_vector(upper);
    GenericRegressionResult result;
    {
        py::gil_scoped_release release;
        result = fit_generic_least_squares_cpp(
            cpp_args,
            cpp_records,
            cpp_target_kinds,
            cpp_target_indices,
            cpp_target_indices_2,
            cpp_x0,
            cpp_lower,
            cpp_upper,
            multistart,
            max_nfev
        );
    }
    return generic_regression_result_to_dict(result);
}

py::dict evaluate_generic_native_debug_binding(
    const py::list& args_by_record,
    const py::list& records,
    const py::array& target_kinds,
    const py::array& target_indices,
    const py::array& target_indices_2,
    const py::array& x
) {
    auto cpp_args = native_args_from_list(args_by_record);
    auto cpp_records = generic_records_from_list(records);
    auto cpp_target_kinds = array_to_int_vector(target_kinds);
    auto cpp_target_indices = array_to_int_vector(target_indices);
    auto cpp_target_indices_2 = array_to_int_vector(target_indices_2);
    auto cpp_x = array_to_double_vector(x);
    GenericRegressionDebugResult result;
    {
        py::gil_scoped_release release;
        result = evaluate_generic_regression_debug_cpp(
            cpp_args,
            cpp_records,
            cpp_target_kinds,
            cpp_target_indices,
            cpp_target_indices_2,
            cpp_x
        );
    }
    return generic_regression_debug_to_dict(result);
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "pybind11 native backend for epcsaft";

    py::register_exception<ValueError>(m, "NativeValueError");
    py::register_exception<SolutionError>(m, "NativeSolutionError");
    native_solution_error_type = m.attr("NativeSolutionError");

    py::class_<add_args>(m, "NativeArgs")
        .def(py::init<>())
        .def_readwrite("m", &add_args::m)
        .def_readwrite("s", &add_args::s)
        .def_readwrite("e", &add_args::e)
        .def_readwrite("k_ij", &add_args::k_ij)
        .def_readwrite("e_assoc", &add_args::e_assoc)
        .def_readwrite("vol_a", &add_args::vol_a)
        .def_readwrite("z", &add_args::z)
        .def_readwrite("dielc", &add_args::dielc)
        .def_readwrite("mw", &add_args::mw)
        .def_readwrite("mixed_rel_perm_a", &add_args::mixed_rel_perm_a)
        .def_readwrite("mixed_rel_perm_b", &add_args::mixed_rel_perm_b)
        .def_readwrite("mixed_rel_perm_c", &add_args::mixed_rel_perm_c)
        .def_readwrite("mixed_rel_perm_mask", &add_args::mixed_rel_perm_mask)
        .def_readwrite("mixed_rel_perm_water_index", &add_args::mixed_rel_perm_water_index)
        .def_readwrite("dielc_rule", &add_args::dielc_rule)
        .def_readwrite("dielc_diff_mode", &add_args::dielc_diff_mode)
        .def_readwrite("hc_dadx_diff_mode", &add_args::hc_dadx_diff_mode)
        .def_readwrite("disp_dadx_diff_mode", &add_args::disp_dadx_diff_mode)
        .def_readwrite("assoc_dadx_diff_mode", &add_args::assoc_dadx_diff_mode)
        .def_readwrite("d_ion_mode", &add_args::d_ion_mode)
        .def_readwrite("mu_DH_diff_mode", &add_args::mu_DH_diff_mode)
        .def_readwrite("mu_DH_comp_dep_rel_perm", &add_args::mu_DH_comp_dep_rel_perm)
        .def_readwrite("mu_DH_include_sum_term", &add_args::mu_DH_include_sum_term)
        .def_readwrite("include_born_model", &add_args::include_born_model)
        .def_readwrite("d_born_mode", &add_args::d_born_mode)
        .def_readwrite("born_solvation_shell_model", &add_args::born_solvation_shell_model)
        .def_readwrite("born_dielectric_saturation", &add_args::born_dielectric_saturation)
        .def_readwrite("born_bulk_mode", &add_args::born_bulk_mode)
        .def_readwrite("mu_born_diff_mode", &add_args::mu_born_diff_mode)
        .def_readwrite("mu_born_comp_dep_rel_perm", &add_args::mu_born_comp_dep_rel_perm)
        .def_readwrite("mu_born_include_sum_term", &add_args::mu_born_include_sum_term)
        .def_readwrite("mu_born_comp_dep_delta_d", &add_args::mu_born_comp_dep_delta_d)
        .def_readwrite("d_born", &add_args::d_born)
        .def_readwrite("f_solv", &add_args::f_solv)
        .def_readwrite("born_model", &add_args::born_model)
        .def_readwrite("born_radius_model", &add_args::born_radius_model)
        .def_readwrite("born_diff_mode", &add_args::born_diff_mode)
        .def_readwrite("born_eps_mode", &add_args::born_eps_mode)
        .def_readwrite("DH_model", &add_args::DH_model)
        .def_readwrite("assoc_num", &add_args::assoc_num)
        .def_readwrite("assoc_matrix", &add_args::assoc_matrix)
        .def_readwrite("k_hb", &add_args::k_hb)
        .def_readwrite("l_ij", &add_args::l_ij);

    py::class_<ScalarContributionTerms>(m, "ScalarContributionTerms")
        .def_readonly("hc", &ScalarContributionTerms::hc)
        .def_readonly("disp", &ScalarContributionTerms::disp)
        .def_readonly("assoc", &ScalarContributionTerms::assoc)
        .def_readonly("ion", &ScalarContributionTerms::ion)
        .def_readonly("born", &ScalarContributionTerms::born)
        .def_readonly("total", &ScalarContributionTerms::total);

    py::class_<CompressibilityFactorResult>(m, "CompressibilityFactorResult")
        .def_readonly("raw", &CompressibilityFactorResult::raw)
        .def_readonly("terms", &CompressibilityFactorResult::terms);

    py::class_<VectorContributionTerms>(m, "VectorContributionTerms")
        .def_readonly("hc", &VectorContributionTerms::hc)
        .def_readonly("disp", &VectorContributionTerms::disp)
        .def_readonly("assoc", &VectorContributionTerms::assoc)
        .def_readonly("ion", &VectorContributionTerms::ion)
        .def_readonly("born", &VectorContributionTerms::born)
        .def_readonly("total", &VectorContributionTerms::total);

    py::class_<CompositionContributionResult>(m, "CompositionContributionResult")
        .def_readonly("dadx", &CompositionContributionResult::dadx)
        .def_readonly("ares", &CompositionContributionResult::ares)
        .def_readonly("sum_x_dadx", &CompositionContributionResult::sum_x_dadx)
        .def_readonly("z_raw", &CompositionContributionResult::z_raw)
        .def_readonly("z", &CompositionContributionResult::z)
        .def_readonly("derivative_backend", &CompositionContributionResult::derivative_backend)
        .def_readonly("unsupported_derivative_fallback_used", &CompositionContributionResult::unsupported_derivative_fallback_used)
        .def_readonly("unsupported_derivative_fallback_reason", &CompositionContributionResult::unsupported_derivative_fallback_reason);

    py::class_<PressureCompositionDerivativeResult>(m, "PressureCompositionDerivativeResult")
        .def_readonly("dpdx", &PressureCompositionDerivativeResult::dpdx)
        .def_readonly("pressure", &PressureCompositionDerivativeResult::pressure)
        .def_readonly("supported", &PressureCompositionDerivativeResult::supported)
        .def_readonly("derivative_backend", &PressureCompositionDerivativeResult::derivative_backend)
        .def_readonly("unsupported_derivative_fallback_used", &PressureCompositionDerivativeResult::unsupported_derivative_fallback_used)
        .def_readonly("unsupported_derivative_fallback_reason", &PressureCompositionDerivativeResult::unsupported_derivative_fallback_reason);

    py::class_<PressureDensityDerivativeResult>(m, "PressureDensityDerivativeResult")
        .def_readonly("pressure", &PressureDensityDerivativeResult::pressure)
        .def_readonly("dpdrho", &PressureDensityDerivativeResult::dpdrho)
        .def_readonly("supported", &PressureDensityDerivativeResult::supported)
        .def_readonly("derivative_backend", &PressureDensityDerivativeResult::derivative_backend)
        .def_readonly("unsupported_derivative_fallback_used", &PressureDensityDerivativeResult::unsupported_derivative_fallback_used)
        .def_readonly("unsupported_derivative_fallback_reason", &PressureDensityDerivativeResult::unsupported_derivative_fallback_reason);

    py::class_<LnfugCompositionDerivativeResult>(m, "LnfugCompositionDerivativeResult")
        .def_readonly("lnfug", &LnfugCompositionDerivativeResult::lnfug)
        .def_readonly("dlnfugdx_row_major", &LnfugCompositionDerivativeResult::dlnfugdx_row_major)
        .def_readonly("rows", &LnfugCompositionDerivativeResult::rows)
        .def_readonly("cols", &LnfugCompositionDerivativeResult::cols)
        .def_readonly("supported", &LnfugCompositionDerivativeResult::supported)
        .def_readonly("derivative_backend", &LnfugCompositionDerivativeResult::derivative_backend)
        .def_readonly("unsupported_derivative_fallback_used", &LnfugCompositionDerivativeResult::unsupported_derivative_fallback_used)
        .def_readonly("unsupported_derivative_fallback_reason", &LnfugCompositionDerivativeResult::unsupported_derivative_fallback_reason);

    py::class_<ResidualChemicalPotentialResult>(m, "ResidualChemicalPotentialResult")
        .def_readonly("mu", &ResidualChemicalPotentialResult::mu)
        .def_readonly("composition", &ResidualChemicalPotentialResult::composition);

    py::class_<FugacityContributionResult>(m, "FugacityContributionResult")
        .def_readonly("mu", &FugacityContributionResult::mu)
        .def_readonly("lnfugcoef", &FugacityContributionResult::lnfugcoef)
        .def_readonly("composition", &FugacityContributionResult::composition);

    py::class_<ActivityCoefficientNative>(m, "ActivityCoefficientNative")
        .def_readonly("component_activity_coefficients", &ActivityCoefficientNative::component_activity_coefficients)
        .def_readonly("mean_ionic_activity_coefficients_mole_fraction", &ActivityCoefficientNative::mean_ionic_activity_coefficients_mole_fraction)
        .def_readonly("mean_ionic_activity_coefficients_molality", &ActivityCoefficientNative::mean_ionic_activity_coefficients_molality)
        .def_readonly("solvation_free_energy", &ActivityCoefficientNative::solvation_free_energy)
        .def_readonly("pair_molality", &ActivityCoefficientNative::pair_molality)
        .def_readonly("pair_conversion_factor", &ActivityCoefficientNative::pair_conversion_factor)
        .def_readonly("cation_indices", &ActivityCoefficientNative::cation_indices)
        .def_readonly("anion_indices", &ActivityCoefficientNative::anion_indices)
        .def_readonly("solvent_indices", &ActivityCoefficientNative::solvent_indices)
        .def_readonly("pair_cation_indices", &ActivityCoefficientNative::pair_cation_indices)
        .def_readonly("pair_anion_indices", &ActivityCoefficientNative::pair_anion_indices)
        .def_readonly("pair_nu_cation", &ActivityCoefficientNative::pair_nu_cation)
        .def_readonly("pair_nu_anion", &ActivityCoefficientNative::pair_nu_anion)
        .def_readonly("solvent_index", &ActivityCoefficientNative::solvent_index)
        .def_readonly("osmotic_coefficient", &ActivityCoefficientNative::osmotic_coefficient);

    py::class_<ePCSAFTMixtureNative, std::shared_ptr<ePCSAFTMixtureNative>>(m, "NativeMixture")
        .def(py::init<const add_args&>())
        .def("ncomp", &ePCSAFTMixtureNative::ncomp)
        .def("clear_runtime_caches", &ePCSAFTMixtureNative::clear_runtime_caches)
        .def("reset_runtime_cache_stats", &ePCSAFTMixtureNative::reset_runtime_cache_stats)
        .def("reference_state_cache_hits", &ePCSAFTMixtureNative::reference_state_cache_hits)
        .def("reference_state_cache_misses", &ePCSAFTMixtureNative::reference_state_cache_misses)
        .def("density_warm_start_hits", &ePCSAFTMixtureNative::density_warm_start_hits)
        .def("density_warm_start_fallbacks", &ePCSAFTMixtureNative::density_warm_start_fallbacks)
        .def("last_density_diagnostics", [](const ePCSAFTMixtureNative& mixture) {
            return native_density_failure_payload(mixture.last_density_diagnostics());
        });

    py::class_<ePCSAFTStateNative, std::shared_ptr<ePCSAFTStateNative>>(m, "NativeState")
        .def(py::init<
             std::shared_ptr<ePCSAFTMixtureNative>,
             double,
             std::vector<double>,
             int,
             bool,
             double,
             bool,
             double,
             bool,
             double>())
        .def("temperature", &ePCSAFTStateNative::temperature)
        .def("phase", &ePCSAFTStateNative::phase)
        .def("composition", &ePCSAFTStateNative::composition)
        .def("pressure", &ePCSAFTStateNative::pressure)
        .def("density", &ePCSAFTStateNative::density)
        .def("compressibility_factor", &ePCSAFTStateNative::compressibility_factor)
        .def("compressibility_factor_result", &ePCSAFTStateNative::compressibility_factor_result)
        .def("residual_helmholtz", &ePCSAFTStateNative::residual_helmholtz)
        .def("residual_helmholtz_result", &ePCSAFTStateNative::residual_helmholtz_result)
        .def("temperature_derivative_residual_helmholtz", &ePCSAFTStateNative::temperature_derivative_residual_helmholtz)
        .def("temperature_derivative_residual_helmholtz_result", &ePCSAFTStateNative::temperature_derivative_residual_helmholtz_result)
        .def("residual_enthalpy", &ePCSAFTStateNative::residual_enthalpy)
        .def("residual_entropy", &ePCSAFTStateNative::residual_entropy)
        .def("residual_gibbs", &ePCSAFTStateNative::residual_gibbs)
        .def("residual_chemical_potential", &ePCSAFTStateNative::residual_chemical_potential)
        .def("residual_chemical_potential_result", &ePCSAFTStateNative::residual_chemical_potential_result)
        .def("composition_derivative_residual_helmholtz_result", &ePCSAFTStateNative::composition_derivative_residual_helmholtz_result)
        .def("pressure_composition_derivative_result", &ePCSAFTStateNative::pressure_composition_derivative_result)
        .def("pressure_density_derivative_result", &ePCSAFTStateNative::pressure_density_derivative_result)
        .def("lnfug_composition_derivative_result", &ePCSAFTStateNative::lnfug_composition_derivative_result)
        .def("ln_fugacity_coefficient", &ePCSAFTStateNative::ln_fugacity_coefficient)
        .def("fugacity_coefficient", &ePCSAFTStateNative::fugacity_coefficient)
        .def("fugacity_coefficient_result", &ePCSAFTStateNative::fugacity_coefficient_result)
        .def("relative_permittivity", &ePCSAFTStateNative::relative_permittivity)
        .def("osmotic_coefficient", &ePCSAFTStateNative::osmotic_coefficient)
        .def("solvation_free_energy", &ePCSAFTStateNative::solvation_free_energy)
        .def(
            "activity_coefficient_native",
            &ePCSAFTStateNative::activity_coefficient_native,
            py::arg("include_aux") = true,
            py::arg("has_solvent_override") = false,
            py::arg("solvent_override_index") = -1
        );

    m.def("_fit_pure_neutral_native_least_squares", &fit_pure_neutral_native_least_squares_binding);
    m.def("_fit_pure_neutral_native_debug", &evaluate_pure_neutral_objective_debug_binding);
    m.def("_fit_generic_native_least_squares", &fit_generic_native_least_squares_binding);
    m.def("_evaluate_generic_native_debug", &evaluate_generic_native_debug_binding);
    m.def("_native_regression_contract_schema", []() {
        return native_regression_contract_to_dict(native_regression_contract_schema());
    });
    m.def("_native_autodiff_derivative_checks", []() {
        return native_autodiff_derivative_check_to_dict(epcsaft::autodiff::native_autodiff_derivative_checks());
    });
    m.def(
        "_solve_native_implicit_sensitivity",
        [](const std::vector<double>& residual_jacobian_u,
           int residual_rows,
           int state_cols,
           const std::vector<double>& residual_jacobian_theta,
           int parameter_cols) {
            return native_implicit_sensitivity_to_dict(solve_native_implicit_sensitivity(
                residual_jacobian_u,
                residual_rows,
                state_cols,
                residual_jacobian_theta,
                parameter_cols
            ));
        },
        py::arg("residual_jacobian_u"),
        py::arg("residual_rows"),
        py::arg("state_cols"),
        py::arg("residual_jacobian_theta"),
        py::arg("parameter_cols")
    );
    m.def(
        "_evaluate_native_regression_residual_records",
        &evaluate_native_regression_residual_records_binding,
        py::arg("records"),
        py::arg("penalty_residual") = 1.0e6
    );
    m.def(
        "_solve_native_regression_residual_records",
        &solve_native_regression_residual_records_binding,
        py::arg("records"),
        py::arg("parameters"),
        py::arg("options") = py::dict()
    );
    m.def(
        "_evaluate_native_thermo_regression_rows",
        &evaluate_native_thermo_regression_rows_binding,
        py::arg("mixture"),
        py::arg("request")
    );
    m.def(
        "_fit_native_thermo_regression",
        &fit_native_thermo_regression_binding,
        py::arg("mixture"),
        py::arg("request")
    );
    m.def("_solve_equilibrium_native", &solve_equilibrium_native_binding);
    m.def("_evaluate_electrolyte_lle_residual_native", &evaluate_electrolyte_lle_residual_native_binding);
    m.def("_solve_electrolyte_bubble_native", &solve_electrolyte_bubble_native_binding);
    m.def("_solve_chemical_equilibrium_native", &solve_chemical_equilibrium_native_binding);
    m.def("_evaluate_chemical_equilibrium_residual_native", &evaluate_chemical_equilibrium_residual_native_binding);
}



