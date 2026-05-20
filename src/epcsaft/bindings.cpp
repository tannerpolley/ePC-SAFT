#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "association_block.h"
#include "epcsaft_chemical_equilibrium.h"
#include "epcsaft_equilibrium.h"
#include "cppad_smoke_checks.h"
#include "electrolyte_block.h"
#include "eos_phase_block.h"
#include "gibbs_blocks.h"
#include "ipopt_adapter.h"
#include "reaction_block.h"
#include "result_builder.h"
#include "route_builders.h"
#include "route_result_bridge.h"
#include "second_order.h"
#include "stability_route_builders.h"

epcsaft::native::cppad_support::CppADDerivativeResult cppad_eos_contribution_derivatives_cpp(
    double t,
    double rho,
    const std::vector<double>& x,
    const add_args& cppargs
);
epcsaft::native::cppad_support::CppADDerivativeResult cppad_pressure_density_derivative_cpp(
    double t,
    double rho,
    const std::vector<double>& x,
    const add_args& cppargs
);
PhaseStateCompositionSensitivityResult phase_state_ln_fugacity_composition_sensitivity_cpp(
    double t,
    double p,
    std::vector<double> x,
    int phase,
    const add_args& cppargs
);
epcsaft::native::cppad_support::CppADDerivativeResult cppad_pure_neutral_parameter_derivatives_cpp(
    double t,
    double rho,
    const add_args& cppargs
);
NeutralBinaryKijPhaseDerivatives neutral_binary_pair_parameter_phase_derivatives_cpp(
    double t,
    double rho,
    const std::vector<double>& x,
    const add_args& cppargs,
    int parameter_index,
    const std::string& parameter_name
);

namespace py = pybind11;

namespace {

using epcsaft::native::equilibrium_nlp::route_result_bridge::apply_eos_route_metadata_fields;
using epcsaft::native::equilibrium_nlp::route_result_bridge::apply_ipopt_route_solution_fields;
using epcsaft::native::equilibrium_nlp::route_result_bridge::apply_ipopt_route_status_fields;

py::dict cppad_smoke_to_dict(const epcsaft::native::cppad_support::CppADDerivativeResult& result) {
    py::dict out;
    out["cppad_compiled"] = epcsaft::native::cppad_support::cppad_compiled();
    out["supported"] = result.supported;
    out["cppad_used"] = result.supported && result.backend == "cppad";
    out["status"] = epcsaft::native::cppad_support::cppad_build_status();
    out["derivative_backend"] = result.backend;
    out["message"] = result.message;
    out["value"] = result.value;
    out["jacobian_row_major"] = result.jacobian_row_major;
    out["outputs"] = result.outputs;
    out["variables"] = result.variables;
    out["shape"] = py::make_tuple(result.rows, result.cols);
    return out;
}

py::dict phase_state_sensitivity_to_dict(const PhaseStateCompositionSensitivityResult& result) {
    py::dict out;
    out["supported"] = result.supported;
    out["backend"] = result.backend;
    out["derivative_backend"] = result.backend;
    out["density_backend"] = result.density_backend;
    out["message"] = result.message;
    out["temperature"] = result.temperature;
    out["pressure"] = result.pressure;
    out["density"] = result.density;
    out["pressure_density_derivative"] = result.pressure_density_derivative;
    out["pressure_density_second_derivative"] = result.pressure_density_second_derivative;
    out["shape"] = py::make_tuple(result.rows, result.cols);
    out["composition"] = result.composition;
    out["ln_fugacity"] = result.ln_fugacity;
    out["density_composition_derivative"] = result.density_composition_derivative;
    out["density_composition_hessian_row_major"] = result.density_composition_hessian_row_major;
    out["pressure_composition_fixed_density_derivative"] = result.pressure_composition_fixed_density_derivative;
    out["pressure_density_composition_cross_derivative"] = result.pressure_density_composition_cross_derivative;
    out["pressure_composition_fixed_density_hessian_row_major"] =
        result.pressure_composition_fixed_density_hessian_row_major;
    out["ln_fugacity_density_derivative"] = result.ln_fugacity_density_derivative;
    out["fixed_density_jacobian_row_major"] = result.fixed_density_jacobian_row_major;
    out["fixed_density_hessian_tensor_row_major"] = result.fixed_density_hessian_tensor_row_major;
    out["jacobian_row_major"] = result.jacobian_row_major;
    out["hessian_tensor_row_major"] = result.hessian_tensor_row_major;
    return out;
}

py::dict born_ssmds_derivative_to_dict(const BornSSMDSDerivativeResult& result) {
    py::dict out;
    out["supported"] = result.supported;
    out["backend"] = result.backend;
    out["message"] = result.message;
    out["ncomp"] = result.ncomp;
    out["shape"] = py::make_tuple(result.ncomp, result.ncomp);
    out["a_born_d_d_born"] = result.a_born_d_d_born;
    out["a_born_d_f_solv"] = result.a_born_d_f_solv;
    out["mu_res_d_d_born_row_major"] = result.mu_res_d_d_born_row_major;
    out["mu_res_d_f_solv_row_major"] = result.mu_res_d_f_solv_row_major;
    out["lnfug_d_d_born_row_major"] = result.lnfug_d_d_born_row_major;
    out["lnfug_d_f_solv_row_major"] = result.lnfug_d_f_solv_row_major;
    out["lngamma_d_d_born_row_major"] = result.lngamma_d_d_born_row_major;
    out["lngamma_d_f_solv_row_major"] = result.lngamma_d_f_solv_row_major;
    return out;
}

py::dict neutral_binary_kij_property_derivatives_to_dict(
    const NeutralBinaryKijPhaseDerivatives& forward,
    const NeutralBinaryKijPhaseDerivatives& reverse
) {
    if (forward.lnphi.size() != reverse.lnphi.size()
        || forward.dlnphi_dk_fixed_rho.size() != reverse.dlnphi_dk_fixed_rho.size()
        || forward.mu_res.size() != reverse.mu_res.size()
        || forward.dmu_res_dk_fixed_rho.size() != reverse.dmu_res_dk_fixed_rho.size()) {
        throw ValueError("Neutral binary k_ij derivative payloads have inconsistent sizes.");
    }
    std::vector<double> dlnphi_dk;
    std::vector<double> dmu_dk;
    dlnphi_dk.reserve(forward.dlnphi_dk_fixed_rho.size());
    dmu_dk.reserve(forward.dmu_res_dk_fixed_rho.size());
    for (std::size_t i = 0; i < forward.dlnphi_dk_fixed_rho.size(); ++i) {
        dlnphi_dk.push_back(forward.dlnphi_dk_fixed_rho[i] + reverse.dlnphi_dk_fixed_rho[i]);
    }
    for (std::size_t i = 0; i < forward.dmu_res_dk_fixed_rho.size(); ++i) {
        dmu_dk.push_back(forward.dmu_res_dk_fixed_rho[i] + reverse.dmu_res_dk_fixed_rho[i]);
    }
    py::dict out;
    out["supported"] = true;
    out["backend"] = "cppad";
    out["message"] = "CppAD neutral binary k_ij property derivatives available";
    out["pressure"] = forward.pressure;
    out["pressure_d_kij"] = forward.dpdk + reverse.dpdk;
    out["residual_chemical_potential"] = forward.mu_res;
    out["residual_chemical_potential_d_kij_fixed_rho"] = dmu_dk;
    out["ln_fugacity"] = forward.lnphi;
    out["ln_fugacity_d_kij_fixed_rho"] = dlnphi_dk;
    return out;
}

void append_pair_parameter_derivatives(
    py::dict& out,
    const std::string& prefix,
    const NeutralBinaryKijPhaseDerivatives& forward,
    const NeutralBinaryKijPhaseDerivatives& reverse
) {
    if (forward.lnphi.size() != reverse.lnphi.size()
        || forward.dlnphi_dk_fixed_rho.size() != reverse.dlnphi_dk_fixed_rho.size()
        || forward.mu_res.size() != reverse.mu_res.size()
        || forward.dmu_res_dk_fixed_rho.size() != reverse.dmu_res_dk_fixed_rho.size()) {
        throw ValueError("Neutral binary pair-parameter derivative payloads have inconsistent sizes.");
    }
    std::vector<double> dlnphi;
    std::vector<double> dmu;
    dlnphi.reserve(forward.dlnphi_dk_fixed_rho.size());
    dmu.reserve(forward.dmu_res_dk_fixed_rho.size());
    for (std::size_t i = 0; i < forward.dlnphi_dk_fixed_rho.size(); ++i) {
        dlnphi.push_back(forward.dlnphi_dk_fixed_rho[i] + reverse.dlnphi_dk_fixed_rho[i]);
    }
    for (std::size_t i = 0; i < forward.dmu_res_dk_fixed_rho.size(); ++i) {
        dmu.push_back(forward.dmu_res_dk_fixed_rho[i] + reverse.dmu_res_dk_fixed_rho[i]);
    }
    out[(prefix + "_pressure").c_str()] = forward.pressure;
    out[(prefix + "_pressure_derivative").c_str()] = forward.dpdk + reverse.dpdk;
    out[(prefix + "_residual_chemical_potential").c_str()] = forward.mu_res;
    out[(prefix + "_residual_chemical_potential_derivative").c_str()] = dmu;
    out[(prefix + "_ln_fugacity").c_str()] = forward.lnphi;
    out[(prefix + "_ln_fugacity_derivative").c_str()] = dlnphi;
}

py::dict neutral_binary_pair_property_derivatives_to_dict(
    const NeutralBinaryKijPhaseDerivatives& kij_forward,
    const NeutralBinaryKijPhaseDerivatives& kij_reverse,
    const NeutralBinaryKijPhaseDerivatives* lij_forward,
    const NeutralBinaryKijPhaseDerivatives* lij_reverse
) {
    py::dict out;
    out["supported"] = true;
    bool uses_implicit = kij_forward.backend == "cppad_implicit" || kij_reverse.backend == "cppad_implicit";
    if (lij_forward != nullptr && lij_reverse != nullptr) {
        uses_implicit = uses_implicit
            || lij_forward->backend == "cppad_implicit"
            || lij_reverse->backend == "cppad_implicit";
    }
    out["backend"] = uses_implicit ? "cppad_implicit" : "cppad";
    out["message"] = uses_implicit
        ? "CppAD binary pair-parameter derivatives with implicit association value routing available"
        : "CppAD neutral binary pair-parameter property derivatives available";
    append_pair_parameter_derivatives(out, "k_ij", kij_forward, kij_reverse);
    out["parameter_names"] = std::vector<std::string>{"k_ij"};
    if (lij_forward != nullptr && lij_reverse != nullptr) {
        append_pair_parameter_derivatives(out, "l_ij", *lij_forward, *lij_reverse);
        out["parameter_names"] = std::vector<std::string>{"k_ij", "l_ij"};
    }
    return out;
}

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
    out["optimizer_backend"] = result.optimizer_backend.empty() ? result.backend : result.optimizer_backend;
    out["derivative_backend"] = result.derivative_backend.empty() ? result.jacobian_backend : result.derivative_backend;
    out["objective_initial"] = result.initial_cost;
    out["objective_final"] = result.cost;
    out["residual_norm_initial"] = std::sqrt(std::max(0.0, 2.0 * result.initial_cost));
    out["residual_norm_final"] = result.residual_norm;
    out["n_residual_evaluations"] = result.objective_evaluations;
    out["n_jacobian_evaluations"] = result.gradient_evaluations;
    out["gradient_norm"] = result.gradient_norm;
    out["step_norm"] = result.step_norm;
    out["python_objective_used"] = false;
    out["jacobian_available"] = result.jacobian_available;
    out["jacobian_backend"] = result.jacobian_backend;
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
    out["optimizer_backend"] = result.optimizer_backend;
    out["derivative_backend"] = result.derivative_backend;
    out["jacobian_available"] = result.jacobian_available;
    out["jacobian_backend"] = result.jacobian_backend;
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
    return out;
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
    return out;
}

std::string diagnostic_string_or(
    const epcsaft::native::equilibrium_nlp::IpoptSolveResult& result,
    const std::string& key,
    const std::string& default_value
) {
    const auto item = result.diagnostics_string.find(key);
    return item == result.diagnostics_string.end() ? default_value : item->second;
}

int diagnostic_int_or(
    const epcsaft::native::equilibrium_nlp::IpoptSolveResult& result,
    const std::string& key,
    int default_value
) {
    const auto item = result.diagnostics_int.find(key);
    return item == result.diagnostics_int.end() ? default_value : item->second;
}

double diagnostic_double_or(
    const epcsaft::native::equilibrium_nlp::IpoptSolveResult& result,
    const std::string& key,
    double default_value
) {
    const auto item = result.diagnostics_double.find(key);
    return item == result.diagnostics_double.end() ? default_value : item->second;
}

bool diagnostic_bool_or(
    const epcsaft::native::equilibrium_nlp::IpoptSolveResult& result,
    const std::string& key,
    bool default_value
) {
    const auto item = result.diagnostics_bool.find(key);
    return item == result.diagnostics_bool.end() ? default_value : item->second;
}

py::list ipopt_iteration_history_to_list(
    const std::vector<epcsaft::native::equilibrium_nlp::IpoptIterationRecord>& history
) {
    py::list out;
    for (const auto& record : history) {
        py::dict row;
        row["iteration"] = record.iteration;
        row["objective"] = record.objective;
        row["primal_infeasibility"] = record.primal_infeasibility;
        row["dual_infeasibility"] = record.dual_infeasibility;
        row["barrier_parameter"] = record.barrier_parameter;
        row["step_size_primal"] = record.step_size_primal;
        row["step_size_dual"] = record.step_size_dual;
        row["regularization_size"] = record.regularization_size;
        row["step_trial_count"] = record.step_trial_count;
        row["restoration_phase"] = record.restoration_phase;
        out.append(row);
    }
    return out;
}

py::list route_seed_attempts_to_list(
    const std::vector<epcsaft::native::equilibrium_nlp::RouteSeedAttempt>& attempts
) {
    py::list out;
    for (const auto& attempt : attempts) {
        py::dict row;
        row["seed_name"] = attempt.seed_name;
        row["status"] = attempt.status;
        row["solver_status"] = attempt.solver_status;
        row["application_status"] = attempt.application_status;
        row["solver_accepted"] = attempt.solver_accepted;
        row["accepted"] = attempt.accepted;
        row["stable"] = attempt.stable;
        row["iteration_count"] = attempt.iteration_count;
        row["objective"] = attempt.objective;
        row["phase_distance"] = attempt.phase_distance;
        row["material_balance_norm"] = attempt.material_balance_norm;
        row["conserved_balance_norm"] = attempt.conserved_balance_norm;
        row["charge_balance_norm"] = attempt.charge_balance_norm;
        row["pressure_consistency_norm"] = attempt.pressure_consistency_norm;
        row["chemical_potential_consistency_norm"] = attempt.chemical_potential_consistency_norm;
        row["phase_equilibrium_norm"] = attempt.phase_equilibrium_norm;
        row["reaction_stationarity_norm"] = attempt.reaction_stationarity_norm;
        row["min_tpd"] = attempt.min_tpd;
        out.append(row);
    }
    return out;
}

py::dict ipopt_continuation_state_to_dict(
    const epcsaft::native::equilibrium_nlp::IpoptSolveResult& result
) {
    py::dict out;
    out["variables"] = result.variables;
    out["bound_lower_multipliers"] = result.bound_lower_multipliers;
    out["bound_upper_multipliers"] = result.bound_upper_multipliers;
    out["constraint_multipliers"] = result.constraint_multipliers;
    return out;
}

py::dict ipopt_continuation_state_to_dict(
    const std::vector<double>& variables,
    const std::vector<double>& bound_lower_multipliers,
    const std::vector<double>& bound_upper_multipliers,
    const std::vector<double>& constraint_multipliers
) {
    py::dict out;
    out["variables"] = variables;
    out["bound_lower_multipliers"] = bound_lower_multipliers;
    out["bound_upper_multipliers"] = bound_upper_multipliers;
    out["constraint_multipliers"] = constraint_multipliers;
    return out;
}

void apply_ipopt_continuation_state(
    epcsaft::native::equilibrium_nlp::IpoptSolveOptions& options,
    const py::object& continuation_state
) {
    if (continuation_state.is_none()) {
        return;
    }
    const py::dict state = continuation_state.cast<py::dict>();
    if (state.contains("variables")) {
        options.initial_variables = state["variables"].cast<std::vector<double>>();
    }
    if (state.contains("bound_lower_multipliers")) {
        options.initial_bound_lower_multipliers =
            state["bound_lower_multipliers"].cast<std::vector<double>>();
    }
    if (state.contains("bound_upper_multipliers")) {
        options.initial_bound_upper_multipliers =
            state["bound_upper_multipliers"].cast<std::vector<double>>();
    }
    if (state.contains("constraint_multipliers")) {
        options.initial_constraint_multipliers =
            state["constraint_multipliers"].cast<std::vector<double>>();
    }
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
    out["gres"] = json_safe_native_double(diagnostics.best_near_root.gres);
    out["rejection_reason"] = diagnostics.rejection_reason;
    out["density_best_candidate_refinement_used"] = diagnostics.best_candidate_refinement_used;
    out["density_best_candidate_rejection_reason"] = diagnostics.best_candidate_rejection_reason;
    out["density_warm_start_source"] = diagnostics.warm_start_source;
    out["density_validity_gate"] = diagnostics.validity_gate;
    py::list roots;
    for (const auto& candidate : diagnostics.candidate_roots) {
        roots.append(native_density_candidate_to_dict(candidate));
    }
    out["density_candidate_roots"] = roots;
    return out;
}

py::dict eos_phase_block_to_dict(const epcsaft::native::equilibrium_nlp::EosPhaseBlockResult& result) {
    py::dict out;
    out["block"] = result.block;
    out["derivative_backend"] = result.derivative_backend;
    out["variable_names"] = result.variable_names;
    out["constraint_names"] = result.constraint_names;
    out["temperature"] = result.temperature;
    out["target_pressure"] = result.target_pressure;
    out["gas_constant_temperature"] = result.gas_constant_temperature;
    out["total_amount"] = result.total_amount;
    out["volume"] = result.volume;
    out["density"] = result.density;
    out["composition"] = result.composition;
    out["residual_helmholtz"] = result.residual_helmholtz;
    out["eos_pressure"] = result.eos_pressure;
    out["compressibility_factor"] = result.compressibility_factor;
    py::dict objective_terms;
    objective_terms["ideal_helmholtz"] = result.ideal_helmholtz;
    objective_terms["residual_helmholtz"] = result.residual_helmholtz_term;
    objective_terms["pressure_work"] = result.pressure_work;
    out["objective_terms"] = objective_terms;
    py::dict electrolyte_terms;
    electrolyte_terms["block"] = result.electrolyte_contribution.block;
    electrolyte_terms["value_backend"] = result.electrolyte_contribution.value_backend;
    electrolyte_terms["term_basis"] = result.electrolyte_contribution.term_basis;
    electrolyte_terms["active"] = result.electrolyte_contribution.active;
    electrolyte_terms["temperature"] = result.electrolyte_contribution.temperature;
    electrolyte_terms["density"] = result.electrolyte_contribution.density;
    electrolyte_terms["composition"] = result.electrolyte_contribution.composition;
    electrolyte_terms["amounts"] = result.electrolyte_contribution.amounts;
    electrolyte_terms["charges"] = result.electrolyte_contribution.charges;
    electrolyte_terms["phase_charge_residual"] = result.electrolyte_contribution.phase_charge_residual;
    electrolyte_terms["ion_residual_helmholtz"] = result.electrolyte_contribution.ion_residual_helmholtz;
    electrolyte_terms["born_residual_helmholtz"] = result.electrolyte_contribution.born_residual_helmholtz;
    electrolyte_terms["electrolyte_residual_helmholtz"] =
        result.electrolyte_contribution.electrolyte_residual_helmholtz;
    electrolyte_terms["total_residual_helmholtz"] = result.electrolyte_contribution.total_residual_helmholtz;
    out["electrolyte_terms"] = electrolyte_terms;
    out["objective"] = result.objective;
    out["gradient"] = result.gradient;
    out["objective_curvature_backend"] = result.objective_curvature_backend;
    out["objective_curvature_shape"] =
        py::make_tuple(result.objective_curvature_rows, result.objective_curvature_cols);
    out["objective_curvature_row_major"] = result.objective_curvature_row_major;
    out["pressure_hessian_backend"] = result.pressure_hessian_backend;
    out["pressure_hessian_shape"] = py::make_tuple(result.pressure_hessian_rows, result.pressure_hessian_cols);
    out["pressure_hessian_row_major"] = result.pressure_hessian_row_major;
    out["pressure_consistency_residual"] = result.pressure_consistency_residual;
    out["constraint_jacobian_backend"] = result.constraint_jacobian_backend;
    out["constraint_jacobian_shape"] = py::make_tuple(result.constraint_jacobian_rows, result.constraint_jacobian_cols);
    out["constraint_jacobian_row_major"] = result.constraint_jacobian_row_major;
    out["pressure_jacobian_backend"] = result.pressure_jacobian_backend;
    out["pressure_jacobian_shape"] = py::make_tuple(result.pressure_jacobian_rows, result.pressure_jacobian_cols);
    out["pressure_jacobian"] = result.pressure_jacobian_row_major;
    out["pressure_density_derivative"] = result.pressure_density_derivative;
    return out;
}

py::dict electrolyte_contribution_block_to_dict(
    const epcsaft::native::equilibrium_nlp::ElectrolyteContributionBlockResult& result
) {
    py::dict out;
    out["block"] = result.block;
    out["value_backend"] = result.value_backend;
    out["term_basis"] = result.term_basis;
    out["active"] = result.active;
    out["temperature"] = result.temperature;
    out["density"] = result.density;
    out["composition"] = result.composition;
    out["amounts"] = result.amounts;
    out["charges"] = result.charges;
    out["phase_charge_residual"] = result.phase_charge_residual;
    out["ion_residual_helmholtz"] = result.ion_residual_helmholtz;
    out["born_residual_helmholtz"] = result.born_residual_helmholtz;
    out["electrolyte_residual_helmholtz"] = result.electrolyte_residual_helmholtz;
    out["total_residual_helmholtz"] = result.total_residual_helmholtz;
    return out;
}

py::dict eos_phase_system_to_dict(const epcsaft::native::equilibrium_nlp::EosPhaseSystemResult& result) {
    py::dict out;
    out["block"] = result.block;
    out["derivative_backend"] = result.derivative_backend;
    out["phase_count"] = result.phase_count;
    out["species_count"] = result.species_count;
    out["variable_names"] = result.variable_names;
    out["constraint_names"] = result.constraint_names;
    out["temperature"] = result.temperature;
    out["target_pressure"] = result.target_pressure;
    out["feed_amounts"] = result.feed_amounts;
    out["objective"] = result.objective;
    out["association_objective"] = result.association_objective;
    out["phase_association_objectives"] = result.phase_association_objectives;
    out["gradient"] = result.gradient;
    out["objective_hessian_backend"] = result.objective_hessian_backend;
    out["objective_hessian_shape"] = py::make_tuple(result.objective_hessian_rows, result.objective_hessian_cols);
    out["objective_hessian_row_major"] = result.objective_hessian_row_major;
    out["constraints"] = result.constraints;
    out["phase_charge_residuals"] = result.phase_charge_residuals;
    out["phase_association_residuals"] = result.phase_association_residuals;
    out["constraint_jacobian_backend"] = result.constraint_jacobian_backend;
    out["constraint_jacobian_shape"] =
        py::make_tuple(result.constraint_jacobian_rows, result.constraint_jacobian_cols);
    out["constraint_jacobian_row_major"] = result.constraint_jacobian_row_major;
    out["constraint_hessian_backend"] = result.constraint_hessian_backend;
    out["constraint_hessian_shape"] =
        py::make_tuple(result.constraint_hessian_rows, result.constraint_hessian_cols);
    out["constraint_hessian_tensor_row_major"] = result.constraint_hessian_tensor_row_major;
    out["constraint_has_hessian"] = result.constraint_has_hessian;
    py::list phase_blocks;
    for (const auto& block : result.phase_blocks) {
        phase_blocks.append(eos_phase_block_to_dict(block));
    }
    out["phase_blocks"] = phase_blocks;
    return out;
}

py::dict association_mass_action_block_to_dict(
    const epcsaft::native::equilibrium_nlp::AssociationMassActionBlockResult& result
) {
    py::dict out;
    out["block"] = result.block;
    out["derivative_backend"] = result.derivative_backend;
    out["site_count"] = result.site_count;
    out["constraint_names"] = result.constraint_names;
    out["residuals"] = result.residuals;
    out["density_derivative"] = result.density_derivative;
    out["jacobian_shape"] = py::make_tuple(result.jacobian_rows, result.jacobian_cols);
    out["site_fraction_jacobian_row_major"] = result.site_fraction_jacobian_row_major;
    out["site_composition_jacobian_row_major"] = result.site_composition_jacobian_row_major;
    return out;
}

py::dict neutral_two_phase_eos_nlp_contract_to_dict(
    const epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosNlpContract& result
) {
    py::dict out;
    out["problem_name"] = result.problem_name;
    out["derivative_backend"] = result.derivative_backend;
    out["variable_model"] = result.variable_model;
    out["density_backend"] = result.density_backend;
    out["residual_families"] = result.residual_families;
    out["constraint_families"] = result.constraint_families;
    out["phase_count"] = result.phase_count;
    out["species_count"] = result.species_count;
    out["balance_row_count"] = result.balance_row_count;
    out["reaction_count"] = result.reaction_count;
    out["variable_count"] = result.variable_count;
    out["constraint_count"] = result.constraint_count;
    out["jacobian_nonzero_count"] = result.jacobian_nonzero_count;
    out["standard_mu_rt"] = result.standard_mu_rt;
    out["initial_point"] = result.initial_point;
    out["variable_lower_bounds"] = result.variable_lower_bounds;
    out["variable_upper_bounds"] = result.variable_upper_bounds;
    out["constraint_lower_bounds"] = result.constraint_lower_bounds;
    out["constraint_upper_bounds"] = result.constraint_upper_bounds;
    out["objective_at_initial"] = result.objective_at_initial;
    out["gradient_at_initial"] = result.gradient_at_initial;
    out["constraints_at_initial"] = result.constraints_at_initial;
    out["jacobian_rows"] = result.jacobian_rows;
    out["jacobian_cols"] = result.jacobian_cols;
    out["jacobian_values_at_initial"] = result.jacobian_values_at_initial;
    return out;
}

py::dict neutral_two_phase_eos_postsolve_to_dict(
    const epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosPostsolve& result
) {
    py::dict out;
    out["accepted"] = result.accepted;
    out["rejection_reason"] = result.rejection_reason;
    out["derivative_backend"] = result.derivative_backend;
    out["residual_families"] = result.residual_families;
    out["constraint_families"] = result.constraint_families;
    out["phase_count"] = result.phase_count;
    out["species_count"] = result.species_count;
    out["material_balance_norm"] = result.material_balance_norm;
    out["pressure_consistency_norm"] = result.pressure_consistency_norm;
    out["chemical_potential_consistency_norm"] = result.chemical_potential_consistency_norm;
    out["ln_fugacity_consistency_norm"] = result.ln_fugacity_consistency_norm;
    out["charge_balance_norm"] = result.charge_balance_norm;
    out["fixed_composition_norm"] = result.fixed_composition_norm;
    out["phase_amount_total_norm"] = result.phase_amount_total_norm;
    out["phase_distance"] = result.phase_distance;
    out["objective"] = result.objective;
    out["gibbs_feed"] = result.gibbs_feed;
    out["gibbs_split"] = result.gibbs_split;
    out["gibbs_delta"] = result.gibbs_delta;
    out["minimum_phase_fraction"] = result.minimum_phase_fraction;
    out["density_backend"] = result.density_backend;
    out["constraints"] = result.constraints;
    out["phase_amount_totals"] = result.phase_amount_totals;
    out["phase_volumes"] = result.phase_volumes;
    out["phase_densities"] = result.phase_densities;
    out["phase_compositions"] = result.phase_compositions;
    out["phase_ln_fugacity_coefficients"] = result.phase_ln_fugacity_coefficients;
    return out;
}

py::dict neutral_two_phase_eos_route_result_to_dict(
    const epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosRouteResult& result
) {
    py::dict out;
    apply_ipopt_route_status_fields(out, result);
    apply_eos_route_metadata_fields(out, result);
    out["solver_feasible_point"] = result.solver_feasible_point;
    apply_ipopt_route_solution_fields(out, result);
    out["continuation_state"] = ipopt_continuation_state_to_dict(
        result.variables,
        result.bound_lower_multipliers,
        result.bound_upper_multipliers,
        result.constraint_multipliers
    );
    out["iteration_history"] = ipopt_iteration_history_to_list(result.iteration_history);
    out["seed_attempts"] = route_seed_attempts_to_list(result.seed_attempts);
    out["phase_amounts"] = result.phase_amounts;
    out["phase_volumes"] = result.phase_volumes;
    out["postsolve"] = neutral_two_phase_eos_postsolve_to_dict(result.postsolve);
    return out;
}

py::dict reactive_two_phase_eos_postsolve_to_dict(
    const epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosPostsolve& result
) {
    py::dict out;
    out["accepted"] = result.accepted;
    out["rejection_reason"] = result.rejection_reason;
    out["derivative_backend"] = result.derivative_backend;
    out["density_backend"] = result.density_backend;
    out["residual_families"] = result.residual_families;
    out["constraint_families"] = result.constraint_families;
    out["phase_count"] = result.phase_count;
    out["species_count"] = result.species_count;
    out["balance_row_count"] = result.balance_row_count;
    out["reaction_count"] = result.reaction_count;
    out["conserved_balance_norm"] = result.conserved_balance_norm;
    out["charge_balance_norm"] = result.charge_balance_norm;
    out["pressure_consistency_norm"] = result.pressure_consistency_norm;
    out["phase_equilibrium_norm"] = result.phase_equilibrium_norm;
    out["reaction_stationarity_norm"] = result.reaction_stationarity_norm;
    out["phase_distance"] = result.phase_distance;
    out["objective"] = result.objective;
    out["standard_mu_rt"] = result.standard_mu_rt;
    out["constraints"] = result.constraints;
    out["reaction_stationarity_residuals"] = result.reaction_stationarity_residuals;
    out["phase_equilibrium_residuals"] = result.phase_equilibrium_residuals;
    out["phase_charge_residuals"] = result.phase_charge_residuals;
    out["phase_amount_totals"] = result.phase_amount_totals;
    out["phase_volumes"] = result.phase_volumes;
    out["phase_densities"] = result.phase_densities;
    out["phase_compositions"] = result.phase_compositions;
    out["phase_ln_fugacity_coefficients"] = result.phase_ln_fugacity_coefficients;
    return out;
}

py::dict reactive_two_phase_eos_route_result_to_dict(
    const epcsaft::native::equilibrium_nlp::ReactiveTwoPhaseEosRouteResult& result
) {
    py::dict out;
    apply_ipopt_route_status_fields(out, result);
    apply_eos_route_metadata_fields(out, result);
    out["phase_count"] = result.phase_count;
    out["species_count"] = result.species_count;
    out["balance_row_count"] = result.balance_row_count;
    out["reaction_count"] = result.reaction_count;
    out["standard_mu_rt"] = result.standard_mu_rt;
    apply_ipopt_route_solution_fields(out, result);
    out["continuation_state"] = ipopt_continuation_state_to_dict(
        result.variables,
        result.bound_lower_multipliers,
        result.bound_upper_multipliers,
        result.constraint_multipliers
    );
    out["iteration_history"] = ipopt_iteration_history_to_list(result.iteration_history);
    out["seed_attempts"] = route_seed_attempts_to_list(result.seed_attempts);
    out["phase_amounts"] = result.phase_amounts;
    out["phase_volumes"] = result.phase_volumes;
    out["postsolve"] = reactive_two_phase_eos_postsolve_to_dict(result.postsolve);
    return out;
}

int stability_phase_token_to_int(const std::string& phase) {
    if (phase == "liq" || phase == "liquid") {
        return 0;
    }
    if (phase == "vap" || phase == "vapor") {
        return 1;
    }
    throw ValueError("stability phase must be 'liq' or 'vap'.");
}

py::dict stability_nlp_contract_to_dict(
    const epcsaft::native::equilibrium_nlp::StabilityNlpContract& result
) {
    py::dict out;
    out["problem_name"] = result.problem_name;
    out["derivative_backend"] = result.derivative_backend;
    out["species_count"] = result.species_count;
    out["variable_count"] = result.variable_count;
    out["constraint_count"] = result.constraint_count;
    out["jacobian_nonzero_count"] = result.jacobian_nonzero_count;
    out["parent_phase"] = result.parent_phase;
    out["trial_phase"] = result.trial_phase;
    out["feed_composition"] = result.feed_composition;
    out["parent_reduced_potential"] = result.parent_reduced_potential;
    out["initial_point"] = result.initial_point;
    out["variable_lower_bounds"] = result.variable_lower_bounds;
    out["variable_upper_bounds"] = result.variable_upper_bounds;
    out["constraint_lower_bounds"] = result.constraint_lower_bounds;
    out["constraint_upper_bounds"] = result.constraint_upper_bounds;
    out["objective_at_initial"] = result.objective_at_initial;
    out["gradient_at_initial"] = result.gradient_at_initial;
    out["constraints_at_initial"] = result.constraints_at_initial;
    out["jacobian_rows"] = result.jacobian_rows;
    out["jacobian_cols"] = result.jacobian_cols;
    out["jacobian_values_at_initial"] = result.jacobian_values_at_initial;
    return out;
}

py::dict stability_route_result_to_dict(
    const epcsaft::native::equilibrium_nlp::StabilityRouteResult& result
) {
    py::dict out;
    apply_ipopt_route_status_fields(out, result);
    out["stable"] = result.stable;
    out["parent_phase"] = result.parent_phase;
    out["trial_phase"] = result.trial_phase;
    out["min_tpd"] = result.min_tpd;
    apply_ipopt_route_solution_fields(out, result);
    out["continuation_state"] = ipopt_continuation_state_to_dict(
        result.variables,
        result.bound_lower_multipliers,
        result.bound_upper_multipliers,
        result.constraint_multipliers
    );
    out["iteration_history"] = ipopt_iteration_history_to_list(result.iteration_history);
    out["seed_attempts"] = route_seed_attempts_to_list(result.seed_attempts);
    out["trial_composition"] = result.trial_composition;
    out["initial_composition"] = result.initial_composition;
    out["parent_reduced_potential"] = result.parent_reduced_potential;
    return out;
}

py::dict neutral_two_phase_eos_phase_payload_to_dict(
    const epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosPhasePayload& phase
) {
    py::dict diagnostics;
    diagnostics["amount_total"] = phase.amount_total;
    diagnostics["volume"] = phase.volume;
    diagnostics["eos_pressure"] = phase.eos_pressure;
    diagnostics["pressure_consistency_residual"] = phase.pressure_consistency_residual;
    diagnostics["compressibility_factor"] = phase.compressibility_factor;

    py::dict out;
    out["label"] = phase.label;
    out["composition"] = phase.composition;
    out["density"] = phase.density;
    out["temperature"] = phase.temperature;
    out["pressure"] = phase.pressure;
    out["phase_fraction"] = phase.phase_fraction;
    out["ln_fugacity_coefficient"] = phase.ln_fugacity_coefficient;
    out["fugacity_coefficient"] = phase.fugacity_coefficient;
    out["amount_total"] = phase.amount_total;
    out["volume"] = phase.volume;
    out["eos_pressure"] = phase.eos_pressure;
    out["diagnostics"] = diagnostics;
    return out;
}

py::dict neutral_two_phase_eos_result_payload_to_dict(
    const epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosResultPayload& result
) {
    py::list phases;
    py::list labels;
    for (const auto& phase : result.phases) {
        phases.append(neutral_two_phase_eos_phase_payload_to_dict(phase));
        labels.append(phase.label);
    }

    py::dict diagnostics;
    diagnostics["derivative_backend"] = result.derivative_backend;
    diagnostics["rejection_reason"] = result.rejection_reason;
    diagnostics["objective"] = result.objective;
    diagnostics["material_balance_norm"] = result.material_balance_norm;
    diagnostics["pressure_consistency_norm"] = result.pressure_consistency_norm;
    diagnostics["chemical_potential_consistency_norm"] = result.chemical_potential_consistency_norm;
    diagnostics["ln_fugacity_consistency_norm"] = result.ln_fugacity_consistency_norm;
    diagnostics["phase_distance"] = result.phase_distance;
    diagnostics["constraints"] = result.constraints;
    diagnostics["feed_amounts"] = result.feed_amounts;

    py::dict out;
    out["accepted"] = result.accepted;
    out["backend"] = result.backend;
    out["problem_kind"] = result.problem_kind;
    out["phase_labels"] = labels;
    out["phases"] = phases;
    out["stable"] = result.stable;
    out["split_detected"] = result.split_detected;
    out["derivative_backend"] = result.derivative_backend;
    out["rejection_reason"] = result.rejection_reason;
    out["objective"] = result.objective;
    out["material_balance_norm"] = result.material_balance_norm;
    out["pressure_consistency_norm"] = result.pressure_consistency_norm;
    out["chemical_potential_consistency_norm"] = result.chemical_potential_consistency_norm;
    out["ln_fugacity_consistency_norm"] = result.ln_fugacity_consistency_norm;
    out["phase_distance"] = result.phase_distance;
    out["constraints"] = result.constraints;
    out["diagnostics"] = diagnostics;
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
    out["density_best_candidate_refinement_used"] = diagnostics.best_candidate_refinement_used;
    out["density_best_candidate_rejection_reason"] = diagnostics.best_candidate_rejection_reason;
    out["density_warm_start_source"] = diagnostics.warm_start_source;
    out["density_validity_gate"] = diagnostics.validity_gate;
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
    diagnostics["continuation_state"] = ipopt_continuation_state_to_dict(
        result.continuation_variables,
        result.continuation_bound_lower_multipliers,
        result.continuation_bound_upper_multipliers,
        result.continuation_constraint_multipliers
    );
    diagnostics["iteration_history"] = ipopt_iteration_history_to_list(result.iteration_history);
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
        : "unspecified";
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
    out["density_backend"] = result.density_backend;
    out["residual_families"] = result.residual_families;
    out["constraint_families"] = result.constraint_families;
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
        : "unspecified";
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

py::dict native_reactive_phase_residual_evaluation_to_dict(const ReactivePhaseResidualEvaluationNative& result) {
    py::dict out;
    out["variable_model"] = result.variable_model;
    out["density_backend"] = result.density_backend;
    out["residual_families"] = result.residual_families;
    out["constraint_families"] = result.constraint_families;
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
        : "unspecified";
    out["phase1_composition"] = result.phase1_composition;
    out["phase2_composition"] = result.phase2_composition;
    out["phase1_amounts"] = result.phase1_amounts;
    out["phase2_amounts"] = result.phase2_amounts;
    out["phase1_ln_fugacity_coefficient"] = result.phase1_ln_fugacity_coefficient;
    out["phase2_ln_fugacity_coefficient"] = result.phase2_ln_fugacity_coefficient;
    out["phase1_density"] = result.phase1_density;
    out["phase2_density"] = result.phase2_density;
    out["phase_fraction_phase2"] = result.phase_fraction_phase2;
    out["element_balance_residuals"] = result.element_balance_residuals;
    out["reaction_residuals_phase1"] = result.reaction_residuals_phase1;
    out["reaction_residuals_phase2"] = result.reaction_residuals_phase2;
    out["reaction_residuals_cross_phase"] = result.reaction_residuals_cross_phase;
    out["neutral_phase_equilibrium_residuals"] = result.neutral_phase_equilibrium_residuals;
    out["ionic_equilibrium_residuals"] = result.ionic_equilibrium_residuals;
    out["phase_charge_residuals"] = result.phase_charge_residuals;
    out["phase_distance"] = result.phase_distance;
    out["diagnostics"] = native_diagnostics_to_dict(
        result.diagnostics_double,
        result.diagnostics_int,
        result.diagnostics_bool,
        result.diagnostics_string,
        result.diagnostics_vector
    );
    return out;
}

EquilibriumOptionsNative options_from_request(const py::dict& request) {
    EquilibriumOptionsNative options;
    if (!request.contains("options") || request["options"].is_none()) {
        return options;
    }
    py::dict input = request["options"].cast<py::dict>();
    if (input.contains("min_composition")) {
        options.min_composition = input["min_composition"].cast<double>();
    }
    if (input.contains("jacobian_backend")) {
        options.jacobian_backend = input["jacobian_backend"].cast<std::string>();
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
    if (input.contains("min_mole_fraction")) {
        options.min_mole_fraction = input["min_mole_fraction"].cast<double>();
    }
    if (input.contains("jacobian_backend")) {
        options.jacobian_backend = input["jacobian_backend"].cast<std::string>();
    }
    if (input.contains("solver_backend")) {
        options.solver_backend = input["solver_backend"].cast<std::string>();
    }
    if (input.contains("hessian_mode")) {
        options.hessian_mode = input["hessian_mode"].cast<std::string>();
    }
    if (input.contains("iteration_history_limit")) {
        options.iteration_history_limit = input["iteration_history_limit"].cast<int>();
    }
    if (input.contains("linear_solver")) {
        options.linear_solver = input["linear_solver"].cast<std::string>();
    }
    if (input.contains("acceptable_tolerance")) {
        options.acceptable_tolerance = input["acceptable_tolerance"].cast<double>();
    }
    if (input.contains("constraint_violation_tolerance")) {
        options.constraint_violation_tolerance = input["constraint_violation_tolerance"].cast<double>();
    }
    if (input.contains("dual_infeasibility_tolerance")) {
        options.dual_infeasibility_tolerance = input["dual_infeasibility_tolerance"].cast<double>();
    }
    if (input.contains("complementarity_tolerance")) {
        options.complementarity_tolerance = input["complementarity_tolerance"].cast<double>();
    }
    if (input.contains("continuation_state") && !input["continuation_state"].is_none()) {
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions solve_options;
        apply_ipopt_continuation_state(solve_options, input["continuation_state"]);
        options.initial_variables = solve_options.initial_variables;
        options.initial_bound_lower_multipliers = solve_options.initial_bound_lower_multipliers;
        options.initial_bound_upper_multipliers = solve_options.initial_bound_upper_multipliers;
        options.initial_constraint_multipliers = solve_options.initial_constraint_multipliers;
    }
    if (input.contains("phase")) {
        options.phase = input["phase"].cast<std::string>();
    }
    if (input.contains("activity_output")) {
        options.activity_output = input["activity_output"].cast<std::string>();
    }
    return options;
}

epcsaft::native::equilibrium_nlp::IpoptSolveOptions ipopt_solve_options_from_scalars(
    int max_iterations,
    double tolerance,
    double timeout_seconds,
    const std::string& hessian_mode = "auto",
    int iteration_history_limit = 20,
    const std::string& linear_solver = "auto",
    double acceptable_tolerance = 0.0,
    double constraint_violation_tolerance = 0.0,
    double dual_infeasibility_tolerance = 0.0,
    double complementarity_tolerance = 0.0
) {
    epcsaft::native::equilibrium_nlp::IpoptSolveOptions options;
    options.max_iterations = max_iterations;
    options.tolerance = tolerance;
    options.acceptable_tolerance = acceptable_tolerance;
    options.constraint_violation_tolerance = constraint_violation_tolerance;
    options.dual_infeasibility_tolerance = dual_infeasibility_tolerance;
    options.complementarity_tolerance = complementarity_tolerance;
    options.max_wall_time_seconds = timeout_seconds;
    options.hessian_mode = hessian_mode;
    options.iteration_history_limit = iteration_history_limit;
    options.linear_solver = linear_solver;
    return options;
}

void apply_ipopt_control_kwargs(
    epcsaft::native::equilibrium_nlp::IpoptSolveOptions& options,
    const py::kwargs& kwargs
) {
    for (const auto& item : kwargs) {
        const std::string key = py::cast<std::string>(item.first);
        if (key == "linear_solver") {
            options.linear_solver = py::cast<std::string>(item.second);
            continue;
        }
        if (key == "acceptable_tolerance") {
            options.acceptable_tolerance = py::cast<double>(item.second);
            continue;
        }
        if (key == "constraint_violation_tolerance") {
            options.constraint_violation_tolerance = py::cast<double>(item.second);
            continue;
        }
        if (key == "dual_infeasibility_tolerance") {
            options.dual_infeasibility_tolerance = py::cast<double>(item.second);
            continue;
        }
        if (key == "complementarity_tolerance") {
            options.complementarity_tolerance = py::cast<double>(item.second);
            continue;
        }
        if (key == "print_level") {
            options.print_level = py::cast<int>(item.second);
            continue;
        }
        throw ValueError("Unknown Ipopt solve option keyword: " + key);
    }
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

std::vector<std::string> string_vector_from_request(const py::dict& request, const char* key, std::vector<std::string> default_value) {
    if (!request.contains(key) || request[key].is_none()) {
        return default_value;
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

py::dict evaluate_reactive_phase_equilibrium_residual_native_binding(
    const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
    const py::dict& request
) {
    double t = request["T"].cast<double>();
    double p = request["P"].cast<double>();
    std::vector<double> feed = request["z"].cast<std::vector<double>>();
    EquilibriumOptionsNative options = options_from_request(request);
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
    std::vector<double> reaction_phase_stoichiometry;
    if (request.contains("reaction_phase_stoichiometry") && !request["reaction_phase_stoichiometry"].is_none()) {
        reaction_phase_stoichiometry = request["reaction_phase_stoichiometry"].cast<std::vector<double>>();
    }
    std::vector<double> variables;
    bool has_variables = false;
    if (request.contains("variables") && !request["variables"].is_none()) {
        variables = request["variables"].cast<std::vector<double>>();
        has_variables = true;
    }
    std::vector<double> phase1;
    std::vector<double> phase2;
    double beta2 = 0.5;
    bool has_initial = false;
    if (request.contains("initial_phases") && !request["initial_phases"].is_none()) {
        py::dict initial = request["initial_phases"].cast<py::dict>();
        if (initial.contains("liq1")) {
            phase1 = initial["liq1"].cast<std::vector<double>>();
        } else if (initial.contains("aq")) {
            phase1 = initial["aq"].cast<std::vector<double>>();
        } else {
            throw ValueError("initial reactive phases require liq1/liq2 or aq/org keys.");
        }
        if (initial.contains("liq2")) {
            phase2 = initial["liq2"].cast<std::vector<double>>();
        } else if (initial.contains("org")) {
            phase2 = initial["org"].cast<std::vector<double>>();
        } else {
            throw ValueError("initial reactive phases require liq1/liq2 or aq/org keys.");
        }
        beta2 = initial["phase_fraction"].cast<double>();
        has_initial = true;
    }
    ReactivePhaseResidualEvaluationNative result;
    {
        py::gil_scoped_release release;
        result = evaluate_reactive_phase_equilibrium_residual_native(
            mixture,
            t,
            p,
            feed,
            options,
            balance_matrix,
            balance_rows,
            total_vector,
            reaction_stoichiometry,
            reaction_rows,
            log_equilibrium_constants,
            reaction_standard_states,
            reaction_phase_stoichiometry,
            variables,
            has_variables,
            phase1,
            phase2,
            beta2,
            has_initial
        );
    }
    return native_reactive_phase_residual_evaluation_to_dict(result);
}

py::dict fit_pure_neutral_native_ceres_binding(
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
    const py::array& upper
) {
    auto density_records = density_records_from_arrays(density_t, density_p, density_rho_exp, density_phase);
    auto pure_vle_records = vle_records_from_arrays(vle_t, vle_p);
    auto cpp_x0 = array_to_double_vector(x0);
    auto cpp_lower = array_to_double_vector(lower);
    auto cpp_upper = array_to_double_vector(upper);
    PureNeutralRegressionResult result;
    {
        py::gil_scoped_release release;
        result = fit_pure_neutral_ceres_cpp(
            args,
            density_records,
            density_scale,
            pure_vle_records,
            pure_vle_scale,
            cpp_x0,
            cpp_lower,
            cpp_upper
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

py::dict fit_generic_native_ceres_binding(
    const py::list& args_by_record,
    const py::list& records,
    const py::array& target_kinds,
    const py::array& target_indices,
    const py::array& target_indices_2,
    const py::array& x0,
    const py::array& lower,
    const py::array& upper,
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
        result = fit_generic_ceres_cpp(
            cpp_args,
            cpp_records,
            cpp_target_kinds,
            cpp_target_indices,
            cpp_target_indices_2,
            cpp_x0,
            cpp_lower,
            cpp_upper,
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

    m.def("_native_cppad_smoke", []() {
        return cppad_smoke_to_dict(epcsaft::native::cppad_support::cppad_square_smoke_derivative(3.0));
    });
    m.def("_native_ceres_smoke", []() {
        py::dict out;
#ifdef EPCSAFT_HAS_CERES
        const bool compiled = true;
#else
        const bool compiled = false;
#endif
        out["backend"] = "ceres";
        out["compiled"] = compiled;
        out["available"] = compiled;
        out["status"] = compiled ? "enabled_available" : "disabled";
        return out;
    });
    m.def("_native_second_order_assembly_smoke", []() {
        namespace nlp = epcsaft::native::equilibrium_nlp;
        py::dict out;

        nlp::ObjectiveSecondOrderData objective;
        objective.variable_count = 2;
        objective.hessian_row_major = {2.0, 1.0, 1.0, 4.0};
        objective.backend = "analytic";

        nlp::ConstraintSecondOrderData constraints;
        constraints.constraint_count = 2;
        constraints.variable_count = 2;
        constraints.hessian_tensor_row_major = {
            0.0, 2.0,
            2.0, 0.0,
            0.0, 0.0,
            0.0, 0.0
        };
        constraints.has_hessian = {true, false};
        constraints.backend = "analytic";

        const nlp::LagrangianHessianAssembler assembler(2);
        out["nonzero_count"] = assembler.nonzero_count();
        out["lagrangian_lower"] = assembler.values(0.5, objective, constraints, {3.0, 100.0});

        nlp::ConstraintSecondOrderData nearly_symmetric_constraints;
        nearly_symmetric_constraints.constraint_count = 1;
        nearly_symmetric_constraints.variable_count = 2;
        nearly_symmetric_constraints.hessian_tensor_row_major = {
            0.0, 1.0 + 5.0e-9,
            1.0 - 5.0e-9, 0.0,
        };
        nearly_symmetric_constraints.has_hessian = {true};
        nearly_symmetric_constraints.backend = "analytic";
        out["weighted_symmetry_lower"] =
            assembler.values(0.0, objective, nearly_symmetric_constraints, {1.0e12});

        nlp::ResidualSecondOrderData residuals;
        residuals.residual_count = 2;
        residuals.variable_count = 2;
        residuals.residuals = {2.0, -1.0};
        residuals.jacobian_row_major = {
            1.0, 2.0,
            3.0, 4.0
        };
        residuals.residual_hessian_tensor_row_major = {
            1.0, 0.0,
            0.0, 2.0,
            0.0, 1.0,
            1.0, 0.0
        };
        const nlp::ObjectiveSecondOrderData residual_quadratic =
            nlp::residual_quadratic_objective_second_order(residuals);
        out["residual_quadratic_lower"] =
            nlp::lower_triangle_values(residual_quadratic.hessian_row_major, residual_quadratic.variable_count);

        nlp::VariableTransformSecondOrderData transform;
        transform.input_variable_count = 2;
        transform.output_variable_count = 2;
        transform.jacobian_row_major = {
            4.0, 0.0,
            0.0, 1.0
        };
        transform.output_hessian_tensor_row_major = {
            2.0, 0.0,
            0.0, 0.0,
            0.0, 0.0,
            0.0, 0.0
        };
        transform.backend = "analytic";
        nlp::ObjectiveSecondOrderData output_objective;
        output_objective.variable_count = 2;
        output_objective.gradient = {3.0, 5.0};
        output_objective.hessian_row_major = {
            7.0, 11.0,
            11.0, 13.0
        };
        output_objective.backend = "analytic";
        const nlp::ObjectiveSecondOrderData transformed =
            nlp::transformed_objective_second_order(transform, output_objective);
        out["transformed_lower"] =
            nlp::lower_triangle_values(transformed.hessian_row_major, transformed.variable_count);
        return out;
    });
    m.def("_native_ipopt_smoke", []() {
        py::dict out;
        const auto adapter = epcsaft::native::equilibrium_nlp::native_ipopt_adapter_info();
        const bool compiled = adapter.compiled;
        out["backend"] = "ipopt";
        out["compiled"] = compiled;
        out["available"] = compiled;
        out["adapter_available"] = adapter.adapter_available;
        out["adapter_kind"] = adapter.adapter_kind;
        out["adapter_source_available"] = true;
        out["requires_exact_gradient"] = adapter.exact_gradient_required;
        out["requires_exact_jacobian"] = adapter.exact_jacobian_required;
#ifdef EPCSAFT_IPOPT_STATUS
        out["status"] = EPCSAFT_IPOPT_STATUS;
#else
        out["status"] = adapter.status;
#endif
        return out;
    });
    m.def("_native_ipopt_quadratic_smoke", [](
        const std::string& hessian_mode,
        int iteration_history_limit,
        const std::string& linear_solver,
        double acceptable_tolerance,
        double constraint_violation_tolerance,
        double dual_infeasibility_tolerance,
        double complementarity_tolerance,
        const py::object& continuation_state
    ) {
        py::dict out;
        const auto adapter = epcsaft::native::equilibrium_nlp::native_ipopt_adapter_info();
        out["backend"] = "ipopt";
        out["compiled"] = adapter.compiled;
        out["adapter_available"] = adapter.adapter_available;
        out["adapter_kind"] = adapter.adapter_kind;
        out["problem"] = "quadratic_linear_constraint_smoke";
        out["ran"] = false;
        out["accepted"] = false;
        if (!adapter.compiled) {
            out["status"] = "ipopt_dependency_required";
            return out;
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options;
        options.max_iterations = 50;
        options.tolerance = 1.0e-10;
        options.hessian_mode = hessian_mode;
        options.iteration_history_limit = iteration_history_limit;
        options.linear_solver = linear_solver;
        options.acceptable_tolerance = acceptable_tolerance;
        options.constraint_violation_tolerance = constraint_violation_tolerance;
        options.dual_infeasibility_tolerance = dual_infeasibility_tolerance;
        options.complementarity_tolerance = complementarity_tolerance;
        apply_ipopt_continuation_state(options, continuation_state);
        const auto result = epcsaft::native::equilibrium_nlp::solve_ipopt_quadratic_smoke(options);
        out["ran"] = result.solver_ran;
        out["accepted"] = result.accepted;
        out["status"] = result.solver_status;
        out["application_status"] = result.application_status;
        out["objective"] = result.objective;
        out["variables"] = result.variables;
        out["constraints"] = result.constraints;
        out["continuation_state"] = ipopt_continuation_state_to_dict(result);
        out["iteration_history"] = ipopt_iteration_history_to_list(result.iteration_history);
        out["iteration_count"] = diagnostic_int_or(result, "iteration_count", 0);
        out["iteration_history_limit"] = diagnostic_int_or(result, "iteration_history_limit", 0);
        out["hessian_approximation"] = diagnostic_string_or(result, "hessian_approximation", "");
        out["hessian_backend"] = diagnostic_string_or(result, "hessian_backend", "");
        out["eval_h_calls"] = diagnostic_int_or(result, "eval_h_calls", 0);
        out["scaling_method"] = diagnostic_string_or(result, "scaling_method", "");
        out["variable_scaling_count"] = diagnostic_int_or(result, "variable_scaling_count", 0);
        out["constraint_scaling_count"] = diagnostic_int_or(result, "constraint_scaling_count", 0);
        out["objective_scaling"] = diagnostic_double_or(result, "objective_scaling", 1.0);
        out["acceptable_tolerance"] = diagnostic_double_or(result, "acceptable_tolerance", 0.0);
        out["constraint_violation_tolerance"] = diagnostic_double_or(result, "constraint_violation_tolerance", 0.0);
        out["dual_infeasibility_tolerance"] = diagnostic_double_or(result, "dual_infeasibility_tolerance", 0.0);
        out["complementarity_tolerance"] = diagnostic_double_or(result, "complementarity_tolerance", 0.0);
        out["variable_scaling_min"] = diagnostic_double_or(result, "variable_scaling_min", 0.0);
        out["variable_scaling_max"] = diagnostic_double_or(result, "variable_scaling_max", 0.0);
        out["constraint_scaling_min"] = diagnostic_double_or(result, "constraint_scaling_min", 0.0);
        out["constraint_scaling_max"] = diagnostic_double_or(result, "constraint_scaling_max", 0.0);
        out["linear_solver_requested"] = diagnostic_string_or(result, "linear_solver_requested", "auto");
        out["linear_solver_selected"] = diagnostic_string_or(result, "linear_solver_selected", "default");
        out["warm_start_requested"] = diagnostic_bool_or(result, "warm_start_requested", false);
        out["warm_start_used"] = diagnostic_bool_or(result, "warm_start_used", false);
        out["exact_hessian_available"] = diagnostic_bool_or(result, "exact_hessian_available", false);
        out["exact_gradient_required"] = adapter.exact_gradient_required;
        out["exact_jacobian_required"] = adapter.exact_jacobian_required;
        return out;
    },
        py::arg("hessian_mode") = "auto",
        py::arg("iteration_history_limit") = 20,
        py::arg("linear_solver") = "auto",
        py::arg("acceptable_tolerance") = 1.0e-8,
        py::arg("constraint_violation_tolerance") = 1.0e-10,
        py::arg("dual_infeasibility_tolerance") = 1.0e-10,
        py::arg("complementarity_tolerance") = 1.0e-10,
        py::arg("continuation_state") = py::none()
    );
    m.def("_native_ideal_reaction_smoke", []() {
        const double log_k = std::log(3.0);
        const std::vector<double> stoichiometry = {-1.0, 1.0};
        const std::vector<double> amounts = epcsaft::native::equilibrium_nlp::amounts_from_reaction_extents(
            {1.0, 0.0},
            1,
            stoichiometry,
            {0.75}
        );
        const std::vector<double> standard_mu_rt = epcsaft::native::equilibrium_nlp::standard_mu_rt_from_reactions(
            2,
            1,
            stoichiometry,
            {log_k}
        );
        const auto gibbs = epcsaft::native::equilibrium_nlp::evaluate_ideal_reduced_gibbs(
            amounts,
            standard_mu_rt,
            true
        );
        const auto reactions = epcsaft::native::equilibrium_nlp::evaluate_ideal_reaction_quotients(
            amounts,
            1,
            stoichiometry,
            {log_k}
        );
        py::dict out;
        out["model"] = "homogeneous_ideal_reaction";
        out["amounts"] = amounts;
        out["initial_amounts"] = std::vector<double>{1.0, 0.0};
        out["extents"] = std::vector<double>{0.75};
        out["mole_fractions"] = gibbs.mole_fractions;
        out["reduced_gibbs"] = gibbs.value;
        out["standard_mu_rt"] = standard_mu_rt;
        out["gradient"] = gibbs.gradient;
        out["hessian_row_major"] = gibbs.hessian_row_major;
        out["log_q"] = reactions.log_q;
        out["residuals"] = reactions.residuals;
        out["reaction_jacobian_row_major"] = reactions.jacobian_row_major;
        out["reaction_stationarity"] = gibbs.gradient[1] - gibbs.gradient[0];
        out["convex_kernel_scope"] = "homogeneous_ideal_reaction_validation";
        py::dict phase_residuals;
        phase_residuals["ideal_liquid"] = reactions.residuals[0];
        phase_residuals["ideal_vapor"] = reactions.residuals[0];
        out["phase_validation_residuals"] = phase_residuals;
        return out;
    });
    m.def("_native_eos_phase_block", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& amounts,
        double volume
    ) {
        if (!mixture) {
            throw ValueError("EOS phase block requires a native mixture.");
        }
        return eos_phase_block_to_dict(epcsaft::native::equilibrium_nlp::evaluate_eos_phase_block(
            mixture->args(),
            temperature,
            target_pressure,
            amounts,
            volume
        ));
    });
    m.def("_native_eos_phase_system", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        const std::vector<double>& feed_amounts,
        const std::vector<double>& charges,
        const std::vector<std::vector<double>>& association_site_fractions,
        const std::vector<double>& association_delta_row_major
    ) {
        if (!mixture) {
            throw ValueError("EOS phase system requires a native mixture.");
        }
        return eos_phase_system_to_dict(epcsaft::native::equilibrium_nlp::evaluate_eos_phase_system(
            mixture->args(),
            temperature,
            target_pressure,
            phase_amounts,
            volumes,
            feed_amounts,
            charges,
            association_site_fractions,
            association_delta_row_major
        ));
    }, py::arg("mixture"),
       py::arg("temperature"),
       py::arg("target_pressure"),
       py::arg("phase_amounts"),
       py::arg("volumes"),
       py::arg("feed_amounts"),
       py::arg("charges") = std::vector<double>{},
       py::arg("association_site_fractions") = std::vector<std::vector<double>>{},
       py::arg("association_delta_row_major") = std::vector<double>{});
    m.def("_native_association_mass_action_block", [](
        double density,
        const std::vector<double>& site_fractions,
        const std::vector<double>& site_composition,
        const std::vector<double>& delta_row_major
    ) {
        return association_mass_action_block_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_association_mass_action_block(
                density,
                site_fractions,
                site_composition,
                delta_row_major
            )
        );
    });
    m.def("_native_electrolyte_contribution_block", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double density,
        const std::vector<double>& composition,
        const std::vector<double>& amounts
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte contribution block requires a native mixture.");
        }
        return electrolyte_contribution_block_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_electrolyte_contribution_block(
                mixture->args(),
                temperature,
                density,
                composition,
                amounts
            )
        );
    });
    m.def("_native_neutral_two_phase_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        const std::vector<double>& feed_amounts
    ) {
        if (!mixture) {
            throw ValueError("Neutral two-phase EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_two_phase_eos_nlp_contract(
                mixture->args(),
                temperature,
                target_pressure,
                phase_amounts,
                volumes,
                feed_amounts
            )
        );
    });
    // AlgID: neutral_tp_flash_ipopt
    m.def("_native_neutral_tp_flash_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts
    ) {
        if (!mixture) {
            throw ValueError("Neutral TP flash EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_two_phase_eos_tp_flash_nlp_contract(
                mixture->args(),
                temperature,
                target_pressure,
                feed_amounts
            )
        );
    });
    // AlgID: neutral_lle_ipopt
    m.def("_native_neutral_lle_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts
    ) {
        if (!mixture) {
            throw ValueError("Neutral LLE EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_two_phase_eos_lle_nlp_contract(
                mixture->args(),
                temperature,
                target_pressure,
                feed_amounts
            )
        );
    });
    // AlgID: electrolyte_lle_ipopt
    m.def("_native_electrolyte_lle_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte LLE EOS NLP contract requires a native mixture.");
        }
        EquilibriumOptionsNative options;
        return neutral_two_phase_eos_nlp_contract_to_dict(
            evaluate_electrolyte_lle_liquid_root_nlp_contract_native(
                mixture,
                temperature,
                target_pressure,
                feed_amounts,
                options,
                {},
                1.0e-1
            )
        );
    });
    m.def("_native_reactive_two_phase_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants
    ) {
        if (!mixture) {
            throw ValueError("Reactive two-phase EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_reactive_two_phase_eos_nlp_contract(
                mixture->args(),
                temperature,
                target_pressure,
                phase_amounts,
                volumes,
                balance_rows,
                balance_matrix_row_major,
                total_vector,
                reaction_rows,
                reaction_stoichiometry_row_major,
                log_equilibrium_constants
            )
        );
    });
    // AlgID: reactive_lle_liquid_root_ipopt
    m.def("_native_reactive_lle_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants,
        const std::vector<int>& reaction_standard_states,
        const std::vector<double>& reaction_phase_stoichiometry
    ) {
        if (!mixture) {
            throw ValueError("Reactive LLE EOS NLP contract requires a native mixture.");
        }
        EquilibriumOptionsNative options;
        std::vector<int> standard_states = reaction_standard_states;
        if (standard_states.empty()) {
            standard_states.assign(static_cast<std::size_t>(reaction_rows), 0);
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            evaluate_reactive_phase_liquid_root_nlp_contract_native(
                mixture,
                temperature,
                target_pressure,
                feed_amounts,
                options,
                balance_matrix_row_major,
                balance_rows,
                total_vector,
                reaction_stoichiometry_row_major,
                reaction_rows,
                log_equilibrium_constants,
                standard_states,
                reaction_phase_stoichiometry,
                1.0e-8
            )
        );
    },
        py::arg("mixture"),
        py::arg("temperature"),
        py::arg("target_pressure"),
        py::arg("feed_amounts"),
        py::arg("balance_rows"),
        py::arg("balance_matrix_row_major"),
        py::arg("total_vector"),
        py::arg("reaction_rows"),
        py::arg("reaction_stoichiometry_row_major"),
        py::arg("log_equilibrium_constants"),
        py::arg("reaction_standard_states") = std::vector<int>{},
        py::arg("reaction_phase_stoichiometry") = std::vector<double>{}
    );
    // AlgID: reactive_electrolyte_lle_liquid_root_ipopt
    m.def("_native_reactive_electrolyte_lle_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants
    ) {
        if (!mixture) {
            throw ValueError("Reactive electrolyte LLE EOS NLP contract requires a native mixture.");
        }
        EquilibriumOptionsNative options;
        std::vector<int> reaction_standard_states(static_cast<std::size_t>(reaction_rows), 0);
        return neutral_two_phase_eos_nlp_contract_to_dict(
            evaluate_reactive_phase_liquid_root_nlp_contract_native(
                mixture,
                temperature,
                target_pressure,
                feed_amounts,
                options,
                balance_matrix_row_major,
                balance_rows,
                total_vector,
                reaction_stoichiometry_row_major,
                reaction_rows,
                log_equilibrium_constants,
                reaction_standard_states,
                {},
                1.0e-8
            )
        );
    });
    m.def("_native_reactive_two_phase_eos_postsolve", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants,
        double conserved_balance_tolerance,
        double pressure_tolerance,
        double reaction_stationarity_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Reactive two-phase EOS postsolve requires a native mixture.");
        }
        return reactive_two_phase_eos_postsolve_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_reactive_two_phase_eos_postsolve(
                mixture->args(),
                temperature,
                target_pressure,
                phase_amounts,
                volumes,
                balance_rows,
                balance_matrix_row_major,
                total_vector,
                reaction_rows,
                reaction_stoichiometry_row_major,
                log_equilibrium_constants,
                conserved_balance_tolerance,
                pressure_tolerance,
                reaction_stationarity_tolerance,
                phase_distance_tolerance,
                {}
            )
        );
    });
    m.def("_native_reactive_electrolyte_two_phase_eos_postsolve", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants,
        double conserved_balance_tolerance,
        double pressure_tolerance,
        double reaction_stationarity_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Reactive electrolyte two-phase EOS postsolve requires a native mixture.");
        }
        return reactive_two_phase_eos_postsolve_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_reactive_two_phase_eos_postsolve(
                mixture->args(),
                temperature,
                target_pressure,
                phase_amounts,
                volumes,
                balance_rows,
                balance_matrix_row_major,
                total_vector,
                reaction_rows,
                reaction_stoichiometry_row_major,
                log_equilibrium_constants,
                conserved_balance_tolerance,
                pressure_tolerance,
                reaction_stationarity_tolerance,
                phase_distance_tolerance,
                mixture->args().z
            )
        );
    });
    // AlgID: bubble_dew_ipopt
    m.def("_native_neutral_bubble_p_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        const std::vector<double>& liquid_composition
    ) {
        if (!mixture) {
            throw ValueError("Neutral bubble pressure EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_bubble_p_eos_nlp_contract(
                mixture->args(),
                temperature,
                liquid_composition
            )
        );
    });
    m.def("_native_electrolyte_bubble_p_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        const std::vector<double>& liquid_composition
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte bubble pressure EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_electrolyte_bubble_p_eos_nlp_contract(
                mixture->args(),
                temperature,
                liquid_composition
            )
        );
    });
    m.def("_native_neutral_dew_p_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        const std::vector<double>& vapor_composition
    ) {
        if (!mixture) {
            throw ValueError("Neutral dew pressure EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_dew_p_eos_nlp_contract(
                mixture->args(),
                temperature,
                vapor_composition
            )
        );
    });
    m.def("_native_neutral_bubble_t_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double target_pressure,
        const std::vector<double>& liquid_composition
    ) {
        if (!mixture) {
            throw ValueError("Neutral bubble temperature EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_bubble_t_eos_nlp_contract(
                mixture->args(),
                target_pressure,
                liquid_composition
            )
        );
    });
    m.def("_native_neutral_dew_t_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double target_pressure,
        const std::vector<double>& vapor_composition
    ) {
        if (!mixture) {
            throw ValueError("Neutral dew temperature EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_dew_t_eos_nlp_contract(
                mixture->args(),
                target_pressure,
                vapor_composition
            )
        );
    });
    // AlgID: stability_tpd_ipopt
    m.def("_native_neutral_stability_tpd_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double pressure,
        const std::vector<double>& feed_composition,
        const std::string& parent_phase,
        const std::string& trial_phase
    ) {
        if (!mixture) {
            throw ValueError("Neutral stability TPD NLP contract requires a native mixture.");
        }
        return stability_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_stability_tpd_nlp_contract(
                mixture->args(),
                temperature,
                pressure,
                feed_composition,
                stability_phase_token_to_int(parent_phase),
                stability_phase_token_to_int(trial_phase)
            )
        );
    });
    m.def("_native_electrolyte_stability_tpd_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double pressure,
        const std::vector<double>& feed_composition
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte stability TPD NLP contract requires a native mixture.");
        }
        return stability_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_electrolyte_stability_tpd_nlp_contract(
                mixture->args(),
                temperature,
                pressure,
                feed_composition
            )
        );
    });
    m.def("_native_neutral_two_phase_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        const std::vector<double>& feed_amounts,
        int max_iterations,
        double tolerance,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double material_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Neutral two-phase EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(max_iterations, tolerance, 0.0, hessian_mode, iteration_history_limit);
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return neutral_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_neutral_two_phase_eos_route(
                mixture->args(),
                temperature,
                target_pressure,
                phase_amounts,
                volumes,
                feed_amounts,
                options,
                material_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_reactive_two_phase_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double conserved_balance_tolerance,
        double pressure_tolerance,
        double reaction_stationarity_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Reactive two-phase EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return reactive_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_reactive_two_phase_eos_route(
                mixture->args(),
                temperature,
                target_pressure,
                phase_amounts,
                volumes,
                balance_rows,
                balance_matrix_row_major,
                total_vector,
                reaction_rows,
                reaction_stoichiometry_row_major,
                log_equilibrium_constants,
                options,
                conserved_balance_tolerance,
                pressure_tolerance,
                reaction_stationarity_tolerance,
                phase_distance_tolerance,
                {}
            )
        );
    });
    m.def("_native_reactive_lle_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double conserved_balance_tolerance,
        double phase_equilibrium_tolerance,
        double reaction_stationarity_tolerance,
        double phase_distance_tolerance,
        double min_composition,
        const std::vector<int>& reaction_standard_states,
        const std::vector<double>& reaction_phase_stoichiometry,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Reactive LLE EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        EquilibriumOptionsNative equilibrium_options;
        equilibrium_options.min_composition = min_composition;
        return reactive_two_phase_eos_route_result_to_dict(
            solve_reactive_phase_liquid_root_route_native(
                mixture,
                temperature,
                target_pressure,
                feed_amounts,
                equilibrium_options,
                balance_matrix_row_major,
                balance_rows,
                total_vector,
                reaction_stoichiometry_row_major,
                reaction_rows,
                log_equilibrium_constants,
                reaction_standard_states,
                reaction_phase_stoichiometry,
                options,
                conserved_balance_tolerance,
                reaction_stationarity_tolerance,
                phase_equilibrium_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_reactive_electrolyte_lle_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double conserved_balance_tolerance,
        double phase_equilibrium_tolerance,
        double reaction_stationarity_tolerance,
        double phase_distance_tolerance,
        double min_composition,
        const std::vector<int>& reaction_standard_states,
        const std::vector<double>& reaction_phase_stoichiometry,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Reactive electrolyte LLE EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        EquilibriumOptionsNative equilibrium_options;
        equilibrium_options.min_composition = min_composition;
        return reactive_two_phase_eos_route_result_to_dict(
            solve_reactive_phase_liquid_root_route_native(
                mixture,
                temperature,
                target_pressure,
                feed_amounts,
                equilibrium_options,
                balance_matrix_row_major,
                balance_rows,
                total_vector,
                reaction_stoichiometry_row_major,
                reaction_rows,
                log_equilibrium_constants,
                reaction_standard_states,
                reaction_phase_stoichiometry,
                options,
                conserved_balance_tolerance,
                reaction_stationarity_tolerance,
                phase_equilibrium_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_reactive_electrolyte_lle_phase_model_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        const std::shared_ptr<ePCSAFTMixtureNative>& phase1_mixture,
        const std::shared_ptr<ePCSAFTMixtureNative>& phase2_mixture,
        const std::vector<int>& phase1_global_indices,
        const std::vector<int>& phase2_global_indices,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int balance_rows,
        const std::vector<double>& balance_matrix_row_major,
        const std::vector<double>& total_vector,
        int reaction_rows,
        const std::vector<double>& reaction_stoichiometry_row_major,
        const std::vector<double>& log_equilibrium_constants,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double conserved_balance_tolerance,
        double phase_equilibrium_tolerance,
        double reaction_stationarity_tolerance,
        double phase_distance_tolerance,
        double min_composition,
        const std::vector<int>& reaction_standard_states,
        const std::vector<double>& reaction_phase_stoichiometry,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture || !phase1_mixture || !phase2_mixture) {
            throw ValueError("Reactive electrolyte LLE phase_models route requires global, aqueous, and organic native mixtures.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        EquilibriumOptionsNative equilibrium_options;
        equilibrium_options.min_composition = min_composition;
        return reactive_two_phase_eos_route_result_to_dict(
            solve_reactive_phase_liquid_root_phase_model_route_native(
                mixture,
                phase1_mixture,
                phase2_mixture,
                phase1_global_indices,
                phase2_global_indices,
                temperature,
                target_pressure,
                feed_amounts,
                equilibrium_options,
                balance_matrix_row_major,
                balance_rows,
                total_vector,
                reaction_stoichiometry_row_major,
                reaction_rows,
                log_equilibrium_constants,
                reaction_standard_states,
                reaction_phase_stoichiometry,
                options,
                conserved_balance_tolerance,
                reaction_stationarity_tolerance,
                phase_equilibrium_tolerance,
                phase_distance_tolerance
            )
        );
    });
    // AlgID: neutral_tp_flash_ipopt
    m.def("_native_neutral_tp_flash_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double material_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Neutral TP flash EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return neutral_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_neutral_two_phase_eos_tp_flash_route(
                mixture->args(),
                temperature,
                target_pressure,
                feed_amounts,
                options,
                material_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    // AlgID: neutral_lle_ipopt
    m.def("_native_neutral_lle_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double material_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Neutral LLE EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return neutral_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_neutral_two_phase_eos_lle_route(
                mixture->args(),
                temperature,
                target_pressure,
                feed_amounts,
                options,
                material_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    // AlgID: electrolyte_lle_ipopt
    m.def("_native_electrolyte_lle_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double material_tolerance,
        double pressure_tolerance,
        double charge_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte LLE EOS route result requires a native mixture.");
        }
        (void)pressure_tolerance;
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        EquilibriumOptionsNative equilibrium_options;
        return neutral_two_phase_eos_route_result_to_dict(
            solve_electrolyte_lle_liquid_root_route_native(
                mixture,
                temperature,
                target_pressure,
                feed_amounts,
                equilibrium_options,
                {},
                options,
                material_tolerance,
                charge_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    // AlgID: bubble_dew_ipopt
    m.def("_native_neutral_bubble_p_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        const std::vector<double>& liquid_composition,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double phase_total_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Neutral bubble pressure EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return neutral_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_neutral_bubble_p_eos_route(
                mixture->args(),
                temperature,
                liquid_composition,
                options,
                phase_total_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_electrolyte_bubble_p_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        const std::vector<double>& liquid_composition,
        int max_iterations,
        double tolerance,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double phase_total_tolerance,
        double pressure_tolerance,
        double charge_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte bubble pressure EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(max_iterations, tolerance, 0.0, hessian_mode, iteration_history_limit);
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return neutral_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_electrolyte_bubble_p_eos_route(
                mixture->args(),
                temperature,
                liquid_composition,
                options,
                phase_total_tolerance,
                pressure_tolerance,
                charge_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_neutral_dew_p_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        const std::vector<double>& vapor_composition,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double phase_total_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Neutral dew pressure EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return neutral_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_neutral_dew_p_eos_route(
                mixture->args(),
                temperature,
                vapor_composition,
                options,
                phase_total_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_neutral_bubble_t_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double target_pressure,
        const std::vector<double>& liquid_composition,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double phase_total_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Neutral bubble temperature EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return neutral_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_neutral_bubble_t_eos_route(
                mixture->args(),
                target_pressure,
                liquid_composition,
                options,
                phase_total_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_neutral_dew_t_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double target_pressure,
        const std::vector<double>& vapor_composition,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double phase_total_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Neutral dew temperature EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return neutral_two_phase_eos_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_neutral_dew_t_eos_route(
                mixture->args(),
                target_pressure,
                vapor_composition,
                options,
                phase_total_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    // AlgID: stability_tpd_ipopt
    m.def("_native_neutral_stability_tpd_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double pressure,
        const std::vector<double>& feed_composition,
        const std::string& parent_phase,
        const std::string& trial_phase,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double stability_tolerance,
        const std::vector<double>& trial_initial_composition,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Neutral stability TPD route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return stability_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_neutral_stability_tpd_route(
                mixture->args(),
                temperature,
                pressure,
                feed_composition,
                stability_phase_token_to_int(parent_phase),
                stability_phase_token_to_int(trial_phase),
                options,
                stability_tolerance,
                trial_initial_composition
            )
        );
    },
        py::arg("mixture"),
        py::arg("temperature"),
        py::arg("pressure"),
        py::arg("feed_composition"),
        py::arg("parent_phase"),
        py::arg("trial_phase"),
        py::arg("max_iterations"),
        py::arg("tolerance"),
        py::arg("timeout_seconds"),
        py::arg("hessian_mode"),
        py::arg("iteration_history_limit"),
        py::arg("stability_tolerance"),
        py::arg("trial_initial_composition") = std::vector<double>{},
        py::arg("continuation_state") = py::none()
    );
    m.def("_native_electrolyte_stability_tpd_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double pressure,
        const std::vector<double>& feed_composition,
        int max_iterations,
        double tolerance,
        double timeout_seconds,
        const std::string& hessian_mode,
        int iteration_history_limit,
        double stability_tolerance,
        const std::vector<double>& trial_initial_composition,
        const py::object& continuation_state,
        const py::kwargs& kwargs
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte stability TPD route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options =
            ipopt_solve_options_from_scalars(
                max_iterations,
                tolerance,
                timeout_seconds,
                hessian_mode,
                iteration_history_limit
            );
        apply_ipopt_control_kwargs(options, kwargs);
        apply_ipopt_continuation_state(options, continuation_state);
        return stability_route_result_to_dict(
            epcsaft::native::equilibrium_nlp::solve_electrolyte_stability_tpd_route(
                mixture->args(),
                temperature,
                pressure,
                feed_composition,
                options,
                stability_tolerance,
                trial_initial_composition
            )
        );
    },
        py::arg("mixture"),
        py::arg("temperature"),
        py::arg("pressure"),
        py::arg("feed_composition"),
        py::arg("max_iterations"),
        py::arg("tolerance"),
        py::arg("timeout_seconds"),
        py::arg("hessian_mode"),
        py::arg("iteration_history_limit"),
        py::arg("stability_tolerance"),
        py::arg("trial_initial_composition") = std::vector<double>{},
        py::arg("continuation_state") = py::none()
    );
    m.def("_native_neutral_two_phase_eos_postsolve", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        const std::vector<double>& feed_amounts,
        double material_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Neutral two-phase EOS postsolve requires a native mixture.");
        }
        return neutral_two_phase_eos_postsolve_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_neutral_two_phase_eos_postsolve(
                mixture->args(),
                temperature,
                target_pressure,
                phase_amounts,
                volumes,
                feed_amounts,
                material_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_neutral_two_phase_eos_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        const std::vector<double>& feed_amounts,
        double material_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Neutral two-phase EOS result builder requires a native mixture.");
        }
        return neutral_two_phase_eos_result_payload_to_dict(
            epcsaft::native::equilibrium_nlp::build_neutral_two_phase_eos_result(
                mixture->args(),
                temperature,
                target_pressure,
                phase_amounts,
                volumes,
                feed_amounts,
                material_tolerance,
                pressure_tolerance,
                chemical_potential_tolerance,
                phase_distance_tolerance
            )
        );
    });
    m.def("_native_cppad_eos_contributions", [](double t, double rho, const std::vector<double>& x, const add_args& args) {
        return cppad_smoke_to_dict(cppad_eos_contribution_derivatives_cpp(t, rho, x, args));
    });
    m.def("_native_cppad_pressure_density", [](double t, double rho, const std::vector<double>& x, const add_args& args) {
        return cppad_smoke_to_dict(cppad_pressure_density_derivative_cpp(t, rho, x, args));
    });
    m.def("_native_phase_state_ln_fugacity_composition_sensitivity", [](
        double t,
        double p,
        const std::vector<double>& x,
        int phase,
        const add_args& args
    ) {
        return phase_state_sensitivity_to_dict(
            phase_state_ln_fugacity_composition_sensitivity_cpp(t, p, x, phase, args)
        );
    });
    m.def("_native_cppad_pure_neutral_parameters", [](double t, double rho, const add_args& args) {
        return cppad_smoke_to_dict(cppad_pure_neutral_parameter_derivatives_cpp(t, rho, args));
    });
    m.def("_native_cppad_neutral_binary_kij_properties", [](double t, double rho, const std::vector<double>& x, const add_args& args) {
        NeutralBinaryKijPhaseDerivatives forward = neutral_binary_pair_parameter_phase_derivatives_cpp(t, rho, x, args, 1, "k_ij");
        NeutralBinaryKijPhaseDerivatives reverse = neutral_binary_pair_parameter_phase_derivatives_cpp(t, rho, x, args, 2, "k_ij");
        return neutral_binary_kij_property_derivatives_to_dict(forward, reverse);
    });
    m.def("_native_cppad_neutral_binary_pair_properties", [](double t, double rho, const std::vector<double>& x, const add_args& args) {
        NeutralBinaryKijPhaseDerivatives kij_forward = neutral_binary_pair_parameter_phase_derivatives_cpp(t, rho, x, args, 1, "k_ij");
        NeutralBinaryKijPhaseDerivatives kij_reverse = neutral_binary_pair_parameter_phase_derivatives_cpp(t, rho, x, args, 2, "k_ij");
        if (args.l_ij.size() != 4) {
            return neutral_binary_pair_property_derivatives_to_dict(kij_forward, kij_reverse, nullptr, nullptr);
        }
        NeutralBinaryKijPhaseDerivatives lij_forward = neutral_binary_pair_parameter_phase_derivatives_cpp(t, rho, x, args, 1, "l_ij");
        NeutralBinaryKijPhaseDerivatives lij_reverse = neutral_binary_pair_parameter_phase_derivatives_cpp(t, rho, x, args, 2, "l_ij");
        return neutral_binary_pair_property_derivatives_to_dict(kij_forward, kij_reverse, &lij_forward, &lij_reverse);
    });

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
        .def_readonly("derivative_available", &CompositionContributionResult::derivative_available);

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
        .def("density_warm_start_rejections", &ePCSAFTMixtureNative::density_warm_start_rejections)
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
        .def("ln_fugacity_coefficient", &ePCSAFTStateNative::ln_fugacity_coefficient)
        .def("fugacity_coefficient", &ePCSAFTStateNative::fugacity_coefficient)
        .def("fugacity_coefficient_result", &ePCSAFTStateNative::fugacity_coefficient_result)
        .def("born_ssmds_liquid_derivatives", [](ePCSAFTStateNative& state) {
            return born_ssmds_derivative_to_dict(state.born_ssmds_liquid_derivatives());
        })
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

    // AlgID: pure_neutral_ceres_regression
    m.def("_fit_pure_neutral_native_ceres", &fit_pure_neutral_native_ceres_binding);
    m.def("_fit_pure_neutral_native_debug", &evaluate_pure_neutral_objective_debug_binding);
    // AlgID: native_ceres_regression_adapter
    // AlgID: pure_ion_ceres_regression
    // AlgID: binary_kij_ceres_regression
    m.def("_fit_generic_native_ceres", &fit_generic_native_ceres_binding);
    m.def("_evaluate_generic_native_debug", &evaluate_generic_native_debug_binding);
    m.def("_evaluate_electrolyte_lle_residual_native", &evaluate_electrolyte_lle_residual_native_binding);
    m.def("_evaluate_reactive_phase_equilibrium_residual_native", &evaluate_reactive_phase_equilibrium_residual_native_binding);
    // AlgID: ideal_speciation_ipopt
    // AlgID: nonideal_speciation_ipopt
    m.def("_solve_chemical_equilibrium_native", &solve_chemical_equilibrium_native_binding);
    m.def("_evaluate_chemical_equilibrium_residual_native", &evaluate_chemical_equilibrium_residual_native_binding);
}
