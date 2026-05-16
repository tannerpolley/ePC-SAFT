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

epcsaft::native::cppad_support::CppADDerivativeResult cppad_eos_contribution_derivatives_cpp(
    double t,
    double rho,
    const std::vector<double>& x,
    const add_args& cppargs
);
epcsaft::native::cppad_support::CppADDerivativeResult cppad_pressure_density_derivative_cpp(
    double t,
    double rho
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

py::dict cppad_smoke_to_dict(const epcsaft::native::cppad_support::CppADDerivativeResult& result) {
    py::dict out;
    out["cppad_compiled"] = epcsaft::native::cppad_support::cppad_compiled();
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
    out["shape"] = py::make_tuple(result.rows, result.cols);
    out["composition"] = result.composition;
    out["ln_fugacity"] = result.ln_fugacity;
    out["density_composition_derivative"] = result.density_composition_derivative;
    out["pressure_composition_fixed_density_derivative"] = result.pressure_composition_fixed_density_derivative;
    out["ln_fugacity_density_derivative"] = result.ln_fugacity_density_derivative;
    out["fixed_density_jacobian_row_major"] = result.fixed_density_jacobian_row_major;
    out["jacobian_row_major"] = result.jacobian_row_major;
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
    out["gradient"] = result.gradient;
    out["constraints"] = result.constraints;
    out["phase_charge_residuals"] = result.phase_charge_residuals;
    out["phase_association_residuals"] = result.phase_association_residuals;
    out["constraint_jacobian_backend"] = result.constraint_jacobian_backend;
    out["constraint_jacobian_shape"] =
        py::make_tuple(result.constraint_jacobian_rows, result.constraint_jacobian_cols);
    out["constraint_jacobian_row_major"] = result.constraint_jacobian_row_major;
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
    out["phase_count"] = result.phase_count;
    out["species_count"] = result.species_count;
    out["variable_count"] = result.variable_count;
    out["constraint_count"] = result.constraint_count;
    out["jacobian_nonzero_count"] = result.jacobian_nonzero_count;
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
    out["constraints"] = result.constraints;
    out["phase_amount_totals"] = result.phase_amount_totals;
    out["phase_volumes"] = result.phase_volumes;
    out["phase_compositions"] = result.phase_compositions;
    return out;
}

py::dict neutral_two_phase_eos_route_result_to_dict(
    const epcsaft::native::equilibrium_nlp::NeutralTwoPhaseEosRouteResult& result
) {
    py::dict out;
    out["backend"] = result.backend;
    out["compiled"] = result.compiled;
    out["adapter_available"] = result.adapter_available;
    out["adapter_kind"] = result.adapter_kind;
    out["problem_name"] = result.problem_name;
    out["derivative_backend"] = result.derivative_backend;
    out["ran"] = result.ran;
    out["solver_accepted"] = result.solver_accepted;
    out["accepted"] = result.accepted;
    out["exact_gradient_required"] = result.exact_gradient_required;
    out["exact_jacobian_required"] = result.exact_jacobian_required;
    out["status"] = result.status;
    out["solver_status"] = result.solver_status;
    out["application_status"] = result.application_status;
    out["hessian_strategy"] = result.hessian_strategy;
    out["objective"] = result.objective;
    out["variables"] = result.variables;
    out["constraints"] = result.constraints;
    out["phase_amounts"] = result.phase_amounts;
    out["phase_volumes"] = result.phase_volumes;
    out["postsolve"] = neutral_two_phase_eos_postsolve_to_dict(result.postsolve);
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
    if (input.contains("phase")) {
        options.phase = input["phase"].cast<std::string>();
    }
    if (input.contains("activity_output")) {
        options.activity_output = input["activity_output"].cast<std::string>();
    }
    return options;
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
        out["hessian_strategy"] = adapter.hessian_strategy;
        out["requires_exact_gradient"] = adapter.exact_gradient_required;
        out["requires_exact_jacobian"] = adapter.exact_jacobian_required;
        out["requires_exact_hessian"] = adapter.exact_hessian_required;
#ifdef EPCSAFT_IPOPT_STATUS
        out["status"] = EPCSAFT_IPOPT_STATUS;
#else
        out["status"] = adapter.status;
#endif
        return out;
    });
    m.def("_native_ipopt_quadratic_smoke", []() {
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
            out["status"] = "requires_ipopt_build";
            return out;
        }
        const auto result = epcsaft::native::equilibrium_nlp::solve_ipopt_quadratic_smoke();
        out["ran"] = result.solver_ran;
        out["accepted"] = result.accepted;
        out["status"] = result.solver_status;
        out["application_status"] = result.application_status;
        out["objective"] = result.objective;
        out["variables"] = result.variables;
        out["constraints"] = result.constraints;
        out["hessian_strategy"] = result.hessian_strategy;
        out["exact_gradient_required"] = adapter.exact_gradient_required;
        out["exact_jacobian_required"] = adapter.exact_jacobian_required;
        return out;
    });
    m.def("_native_ideal_reaction_smoke", []() {
        const double log_k = std::log(3.0);
        const std::vector<double> stoichiometry = {-1.0, 1.0};
        const std::vector<double> amounts = epcsaft::native::equilibrium_nlp::amounts_from_reaction_extents(
            {1.0, 0.0},
            1,
            stoichiometry,
            {0.75}
        );
        const std::vector<double> standard_mu_rt = {0.0, -log_k};
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
    m.def("_native_electrolyte_lle_eos_nlp_contract", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte LLE EOS NLP contract requires a native mixture.");
        }
        return neutral_two_phase_eos_nlp_contract_to_dict(
            epcsaft::native::equilibrium_nlp::evaluate_electrolyte_lle_eos_nlp_contract(
                mixture->args(),
                temperature,
                target_pressure,
                feed_amounts
            )
        );
    });
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
    m.def("_native_neutral_two_phase_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<std::vector<double>>& phase_amounts,
        const std::vector<double>& volumes,
        const std::vector<double>& feed_amounts,
        int max_iterations,
        double tolerance,
        double material_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Neutral two-phase EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        options.acceptable_tolerance = std::max(tolerance * 100.0, 1.0e-10);
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
    m.def("_native_neutral_tp_flash_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int max_iterations,
        double tolerance,
        double material_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Neutral TP flash EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        options.acceptable_tolerance = std::max(tolerance * 100.0, 1.0e-10);
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
    m.def("_native_neutral_lle_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        double target_pressure,
        const std::vector<double>& feed_amounts,
        int max_iterations,
        double tolerance,
        double material_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Neutral LLE EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        options.acceptable_tolerance = std::max(tolerance * 100.0, 1.0e-10);
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
    m.def("_native_neutral_bubble_p_eos_route_result", [](
        const std::shared_ptr<ePCSAFTMixtureNative>& mixture,
        double temperature,
        const std::vector<double>& liquid_composition,
        int max_iterations,
        double tolerance,
        double phase_total_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Neutral bubble pressure EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        options.acceptable_tolerance = std::max(tolerance * 100.0, 1.0e-10);
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
        double phase_total_tolerance,
        double pressure_tolerance,
        double charge_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Electrolyte bubble pressure EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        options.acceptable_tolerance = std::max(tolerance * 100.0, 1.0e-10);
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
        double phase_total_tolerance,
        double pressure_tolerance,
        double chemical_potential_tolerance,
        double phase_distance_tolerance
    ) {
        if (!mixture) {
            throw ValueError("Neutral dew pressure EOS route result requires a native mixture.");
        }
        epcsaft::native::equilibrium_nlp::IpoptSolveOptions options;
        options.max_iterations = max_iterations;
        options.tolerance = tolerance;
        options.acceptable_tolerance = std::max(tolerance * 100.0, 1.0e-10);
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
        (void)x;
        (void)args;
        return cppad_smoke_to_dict(cppad_pressure_density_derivative_cpp(t, rho));
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

    m.def("_fit_pure_neutral_native_ceres", &fit_pure_neutral_native_ceres_binding);
    m.def("_fit_pure_neutral_native_debug", &evaluate_pure_neutral_objective_debug_binding);
    m.def("_fit_generic_native_ceres", &fit_generic_native_ceres_binding);
    m.def("_evaluate_generic_native_debug", &evaluate_generic_native_debug_binding);
    m.def("_evaluate_electrolyte_lle_residual_native", &evaluate_electrolyte_lle_residual_native_binding);
    m.def("_evaluate_reactive_phase_equilibrium_residual_native", &evaluate_reactive_phase_equilibrium_residual_native_binding);
    m.def("_solve_chemical_equilibrium_native", &solve_chemical_equilibrium_native_binding);
    m.def("_evaluate_chemical_equilibrium_residual_native", &evaluate_chemical_equilibrium_residual_native_binding);
}
