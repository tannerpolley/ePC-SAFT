#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <string>
#include <vector>

#include "epcsaft_equilibrium.h"

namespace py = pybind11;

namespace {

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
    return out;
}

py::dict regression_debug_to_dict(const PureNeutralRegressionDebugResult& result) {
    py::dict out;
    out["objective"] = result.objective;
    out["gradient"] = result.gradient;
    out["residuals"] = result.residuals;
    out["jacobian_row_major"] = result.jacobian_row_major;
    out["jacobian_shape"] = py::make_tuple(result.jacobian_rows, result.jacobian_cols);
    out["density_raw_residuals"] = result.density_raw_residuals;
    out["pure_vle_raw_residuals"] = result.pure_vle_raw_residuals;
    out["residual_evaluations"] = result.residual_evaluations;
    out["density_solves"] = result.density_solves;
    out["fused_state_evaluations"] = result.fused_state_evaluations;
    out["callback_wall_time_s"] = result.callback_wall_time_s;
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
    return out;
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
    return options;
}

std::vector<std::string> string_vector_from_request(const py::dict& request, const char* key, std::vector<std::string> fallback) {
    if (!request.contains(key) || request[key].is_none()) {
        return fallback;
    }
    return request[key].cast<std::vector<std::string>>();
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

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "pybind11 native backend for epcsaft";

    py::register_exception<ValueError>(m, "NativeValueError");
    py::register_exception<SolutionError>(m, "NativeSolutionError");

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
        .def_readonly("z", &CompositionContributionResult::z);

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
        .def("density_warm_start_fallbacks", &ePCSAFTMixtureNative::density_warm_start_fallbacks);

    py::class_<ePCSAFTStateNative, std::shared_ptr<ePCSAFTStateNative>>(m, "NativeState")
        .def(py::init<std::shared_ptr<ePCSAFTMixtureNative>, double, std::vector<double>, int, bool, double, bool, double>())
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
    m.def("_solve_equilibrium_native", &solve_equilibrium_native_binding);
}
