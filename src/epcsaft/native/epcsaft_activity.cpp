#include "epcsaft_core_internal.h"

using namespace thermo_detail;

namespace miac_detail {

// EqID: gamma_sym
// EqID: lngamma_sym
vector<double> miac_gamma_vector_cpp(double t, double rho, const vector<double>& x, const add_args& cppargs)
{
    add_args args = cppargs;
    const int ncomp = static_cast<int>(x.size());
    if (args.z.empty() || std::all_of(args.z.begin(), args.z.end(), [](double v) { return std::abs(v) <= 1e-12; })) {
        throw ValueError("miac requires ionic species (non-zero z).");
    }
    if (args.mw.size() != x.size()) {
        throw ValueError("miac requires params['MW'] to be present and aligned with x.");
    }
    ChargeGroups groups = collect_charge_groups(args, x.size());
    if (groups.cations.empty() || groups.anions.empty()) {
        throw ValueError("miac requires at least one cation and one anion.");
    }
    if (groups.solvents.empty()) {
        throw ValueError("miac requires a neutral solvent reference.");
    }

    vector<double> fugcoef = fugcoef_cpp(t, rho, x, args);
    double p = p_cpp(t, rho, x, args);

    const double eps = 1e-12;
    vector<double> x_inf(ncomp, eps);
    vector<double> solvent_ref(groups.solvents.size(), 0.0);
    double solvent_sum = 0.0;
    for (size_t k = 0; k < groups.solvents.size(); ++k) {
        solvent_ref[k] = x[groups.solvents[k]];
        solvent_sum += solvent_ref[k];
    }
    if (solvent_sum <= 0.0) {
        throw ValueError("miac requires a positive solvent fraction.");
    }
    for (size_t k = 0; k < groups.solvents.size(); ++k) {
        x_inf[groups.solvents[k]] = solvent_ref[k] / solvent_sum;
    }
    double solvent_budget = std::max(1.0 - eps * static_cast<double>(ncomp - groups.solvents.size()), eps * static_cast<double>(groups.solvents.size()));
    for (size_t k = 0; k < groups.solvents.size(); ++k) {
        x_inf[groups.solvents[k]] *= solvent_budget;
    }
    double x_inf_sum = 0.0;
    for (double xi : x_inf) {
        x_inf_sum += xi;
    }
    for (double& xi : x_inf) {
        xi /= x_inf_sum;
    }

    double rho_inf = den_cpp(t, p, x_inf, 0, args);
    vector<double> fugcoef_inf = fugcoef_cpp(t, rho_inf, x_inf, args);
    vector<double> gamma_i(ncomp, 1.0);
    for (int i = 0; i < ncomp; ++i) {
        gamma_i[i] = fugcoef[i] / fugcoef_inf[i];
    }
    return gamma_i;
}

// EqID: delta_g_solv_inf_x
// EqID: delta_g_transfer_inf
vector<double> gsolv_values_cpp(double t, double rho, const vector<double>& x, const add_args& cppargs)
{
    add_args args = cppargs;
    const int ncomp = static_cast<int>(x.size());
    if (args.z.empty() || std::all_of(args.z.begin(), args.z.end(), [](double v) { return std::abs(v) <= 1e-12; })) {
        throw ValueError("gsolv requires ionic species in params['z'].");
    }
    if (args.mw.size() != x.size()) {
        throw ValueError("gsolv requires params['MW'] to be present and aligned with x.");
    }
    ChargeGroups groups = collect_charge_groups(args, x.size());
    if (groups.cations.empty() && groups.anions.empty()) {
        throw ValueError("gsolv requires ionic species in params['z'].");
    }
    if (groups.solvents.empty()) {
        throw ValueError("gsolv requires at least one solvent species (z=0).");
    }

    vector<double> x_ref = x;
    vector<int> idx_ion = groups.cations;
    idx_ion.insert(idx_ion.end(), groups.anions.begin(), groups.anions.end());
    for (int i : idx_ion) {
        x_ref[i] = 0.0;
    }
    double solv_sum = 0.0;
    for (int i : groups.solvents) {
        solv_sum += x_ref[i];
    }
    if (solv_sum > 0.0) {
        for (int i : groups.solvents) {
            x_ref[i] /= solv_sum;
        }
    }
    else {
        double equal = 1.0 / static_cast<double>(groups.solvents.size());
        for (int i : groups.solvents) {
            x_ref[i] = equal;
        }
    }

    double p = p_cpp(t, rho, x_ref, args);
    int phase = (rho < 900.0) ? 1 : 0;
    vector<double> result(ncomp, 0.0);
    const double eps = 1e-12;
    for (int i : idx_ion) {
        vector<double> x_inf = x_ref;
        x_inf[i] = eps;
        double sum_inf = 0.0;
        for (double xi : x_inf) {
            sum_inf += xi;
        }
        for (double& xi : x_inf) {
            xi /= sum_inf;
        }
        double rho_inf = den_cpp(t, p, x_inf, phase, args);
        vector<double> lnfug_inf = lnfug_cpp(t, rho_inf, x_inf, args);
        result[i] = 8.31446261815324 * t * lnfug_inf[i];
    }
    return result;
}

int resolve_solvent_index_cpp(
    const vector<int>& solvent_indices,
    const vector<double>& x,
    bool has_solvent_override,
    int solvent_override_index
)
{
    if (solvent_indices.empty()) {
        throw ValueError("activity_coefficient requires at least one neutral solvent species.");
    }
    if (!has_solvent_override || solvent_override_index < 0) {
        for (int idx : solvent_indices) {
            if (idx >= 0 && idx < static_cast<int>(x.size()) && x[idx] > 0.0) {
                return idx;
            }
        }
        return solvent_indices.front();
    }
    if (std::find(solvent_indices.begin(), solvent_indices.end(), solvent_override_index) == solvent_indices.end()) {
        throw ValueError("solvent_override_index must reference a neutral solvent species.");
    }
    return solvent_override_index;
}

double normalize_mw_cpp(double mw)
{
    if (mw > 1.0) {
        mw /= 1000.0;
    }
    return mw;
}

double solvent_pool_mix_mw_cpp(const vector<double>& x, const add_args& args, const vector<int>& solvent_pool)
{
    if (solvent_pool.empty()) {
        throw ValueError("activity_coefficient requires a non-empty solvent pool.");
    }
    vector<double> mass_neutral(solvent_pool.size(), 0.0);
    double mass_neutral_sum = 0.0;
    for (size_t k = 0; k < solvent_pool.size(); ++k) {
        int idx = solvent_pool[k];
        double mw = normalize_mw_cpp(args.mw[idx]);
        if (mw <= 0.0) {
            throw ValueError("Solvent molecular weight must be positive.");
        }
        mass_neutral[k] = x[idx] * mw;
        mass_neutral_sum += mass_neutral[k];
    }
    if (mass_neutral_sum <= 0.0) {
        throw ValueError("Solvent mass is zero; check solvent mole fraction and MW.");
    }
    double mw_mix_inv = 0.0;
    for (size_t k = 0; k < solvent_pool.size(); ++k) {
        int idx = solvent_pool[k];
        double mw = normalize_mw_cpp(args.mw[idx]);
        double w_sf = mass_neutral[k] / mass_neutral_sum;
        mw_mix_inv += w_sf / mw;
    }
    if (mw_mix_inv <= 0.0) {
        throw ValueError("Solvent molecular weight mixture is invalid.");
    }
    return 1.0 / mw_mix_inv;
}

void validate_activity_inputs_cpp(
    const vector<double>& x,
    const add_args& args,
    const vector<int>& cation_indices,
    const vector<int>& anion_indices,
    const vector<int>& pair_cation_indices,
    const vector<int>& pair_anion_indices,
    const vector<int>& pair_nu_cation,
    const vector<int>& pair_nu_anion
) {
    if (args.z.empty() || std::all_of(args.z.begin(), args.z.end(), [](double v) { return std::abs(v) <= 1e-12; })) {
        throw ValueError("activity_coefficient requires ionic species (non-zero z).");
    }
    if (args.mw.size() != x.size()) {
        throw ValueError("activity_coefficient requires params['MW'] to be present and aligned with x.");
    }
    if (cation_indices.empty() || anion_indices.empty()) {
        throw ValueError("activity_coefficient requires at least one cation and one anion.");
    }
    if (pair_cation_indices.size() != pair_anion_indices.size()
        || pair_cation_indices.size() != pair_nu_cation.size()
        || pair_cation_indices.size() != pair_nu_anion.size()) {
        throw ValueError("Invalid ionic pair metadata for activity_coefficient.");
    }
}

void assign_activity_metadata_cpp(
    ActivityCoefficientNative& out,
    const vector<double>& x,
    const vector<int>& cation_indices,
    const vector<int>& anion_indices,
    const vector<int>& solvent_indices,
    const vector<int>& pair_cation_indices,
    const vector<int>& pair_anion_indices,
    const vector<int>& pair_nu_cation,
    const vector<int>& pair_nu_anion,
    bool has_solvent_override,
    int solvent_override_index
) {
    out.cation_indices = cation_indices;
    out.anion_indices = anion_indices;
    out.solvent_indices = solvent_indices;
    out.pair_cation_indices = pair_cation_indices;
    out.pair_anion_indices = pair_anion_indices;
    out.pair_nu_cation = pair_nu_cation;
    out.pair_nu_anion = pair_nu_anion;
    out.solvent_index = resolve_solvent_index_cpp(solvent_indices, x, has_solvent_override, solvent_override_index);
}

void assign_activity_aux_cpp(
    ActivityCoefficientNative& out,
    double t,
    double rho,
    const vector<double>& x,
    const add_args& args,
    bool include_aux
) {
    out.component_activity_coefficients = miac_gamma_vector_cpp(t, rho, x, args);
    if (include_aux) {
        out.solvation_free_energy = gsolv_values_cpp(t, rho, x, args);
        if (out.solvation_free_energy.size() != x.size()) {
            throw ValueError("Unexpected solvation_free_energy payload size in activity_coefficient.");
        }
    } else {
        out.solvation_free_energy.assign(x.size(), std::numeric_limits<double>::quiet_NaN());
    }
}

void assign_pair_activity_cpp(
    ActivityCoefficientNative& out,
    const vector<double>& x,
    const add_args& args,
    const vector<int>& solvent_indices,
    bool has_solvent_override
) {
    vector<int> solvent_pool = has_solvent_override ? vector<int>{out.solvent_index} : solvent_indices;
    double mass_solvent = 0.0;
    for (int idx : solvent_pool) {
        mass_solvent += x[idx] * normalize_mw_cpp(args.mw[idx]);
    }
    if (mass_solvent <= 0.0) {
        throw ValueError("Solvent mass is zero; check solvent mole fraction and MW.");
    }
    double mw_mix = solvent_pool_mix_mw_cpp(x, args, solvent_pool);

    out.mean_ionic_activity_coefficients_mole_fraction.reserve(out.pair_cation_indices.size());
    out.mean_ionic_activity_coefficients_molality.reserve(out.pair_cation_indices.size());
    out.pair_molality.reserve(out.pair_cation_indices.size());
    out.pair_conversion_factor.reserve(out.pair_cation_indices.size());
    for (size_t k = 0; k < out.pair_cation_indices.size(); ++k) {
        int ic = out.pair_cation_indices[k];
        int ia = out.pair_anion_indices[k];
        double nu_cat = static_cast<double>(out.pair_nu_cation[k]);
        double nu_an = static_cast<double>(out.pair_nu_anion[k]);
        double sum_nu = nu_cat + nu_an;
        double ln_gamma_pm = (nu_cat * std::log(std::max(out.component_activity_coefficients[ic], 1e-300))
            + nu_an * std::log(std::max(out.component_activity_coefficients[ia], 1e-300))) / sum_nu;
        double gamma_pm_x = std::exp(ln_gamma_pm);
        double n_salt = 0.5 * (x[ic] / nu_cat + x[ia] / nu_an);
        double m_salt = n_salt / mass_solvent;
        double conversion = 1.0 + mw_mix * m_salt * sum_nu;
        out.mean_ionic_activity_coefficients_mole_fraction.push_back(gamma_pm_x);
        out.mean_ionic_activity_coefficients_molality.push_back(gamma_pm_x / conversion);
        out.pair_molality.push_back(m_salt);
        out.pair_conversion_factor.push_back(conversion);
    }
}

double osmotic_coefficient_cpp(
    double t,
    double rho,
    double p,
    int phase,
    const vector<double>& x,
    const add_args& args,
    int solvent_index
) {
    double mw_solvent = normalize_mw_cpp(args.mw[solvent_index]);
    if (mw_solvent <= 0.0) {
        throw ValueError("Solvent molecular weight must be positive.");
    }
    if (x[solvent_index] <= 0.0) {
        throw ValueError("Selected solvent mole fraction must be positive for osmotic coefficient.");
    }
    vector<double> x0(x.size(), 0.0);
    x0[solvent_index] = 1.0;
    int ref_phase = (phase == 0 || phase == 1) ? phase : ((rho < 900.0) ? 1 : 0);
    vector<double> fugcoef = fugcoef_cpp(t, rho, x, args);
    double rho0 = den_cpp(t, p, x0, ref_phase, args);
    vector<double> fugcoef0 = fugcoef_cpp(t, rho0, x0, args);
    double gamma_solvent = fugcoef[solvent_index] / fugcoef0[solvent_index];
    double molality_sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        if (static_cast<int>(i) == solvent_index) {
            continue;
        }
        molality_sum += x[i] / (x[solvent_index] * mw_solvent);
    }
    if (molality_sum <= 0.0) {
        throw ValueError("Total molality is zero; osmotic coefficient is undefined.");
    }
    return -std::log(x[solvent_index] * gamma_solvent) / (mw_solvent * molality_sum);
}

ActivityCoefficientNative activity_coefficient_values_impl_cpp(
    double t,
    double rho,
    double p,
    int phase,
    const vector<double>& x,
    const add_args& args,
    const vector<int>& cation_indices,
    const vector<int>& anion_indices,
    const vector<int>& solvent_indices,
    const vector<int>& pair_cation_indices,
    const vector<int>& pair_anion_indices,
    const vector<int>& pair_nu_cation,
    const vector<int>& pair_nu_anion,
    bool include_aux,
    bool has_solvent_override,
    int solvent_override_index
) {
    validate_activity_inputs_cpp(
        x,
        args,
        cation_indices,
        anion_indices,
        pair_cation_indices,
        pair_anion_indices,
        pair_nu_cation,
        pair_nu_anion
    );

    ActivityCoefficientNative out;
    assign_activity_metadata_cpp(
        out,
        x,
        cation_indices,
        anion_indices,
        solvent_indices,
        pair_cation_indices,
        pair_anion_indices,
        pair_nu_cation,
        pair_nu_anion,
        has_solvent_override,
        solvent_override_index
    );
    assign_activity_aux_cpp(out, t, rho, x, args, include_aux);
    assign_pair_activity_cpp(out, x, args, solvent_indices, has_solvent_override);

    if (include_aux) {
        out.osmotic_coefficient = osmotic_coefficient_cpp(t, rho, p, phase, x, args, out.solvent_index);
    } else {
        out.osmotic_coefficient = std::numeric_limits<double>::quiet_NaN();
    }
    return out;
}

} // namespace miac_detail

// EqID: gamma_asym_inf
// EqID: lngamma_asym_inf
ActivityCoefficientNative activity_coefficient_values_cpp(
    double t,
    double rho,
    double p,
    int phase,
    const vector<double>& x,
    const add_args& args,
    const vector<int>& cation_indices,
    const vector<int>& anion_indices,
    const vector<int>& solvent_indices,
    const vector<int>& pair_cation_indices,
    const vector<int>& pair_anion_indices,
    const vector<int>& pair_nu_cation,
    const vector<int>& pair_nu_anion,
    bool include_aux,
    bool has_solvent_override,
    int solvent_override_index
) {
    return miac_detail::activity_coefficient_values_impl_cpp(
        t,
        rho,
        p,
        phase,
        x,
        args,
        cation_indices,
        anion_indices,
        solvent_indices,
        pair_cation_indices,
        pair_anion_indices,
        pair_nu_cation,
        pair_nu_anion,
        include_aux,
        has_solvent_override,
        solvent_override_index
    );
}
