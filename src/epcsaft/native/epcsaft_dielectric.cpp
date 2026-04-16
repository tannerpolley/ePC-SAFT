#include "epcsaft_core_internal.h"
#include "epcsaft_autodiff_internal.h"

using namespace thermo_detail;

namespace {

double mixed_dielectric_constant_cpp(const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule=8 requires params['z'] as an array with length equal to ncomp.");
    }
    if (cppargs.mixed_rel_perm_a.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_b.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_c.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_mask.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule=8 requires mixed relative-permittivity arrays with length equal to ncomp.");
    }

    double x_sol = 0.0;
    double x_water = 0.0;
    double x_org = 0.0;
    double eps_org_num = 0.0;
    double a_num = 0.0;
    double b_num = 0.0;
    double c_num = 0.0;
    bool needs_coeffs = false;

    int water_idx = cppargs.mixed_rel_perm_water_index;
    bool has_water_component = (
        water_idx >= 0 &&
        water_idx < ncomp &&
        std::abs(cppargs.z[water_idx]) <= 1e-12
    );

    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) > 1e-12) {
            continue;
        }
        x_sol += x[i];
        if (has_water_component && i == water_idx) {
            x_water += x[i];
            continue;
        }
        x_org += x[i];
        eps_org_num += x[i] * cppargs.dielc[i];
        if (x[i] > 0.0) {
            needs_coeffs = true;
            if (cppargs.mixed_rel_perm_mask[i] == 0) {
                if (x_water > 0.0 || has_water_component) {
                    throw ValueError("dielc_rule=8 is missing mixed relative-permittivity coefficients for an organic solvent.");
                }
            } else {
                a_num += x[i] * cppargs.mixed_rel_perm_a[i];
                b_num += x[i] * cppargs.mixed_rel_perm_b[i];
                c_num += x[i] * cppargs.mixed_rel_perm_c[i];
            }
        }
    }

    if (x_sol <= 0.0) {
        throw ValueError("dielc_rule=8 requires at least one solvent species (z=0).");
    }
    if (x_org <= DBL_EPSILON) {
        if (!has_water_component) {
            throw ValueError("dielc_rule=8 requires at least one organic solvent or water component.");
        }
        return cppargs.dielc[water_idx];
    }
    if (!has_water_component || x_water <= DBL_EPSILON) {
        return eps_org_num / x_org;
    }
    if (x_water >= x_sol - DBL_EPSILON) {
        return cppargs.dielc[water_idx];
    }
    if (!needs_coeffs) {
        throw ValueError("dielc_rule=8 requires mixed relative-permittivity coefficients for the organic solvent phase.");
    }

    double xw_sf = x_water / x_sol;
    double eps_org = eps_org_num / x_org;
    double a_eff = a_num / x_org;
    double b_eff = b_num / x_org;
    double c_eff = c_num / x_org;
    return eps_org + ((a_eff * xw_sf + b_eff) * xw_sf + c_eff) * xw_sf;
}

AutoDual mixed_dielectric_constant_ad_cpp(const vector<AutoDual> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule=8 requires params['z'] as an array with length equal to ncomp.");
    }
    if (cppargs.mixed_rel_perm_a.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_b.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_c.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_mask.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule=8 requires mixed relative-permittivity arrays with length equal to ncomp.");
    }

    AutoDual x_sol = 0.0;
    AutoDual x_water = 0.0;
    AutoDual x_org = 0.0;
    AutoDual eps_org_num = 0.0;
    AutoDual a_num = 0.0;
    AutoDual b_num = 0.0;
    AutoDual c_num = 0.0;
    bool needs_coeffs = false;

    int water_idx = cppargs.mixed_rel_perm_water_index;
    bool has_water_component = (
        water_idx >= 0 &&
        water_idx < ncomp &&
        std::abs(cppargs.z[water_idx]) <= 1e-12
    );

    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) > 1e-12) {
            continue;
        }
        x_sol += x[i];
        if (has_water_component && i == water_idx) {
            x_water += x[i];
            continue;
        }
        x_org += x[i];
        eps_org_num += x[i] * cppargs.dielc[i];
        if (scalar_value(x[i]) > 0.0) {
            needs_coeffs = true;
            if (cppargs.mixed_rel_perm_mask[i] == 0) {
                if (scalar_value(x_water) > 0.0 || has_water_component) {
                    throw ValueError("dielc_rule=8 is missing mixed relative-permittivity coefficients for an organic solvent.");
                }
            } else {
                a_num += x[i] * cppargs.mixed_rel_perm_a[i];
                b_num += x[i] * cppargs.mixed_rel_perm_b[i];
                c_num += x[i] * cppargs.mixed_rel_perm_c[i];
            }
        }
    }

    if (scalar_value(x_sol) <= 0.0) {
        throw ValueError("dielc_rule=8 requires at least one solvent species (z=0).");
    }
    if (scalar_value(x_org) <= DBL_EPSILON) {
        if (!has_water_component) {
            throw ValueError("dielc_rule=8 requires at least one organic solvent or water component.");
        }
        return AutoDual(cppargs.dielc[water_idx], 0.0);
    }
    if (!has_water_component || scalar_value(x_water) <= DBL_EPSILON) {
        return eps_org_num / x_org;
    }
    if (scalar_value(x_water) >= scalar_value(x_sol) - DBL_EPSILON) {
        return AutoDual(cppargs.dielc[water_idx], 0.0);
    }
    if (!needs_coeffs) {
        throw ValueError("dielc_rule=8 requires mixed relative-permittivity coefficients for the organic solvent phase.");
    }

    AutoDual xw_sf = x_water / x_sol;
    AutoDual eps_org = eps_org_num / x_org;
    AutoDual a_eff = a_num / x_org;
    AutoDual b_eff = b_num / x_org;
    AutoDual c_eff = c_num / x_org;
    return eps_org + ((a_eff * xw_sf + b_eff) * xw_sf + c_eff) * xw_sf;
}

void dielectric_inputs_valid_cpp(const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.dielc.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("params['dielc'] must be an array with length equal to ncomp.");
    }
    if (cppargs.dielc_diff_mode != 0 && cppargs.dielc_diff_mode != 1 && cppargs.dielc_diff_mode != 2) {
        throw ValueError("Unknown dielc_diff_mode. Supported values are 0 (analytic), 1 (finite-diff), and 2 (autodiff).");
    }
    if (cppargs.hc_dadx_diff_mode != 0 && cppargs.hc_dadx_diff_mode != 1 && cppargs.hc_dadx_diff_mode != 2) {
        throw ValueError("Unknown hc_model dadx_differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).");
    }
    if (cppargs.disp_dadx_diff_mode != 0 && cppargs.disp_dadx_diff_mode != 1 && cppargs.disp_dadx_diff_mode != 2) {
        throw ValueError("Unknown disp_model dadx_differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).");
    }
    if (cppargs.assoc_dadx_diff_mode != 0 && cppargs.assoc_dadx_diff_mode != 1 && cppargs.assoc_dadx_diff_mode != 2) {
        throw ValueError("Unknown assoc_model dadx_differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).");
    }
    if (cppargs.born_diff_mode != 0 && cppargs.born_diff_mode != 1 && cppargs.born_diff_mode != 2 && cppargs.born_diff_mode != 3 && cppargs.born_diff_mode != 4) {
        throw ValueError("Unknown born_diff_mode. Supported values are 0 (analytic), 1 (finite-diff), 2 (Eq.133-style), 3 (no dielectric-concentration term), and 4 (autodiff).");
    }
    if (cppargs.d_ion_mode < 0 || cppargs.d_ion_mode > 2) {
        throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
    }
    if (cppargs.mu_DH_diff_mode != 0 && cppargs.mu_DH_diff_mode != 1 && cppargs.mu_DH_diff_mode != 2) {
        throw ValueError("Unknown mu_DH differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).");
    }
    if (cppargs.mu_DH_comp_dep_rel_perm != 0 && cppargs.mu_DH_comp_dep_rel_perm != 1) {
        throw ValueError("mu_DH comp_dep_rel_perm must be 0 or 1.");
    }
    if (cppargs.mu_DH_include_sum_term != 0 && cppargs.mu_DH_include_sum_term != 1) {
        throw ValueError("mu_DH include_sum_term must be 0 or 1.");
    }
    if (cppargs.include_born_model != 0 && cppargs.include_born_model != 1) {
        throw ValueError("include_born_model must be 0 or 1.");
    }
    if (cppargs.d_born_mode < 0 || cppargs.d_born_mode > 3) {
        throw ValueError("Unknown d_Born_mode. Supported values are 0, 1, 2, 3.");
    }
    if (cppargs.born_bulk_mode != 0 && cppargs.born_bulk_mode != 1) {
        throw ValueError("Unknown born bulk_mode. Supported values are mix/solvent (0/1).");
    }
    if (cppargs.mu_born_diff_mode != 0 && cppargs.mu_born_diff_mode != 1 && cppargs.mu_born_diff_mode != 2) {
        throw ValueError("Unknown mu_born differential_mode. Supported values are analytical/numerical/autodiff (0/1/2).");
    }
    if (cppargs.born_eps_mode != 0 && cppargs.born_eps_mode != 1) {
        throw ValueError("Unknown born_eps_mode. Supported values are 0 (eps_r,mix) and 1 (eps_r,solvent).");
    }
    if (cppargs.born_model < 0 || cppargs.born_model > 2) {
        throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
    }
    if (cppargs.born_model > 0 && cppargs.z.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("Born contribution requires params['z'] as an array with length equal to ncomp.");
    }
    if (cppargs.born_radius_model < 1 || cppargs.born_radius_model > 5) {
        throw ValueError("Unknown born_radius_model. Supported values are 1, 2, 3, 4, 5.");
    }
    if (cppargs.born_model == 1 && cppargs.born_radius_model == 5) {
        throw ValueError("born_model=1 supports born_radius_model values 1, 2, 3, 4.");
    }
    if (cppargs.born_model > 0 && (cppargs.born_radius_model == 4 || cppargs.born_radius_model == 5)) {
        if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
            throw ValueError("born_radius_model 4/5 requires ionic charge array params['z'] with length ncomp.");
        }
        for (int i = 0; i < ncomp; i++) {
            if (is_ion_species(cppargs, i) &&
                (cppargs.d_born.size() <= static_cast<size_t>(i) || cppargs.d_born[i] <= 0.0)) {
                throw ValueError("born_radius_model 4/5 requires positive ionic params['d_born'] values.");
            }
        }
    }
    int rule = cppargs.dielc_rule;
    if (rule < 0 || rule > 8) {
        throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8.");
    }
    if ((rule == 2 || rule == 3 || rule == 4 || rule == 5) &&
        cppargs.mw.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule requires params['MW'] as an array with length equal to ncomp.");
    }
    if ((rule == 3 || rule == 4 || rule == 5 || rule == 6) &&
        cppargs.z.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule requires params['z'] as an array with length equal to ncomp.");
    }
    if (rule == 8) {
        if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
            throw ValueError("dielc_rule=8 requires params['z'] as an array with length equal to ncomp.");
        }
        if (cppargs.mixed_rel_perm_a.size() != static_cast<size_t>(ncomp) ||
            cppargs.mixed_rel_perm_b.size() != static_cast<size_t>(ncomp) ||
            cppargs.mixed_rel_perm_c.size() != static_cast<size_t>(ncomp) ||
            cppargs.mixed_rel_perm_mask.size() != static_cast<size_t>(ncomp)) {
            throw ValueError("dielc_rule=8 requires mixed relative-permittivity arrays with length equal to ncomp.");
        }
    }
}

}  // namespace

AutoDual dielectric_constant_rule_autodiff_cpp(int rule, const vector<AutoDual> &x, const add_args &cppargs) {
    const double alpha = 7.01;
    int ncomp = static_cast<int>(x.size());
    if (rule == 0) {
        return AutoDual(*std::max_element(cppargs.dielc.begin(), cppargs.dielc.end()), 0.0);
    }
    if (rule == 1) {
        AutoDual eps = 0.0;
        for (int i = 0; i < ncomp; i++) {
            eps += x[i] * cppargs.dielc[i];
        }
        return eps;
    }
    if (rule == 7) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=7 requires at least one solvent species (z=0).");
        }
        if (idx_ion.empty()) {
            throw ValueError("dielc_rule=7 requires at least one ionic species (z!=0).");
        }
        AutoDual x_sol = 0.0;
        AutoDual eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            x_sol += x[idx];
            eps_sol_num += x[idx] * cppargs.dielc[idx];
        }
        AutoDual eps_sol = 0.0;
        if (scalar_value(x_sol) > 1.0e-16) {
            eps_sol = eps_sol_num / x_sol;
        }
        else {
            double eps_sol_const = 0.0;
            for (int idx : idx_sol) eps_sol_const += cppargs.dielc[idx];
            eps_sol = AutoDual(eps_sol_const / static_cast<double>(idx_sol.size()), 0.0);
        }
        double eps_salt = 0.0;
        for (int idx : idx_ion) eps_salt += cppargs.dielc[idx];
        eps_salt /= static_cast<double>(idx_ion.size());
        return eps_sol * x_sol + eps_salt * (1.0 - x_sol);
    }
    if (rule == 8) {
        return mixed_dielectric_constant_ad_cpp(x, cppargs);
    }
    if (rule == 2) {
        AutoDual mw_bar = 0.0;
        AutoDual num = 0.0;
        for (int i = 0; i < ncomp; i++) {
            mw_bar += x[i] * cppargs.mw[i];
            num += x[i] * cppargs.mw[i] * cppargs.dielc[i];
        }
        if (scalar_value(mw_bar) <= 0.0) {
            throw ValueError("Average molecular weight must be positive for dielc_rule=2.");
        }
        return num / mw_bar;
    }
    if (rule == 3) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=3 requires at least one solvent species (z=0).");
        }
        AutoDual mw_sol = 0.0;
        AutoDual eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx] * cppargs.mw[idx];
            eps_sol_num += x[idx] * cppargs.mw[idx] * cppargs.dielc[idx];
        }
        if (scalar_value(mw_sol) <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.");
        }
        AutoDual eps_sol_w = eps_sol_num / mw_sol;
        AutoDual x_sol = 0.0;
        AutoDual eps_ion = 0.0;
        for (int idx : idx_sol) x_sol += x[idx];
        for (int idx : idx_ion) eps_ion += x[idx] * cppargs.dielc[idx];
        return x_sol * eps_sol_w + eps_ion;
    }
    if (rule == 4 || rule == 5) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule requires at least one solvent species (z=0).");
        }
        AutoDual mw_sol = 0.0;
        AutoDual eps_sf_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx] * cppargs.mw[idx];
            eps_sf_num += x[idx] * cppargs.mw[idx] * cppargs.dielc[idx];
        }
        if (scalar_value(mw_sol) <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule.");
        }
        AutoDual eps_sf = eps_sf_num / mw_sol;
        AutoDual x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        return eps_sf / (1.0 + alpha * x_ion);
    }
    if (rule == 6) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=6 requires at least one solvent species (z=0).");
        }
        double eps_sf_const = 0.0;
        for (int idx : idx_sol) eps_sf_const += cppargs.dielc[idx];
        eps_sf_const /= static_cast<double>(idx_sol.size());
        AutoDual x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        return AutoDual(eps_sf_const, 0.0) / (1.0 + alpha * x_ion);
    }
    throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8.");
}

double dielectric_constant_rule_cpp(int rule, const vector<double> &x, const add_args &cppargs) {
    const double alpha = 7.01;
    int ncomp = static_cast<int>(x.size());
    if (rule == 0) {
        return *std::max_element(cppargs.dielc.begin(), cppargs.dielc.end());
    }
    if (rule == 1) {
        double eps = 0.0;
        for (int i = 0; i < ncomp; i++) {
            eps += x[i]*cppargs.dielc[i];
        }
        return eps;
    }
    if (rule == 7) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=7 requires at least one solvent species (z=0).");
        }
        if (idx_ion.empty()) {
            throw ValueError("dielc_rule=7 requires at least one ionic species (z!=0).");
        }
        double x_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            x_sol += x[idx];
            eps_sol_num += x[idx] * cppargs.dielc[idx];
        }
        if (x_sol < 0.0 || x_sol > 1.0) {
            throw ValueError("dielc_rule=7 encountered invalid solvent mole fraction.");
        }
        double eps_sol = 0.0;
        if (x_sol > 1.0e-16) {
            eps_sol = eps_sol_num / x_sol;
        }
        else {
            for (int idx : idx_sol) eps_sol += cppargs.dielc[idx];
            eps_sol /= static_cast<double>(idx_sol.size());
        }
        double eps_salt = 0.0;
        for (int idx : idx_ion) {
            eps_salt += cppargs.dielc[idx];
        }
        eps_salt /= static_cast<double>(idx_ion.size());
        return eps_sol * x_sol + eps_salt * (1.0 - x_sol);
    }
    if (rule == 8) {
        return mixed_dielectric_constant_cpp(x, cppargs);
    }
    if (rule == 2) {
        double mw_bar = 0.0;
        double num = 0.0;
        for (int i = 0; i < ncomp; i++) {
            mw_bar += x[i]*cppargs.mw[i];
            num += x[i]*cppargs.mw[i]*cppargs.dielc[i];
        }
        if (mw_bar <= 0.0) {
            throw ValueError("Average molecular weight must be positive for dielc_rule=2.");
        }
        return num/mw_bar;
    }
    if (rule == 3) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=3 requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx]*cppargs.mw[idx];
            eps_sol_num += x[idx]*cppargs.mw[idx]*cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.");
        }
        double eps_sol_w = eps_sol_num/mw_sol;
        double x_sol = 0.0;
        double eps_ion = 0.0;
        for (int idx : idx_sol) x_sol += x[idx];
        for (int idx : idx_ion) eps_ion += x[idx]*cppargs.dielc[idx];
        return x_sol*eps_sol_w + eps_ion;
    }
    if (rule == 4 || rule == 5) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sf_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx]*cppargs.mw[idx];
            eps_sf_num += x[idx]*cppargs.mw[idx]*cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule.");
        }
        double eps_sf = eps_sf_num/mw_sol;
        double x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        return eps_sf/(1.0 + alpha*x_ion);
    }
    if (rule == 6) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=6 requires at least one solvent species (z=0).");
        }
        double eps_sf_const = 0.0;
        for (int idx : idx_sol) eps_sf_const += cppargs.dielc[idx];
        eps_sf_const /= static_cast<double>(idx_sol.size());
        double x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        return eps_sf_const/(1.0 + alpha*x_ion);
    }
    throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8.");
}

vector<double> dielectric_derivative_rule_cpp(int rule, const vector<double> &x, const add_args &cppargs) {
    const double alpha = 7.01;
    int ncomp = static_cast<int>(x.size());
    vector<double> deps_dx(ncomp, 0.0);
    if (rule == 0) {
        return deps_dx;
    }
    if (rule == 1) {
        return cppargs.dielc;
    }
    if (rule == 7) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=7 requires at least one solvent species (z=0).");
        }
        if (idx_ion.empty()) {
            throw ValueError("dielc_rule=7 requires at least one ionic species (z!=0).");
        }
        double x_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            x_sol += x[idx];
            eps_sol_num += x[idx] * cppargs.dielc[idx];
        }
        double eps_sol = 0.0;
        if (x_sol > 1.0e-16) {
            eps_sol = eps_sol_num / x_sol;
        }
        else {
            for (int idx : idx_sol) eps_sol += cppargs.dielc[idx];
            eps_sol /= static_cast<double>(idx_sol.size());
        }
        double eps_salt = 0.0;
        for (int idx : idx_ion) eps_salt += cppargs.dielc[idx];
        eps_salt /= static_cast<double>(idx_ion.size());
        for (int idx : idx_sol) deps_dx[idx] = eps_sol;
        for (int idx : idx_ion) deps_dx[idx] = eps_salt;
        return deps_dx;
    }
    if (rule == 8) {
        return dielectric_derivative_rule_fd_cpp(rule, x, cppargs);
    }
    if (rule == 2) {
        double mw_bar = 0.0;
        double eps_mix_num = 0.0;
        for (int i = 0; i < ncomp; i++) {
            mw_bar += x[i]*cppargs.mw[i];
            eps_mix_num += x[i]*cppargs.mw[i]*cppargs.dielc[i];
        }
        if (mw_bar <= 0.0) {
            throw ValueError("Average molecular weight must be positive for dielc_rule=2.");
        }
        double eps_mix = eps_mix_num/mw_bar;
        for (int i = 0; i < ncomp; i++) {
            deps_dx[i] = (cppargs.mw[i]/mw_bar)*(cppargs.dielc[i] - eps_mix);
        }
        return deps_dx;
    }
    if (rule == 3) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=3 requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx]*cppargs.mw[idx];
            eps_sol_num += x[idx]*cppargs.mw[idx]*cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.");
        }
        double eps_sol_w = eps_sol_num/mw_sol;
        double x_sol = 0.0;
        for (int idx : idx_sol) x_sol += x[idx];
        for (int idx : idx_sol) {
            deps_dx[idx] = eps_sol_w + x_sol*(cppargs.mw[idx]/mw_sol)*(cppargs.dielc[idx] - eps_sol_w);
        }
        for (int idx : idx_ion) {
            deps_dx[idx] = cppargs.dielc[idx];
        }
        return deps_dx;
    }
    if (rule == 4) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=4 requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sf_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx]*cppargs.mw[idx];
            eps_sf_num += x[idx]*cppargs.mw[idx]*cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=4.");
        }
        double eps_sf = eps_sf_num/mw_sol;
        double x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        double den = 1.0 + alpha*x_ion;
        for (int idx : idx_sol) {
            deps_dx[idx] = (1.0/den)*(cppargs.mw[idx]/mw_sol)*(cppargs.dielc[idx] - eps_sf);
        }
        for (int idx : idx_ion) {
            deps_dx[idx] = -alpha*eps_sf/(den*den);
        }
        return deps_dx;
    }
    if (rule == 5) {
        return cppargs.dielc;
    }
    if (rule == 6) {
        for (int i = 0; i < ncomp; i++) {
            deps_dx[i] = (std::abs(cppargs.z[i]) <= 1e-12) ? cppargs.dielc[i] : 0.0;
        }
        return deps_dx;
    }
    throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8.");
}

vector<double> dielectric_derivative_rule_fd_cpp(int rule, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    vector<double> deps_dx(ncomp, 0.0);
    double f0 = dielectric_constant_rule_cpp(rule, x, cppargs);
    for (int i = 0; i < ncomp; i++) {
        double h = 1e-6*std::max(1.0, std::abs(x[i]));
        vector<double> xp = x;
        xp[i] += h;
        double fp = dielectric_constant_rule_cpp(rule, xp, cppargs);
        if (x[i] - h >= 0.0) {
            vector<double> xm = x;
            xm[i] -= h;
            double fm = dielectric_constant_rule_cpp(rule, xm, cppargs);
            deps_dx[i] = (fp - fm)/(2.0*h);
        }
        else {
            deps_dx[i] = (fp - f0)/h;
        }
        if (!std::isfinite(deps_dx[i])) {
            throw ValueError("Non-finite dielectric finite-difference derivative.");
        }
    }
    return deps_dx;
}

vector<double> dielectric_derivative_rule_ad_cpp(int rule, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    vector<double> deps_dx(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        vector<AutoDual> x_dual(ncomp, AutoDual(0.0, 0.0));
        for (int j = 0; j < ncomp; ++j) {
            x_dual[j] = AutoDual(x[j], (i == j) ? 1.0 : 0.0);
        }
        deps_dx[i] = scalar_derivative(dielectric_constant_rule_autodiff_cpp(rule, x_dual, cppargs));
        if (!std::isfinite(deps_dx[i])) {
            throw ValueError("Non-finite dielectric autodiff derivative.");
        }
    }
    return deps_dx;
}

DielectricState dielectric_state_cpp(const vector<double> &x, const add_args &cppargs) {
    dielectric_inputs_valid_cpp(x, cppargs);
    DielectricState state;
    state.eps = dielectric_constant_rule_cpp(cppargs.dielc_rule, x, cppargs);
    if (cppargs.dielc_diff_mode == 0 && cppargs.dielc_rule != 8) {
        state.deps_dx = dielectric_derivative_rule_cpp(cppargs.dielc_rule, x, cppargs);
    }
    else if (cppargs.dielc_diff_mode == 2) {
        state.deps_dx = dielectric_derivative_rule_ad_cpp(cppargs.dielc_rule, x, cppargs);
    }
    else {
        state.deps_dx = dielectric_derivative_rule_fd_cpp(cppargs.dielc_rule, x, cppargs);
    }
    return state;
}
