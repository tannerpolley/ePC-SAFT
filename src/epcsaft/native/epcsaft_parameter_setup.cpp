#include "epcsaft_core_internal.h"
#include "epcsaft_autodiff_internal.h"
#include "contributions/epcsaft_contrib_internal.h"

using thermo_detail::AresContributionKind;
using thermo_detail::DielectricState;
using thermo_detail::DispersionPolynomialState;
using thermo_detail::MixtureState;

namespace thermo_detail {
namespace parameter_setup_detail {

template <typename Scalar>
// EqID: epsr_mix_rule
// EqID: epsr_mix_suppressed
Scalar mixed_dielectric_constant_scalar_cpp(const vector<Scalar> &x, const add_args &cppargs) {
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

    Scalar x_sol = scalar_constant<Scalar>(0.0);
    Scalar x_water = scalar_constant<Scalar>(0.0);
    Scalar x_org = scalar_constant<Scalar>(0.0);
    Scalar eps_org_num = scalar_constant<Scalar>(0.0);
    Scalar a_num = scalar_constant<Scalar>(0.0);
    Scalar b_num = scalar_constant<Scalar>(0.0);
    Scalar c_num = scalar_constant<Scalar>(0.0);
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
        return scalar_constant<Scalar>(cppargs.dielc[water_idx]);
    }
    if (!has_water_component || scalar_value(x_water) <= DBL_EPSILON) {
        return eps_org_num / x_org;
    }
    if (scalar_value(x_water) >= scalar_value(x_sol) - DBL_EPSILON) {
        return scalar_constant<Scalar>(cppargs.dielc[water_idx]);
    }
    if (!needs_coeffs) {
        throw ValueError("dielc_rule=8 requires mixed relative-permittivity coefficients for the organic solvent phase.");
    }

    Scalar xw_sf = x_water / x_sol;
    Scalar eps_org = eps_org_num / x_org;
    Scalar a_eff = a_num / x_org;
    Scalar b_eff = b_num / x_org;
    Scalar c_eff = c_num / x_org;
    return eps_org + ((a_eff * xw_sf + b_eff) * xw_sf + c_eff) * xw_sf;
}

template <typename Scalar>
// EqID: mw_bar
// EqID: mw_solvent_bar
// EqID: solvent_ion_sets
// EqID: x_solvent_total
// EqID: epsr_solvent_mass
// EqID: epsr_salt_free
// EqID: epsr_sf
// EqID: x_ion_total
Scalar dielectric_constant_rule_scalar_cpp(int rule, const vector<Scalar> &x, const add_args &cppargs) {
    const double alpha = 7.01;
    const Scalar one = scalar_constant<Scalar>(1.0);
    int ncomp = static_cast<int>(x.size());
    if (rule == 0) {
        return scalar_constant<Scalar>(*std::max_element(cppargs.dielc.begin(), cppargs.dielc.end()));
    }
    if (rule == 1) {
        Scalar eps = scalar_constant<Scalar>(0.0);
        for (int i = 0; i < ncomp; i++) {
            eps += x[i] * cppargs.dielc[i];
        }
        return eps;
    }

    vector<int> idx_sol;
    vector<int> idx_ion;
    if (rule == 3 || rule == 4 || rule == 5 || rule == 6 || rule == 7 || rule == 9) {
        idx_sol.reserve(ncomp);
        idx_ion.reserve(ncomp);
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
    }

    if (rule == 7) {
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=7 requires at least one solvent species (z=0).");
        }
        if (idx_ion.empty()) {
            throw ValueError("dielc_rule=7 requires at least one ionic species (z!=0).");
        }
        Scalar x_sol = scalar_constant<Scalar>(0.0);
        Scalar eps_sol_num = scalar_constant<Scalar>(0.0);
        for (int idx : idx_sol) {
            x_sol += x[idx];
            eps_sol_num += x[idx] * cppargs.dielc[idx];
        }
        if (scalar_value(x_sol) < 0.0 || scalar_value(x_sol) > 1.0) {
            throw ValueError("dielc_rule=7 encountered invalid solvent mole fraction.");
        }
        Scalar eps_sol = scalar_constant<Scalar>(0.0);
        if (scalar_value(x_sol) > 1.0e-16) {
            eps_sol = eps_sol_num / x_sol;
        } else {
            double eps_sol_const = 0.0;
            for (int idx : idx_sol) eps_sol_const += cppargs.dielc[idx];
            eps_sol = scalar_constant<Scalar>(eps_sol_const / static_cast<double>(idx_sol.size()));
        }
        double eps_salt = 0.0;
        for (int idx : idx_ion) {
            eps_salt += cppargs.dielc[idx];
        }
        eps_salt /= static_cast<double>(idx_ion.size());
        return eps_sol * x_sol + eps_salt * (one - x_sol);
    }
    if (rule == 8) {
        return mixed_dielectric_constant_scalar_cpp(x, cppargs);
    }
    if (rule == 2) {
        Scalar mw_bar = scalar_constant<Scalar>(0.0);
        Scalar num = scalar_constant<Scalar>(0.0);
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
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=3 requires at least one solvent species (z=0).");
        }
        Scalar mw_sol = scalar_constant<Scalar>(0.0);
        Scalar eps_sol_num = scalar_constant<Scalar>(0.0);
        for (int idx : idx_sol) {
            mw_sol += x[idx] * cppargs.mw[idx];
            eps_sol_num += x[idx] * cppargs.mw[idx] * cppargs.dielc[idx];
        }
        if (scalar_value(mw_sol) <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.");
        }
        Scalar eps_sol_w = eps_sol_num / mw_sol;
        Scalar x_sol = scalar_constant<Scalar>(0.0);
        Scalar eps_ion = scalar_constant<Scalar>(0.0);
        for (int idx : idx_sol) x_sol += x[idx];
        for (int idx : idx_ion) eps_ion += x[idx] * cppargs.dielc[idx];
        return x_sol * eps_sol_w + eps_ion;
    }
    if (rule == 9) {
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=9 requires at least one solvent species (z=0).");
        }
        Scalar mw_sol = scalar_constant<Scalar>(0.0);
        Scalar eps_sol_num = scalar_constant<Scalar>(0.0);
        for (int idx : idx_sol) {
            mw_sol += x[idx] * cppargs.mw[idx];
            eps_sol_num += x[idx] * cppargs.mw[idx] * cppargs.dielc[idx];
        }
        if (scalar_value(mw_sol) <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=9.");
        }
        return eps_sol_num / mw_sol;
    }
    if (rule == 4 || rule == 5) {
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule requires at least one solvent species (z=0).");
        }
        Scalar mw_sol = scalar_constant<Scalar>(0.0);
        Scalar eps_sf_num = scalar_constant<Scalar>(0.0);
        for (int idx : idx_sol) {
            mw_sol += x[idx] * cppargs.mw[idx];
            eps_sf_num += x[idx] * cppargs.mw[idx] * cppargs.dielc[idx];
        }
        if (scalar_value(mw_sol) <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule.");
        }
        Scalar eps_sf = eps_sf_num / mw_sol;
        Scalar x_ion = scalar_constant<Scalar>(0.0);
        for (int idx : idx_ion) x_ion += x[idx];
        return eps_sf / (one + alpha * x_ion);
    }
    if (rule == 6) {
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=6 requires at least one solvent species (z=0).");
        }
        double eps_sf_const = 0.0;
        for (int idx : idx_sol) eps_sf_const += cppargs.dielc[idx];
        eps_sf_const /= static_cast<double>(idx_sol.size());
        Scalar x_ion = scalar_constant<Scalar>(0.0);
        for (int idx : idx_ion) x_ion += x[idx];
        return scalar_constant<Scalar>(eps_sf_const) / (one + alpha * x_ion);
    }
    throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8, 9.");
}

void dielectric_inputs_valid_cpp(const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.dielc.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("params['dielc'] must be an array with length equal to ncomp.");
    }
    if (cppargs.dielc_diff_mode < 0 || cppargs.dielc_diff_mode > 3) {
        throw ValueError("Unknown dielc_diff_mode. Supported values are 0 (analytic), 1 (finite-diff), 2 (autodiff), and 3 (auto).");
    }
    if (cppargs.hc_dadx_diff_mode < 0 || cppargs.hc_dadx_diff_mode > 3) {
        throw ValueError("Unknown hc_model dadx_differential_mode. Supported values are analytical/numerical/autodiff/auto (0/1/2/3).");
    }
    if (cppargs.disp_dadx_diff_mode < 0 || cppargs.disp_dadx_diff_mode > 3) {
        throw ValueError("Unknown disp_model dadx_differential_mode. Supported values are analytical/numerical/autodiff/auto (0/1/2/3).");
    }
    if (cppargs.assoc_dadx_diff_mode < 0 || cppargs.assoc_dadx_diff_mode > 3) {
        throw ValueError("Unknown assoc_model dadx_differential_mode. Supported values are analytical/numerical/autodiff/auto (0/1/2/3).");
    }
    if (cppargs.born_diff_mode < 0 || cppargs.born_diff_mode > 5) {
        throw ValueError("Unknown born_diff_mode. Supported values are 0 (analytic), 1 (finite-diff), 2 (Eq.133-style), 3 (no dielectric-concentration term), 4 (autodiff), and 5 (auto).");
    }
    if (cppargs.d_ion_mode < 0 || cppargs.d_ion_mode > 2) {
        throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
    }
    if (cppargs.mu_DH_diff_mode < 0 || cppargs.mu_DH_diff_mode > 3) {
        throw ValueError("Unknown mu_DH differential_mode. Supported values are analytical/numerical/autodiff/auto (0/1/2/3).");
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
    if (cppargs.mu_born_diff_mode < 0 || cppargs.mu_born_diff_mode > 3) {
        throw ValueError("Unknown mu_born differential_mode. Supported values are analytical/numerical/autodiff/auto (0/1/2/3).");
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
    if (rule < 0 || rule > 9) {
        throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8, 9.");
    }
    if ((rule == 2 || rule == 3 || rule == 4 || rule == 5 || rule == 9) &&
        cppargs.mw.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule requires params['MW'] as an array with length equal to ncomp.");
    }
    if ((rule == 3 || rule == 4 || rule == 5 || rule == 6 || rule == 9) &&
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

struct AutodiffMixtureState {
    vector<double> d;
    vector<double> s_ij;
    vector<double> e_ij;
    double den = 0.0;
    AutoDual m_avg = make_autodiff_scalar(0.0, 0.0);
    AutoDual m2es3 = make_autodiff_scalar(0.0, 0.0);
    AutoDual m2e2s3 = make_autodiff_scalar(0.0, 0.0);
};

struct AutodiffDispersionState {
    std::array<double, 7> a{};
    std::array<double, 7> b{};
    AutoDual I1 = make_autodiff_scalar(0.0, 0.0);
    AutoDual I2 = make_autodiff_scalar(0.0, 0.0);
    AutoDual C1 = make_autodiff_scalar(0.0, 0.0);
};

AutodiffMixtureState mixture_state_autodiff_cpp(double t, double rho, const vector<AutoDual> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    AutodiffMixtureState state;
    state.d.assign(ncomp, 0.0);
    state.e_ij.assign(ncomp * ncomp, 0.0);
    state.s_ij.assign(ncomp * ncomp, 0.0);
    state.den = rho * N_AV / 1.0e30;

    for (int i = 0; i < ncomp; ++i) {
        state.d[i] = cppargs.s[i] * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
        if (!cppargs.z.empty() && std::abs(cppargs.z[i]) > 1e-12) {
            state.d[i] = ion_diameter_cpp(i, t, cppargs);
        }
    }

    for (int i = 0; i < ncomp; ++i) {
        state.m_avg += x[i] * cppargs.m[i];
    }

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            state.s_ij[idx] = pair_sigma_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.e_ij[idx] = pair_epsilon_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.m2es3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * state.e_ij[idx] / t * scalar_pow(state.s_ij[idx], 3);
            state.m2e2s3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * scalar_pow(state.e_ij[idx] / t, 2.0) * scalar_pow(state.s_ij[idx], 3);
        }
    }
    return state;
}

AutodiffDispersionState dispersion_state_autodiff_cpp(const AutoDual &m_avg, const AutoDual &eta) {
    AutodiffDispersionState state;
    for (int i = 0; i < 7; ++i) {
        state.a[i] = kDispersionA0[i] + (1.0 - 1.0 / scalar_value(m_avg)) * kDispersionA1[i]
            + (1.0 - 1.0 / scalar_value(m_avg)) * (1.0 - 2.0 / scalar_value(m_avg)) * kDispersionA2[i];
        state.b[i] = kDispersionB0[i] + (1.0 - 1.0 / scalar_value(m_avg)) * kDispersionB1[i]
            + (1.0 - 1.0 / scalar_value(m_avg)) * (1.0 - 2.0 / scalar_value(m_avg)) * kDispersionB2[i];
        state.I1 += state.a[i] * scalar_pow(eta, i);
        state.I2 += state.b[i] * scalar_pow(eta, i);
    }
    state.C1 = 1.0 / (
        1.0
        + m_avg * (8.0 * eta - 2.0 * eta * eta) / scalar_pow(1.0 - eta, 4)
        + (1.0 - m_avg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * scalar_pow(eta, 3) - 2.0 * scalar_pow(eta, 4))
            / scalar_pow((1.0 - eta) * (2.0 - eta), 2)
    );
    return state;
}

double pair_diameter_cpp(double d_i, double d_j) {
    return d_i * d_j / (d_i + d_j);
}

// EqID: epsilon_assoc_mixing
// EqID: kappa_assoc_mixing
double association_volume_cpp(int comp_i, int comp_j, int ncomp, const vector<double> &s_ij, const add_args &cppargs) {
    int idxi = comp_i * ncomp + comp_i;
    int idxj = comp_j * ncomp + comp_j;
    double volume = std::sqrt(cppargs.vol_a[comp_i] * cppargs.vol_a[comp_j]) * std::pow(
        std::sqrt(s_ij[idxi] * s_ij[idxj]) / (0.5 * (s_ij[idxi] + s_ij[idxj])),
        3.0
    );
    if (!cppargs.k_hb.empty()) {
        volume *= (1.0 - cppargs.k_hb[comp_i * ncomp + comp_j]);
    }
    return volume;
}

double ion_diameter_cpp(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return cppargs.s[i];
    }
    int mode = cppargs.d_ion_mode;
    double sigma_i = cppargs.s[i];
    if (sigma_i <= 0.0) {
        throw ValueError("DH/ion diameter requires positive ionic sigma_i.");
    }
    if (mode == 0) {
        return sigma_i;
    }
    if (mode == 1) {
        return sigma_i * (1.0 - 0.12);
    }
    if (mode == 2) {
        return sigma_i * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
    }
    throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
}

double ion_diameter_cpp_dt(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return 0.0;
    }
    if (cppargs.d_ion_mode == 2) {
        double sigma_i = cppargs.s[i];
        double expo = std::exp(-3.0 * cppargs.e[i] / t);
        return -0.36 * sigma_i * cppargs.e[i] * expo / (t * t);
    }
    return 0.0;
}

// EqID: d_born_rule
double ion_born_radius_cpp(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return cppargs.s[i];
    }
    int mode = cppargs.d_born_mode;
    double sigma_i = cppargs.s[i];
    if (sigma_i <= 0.0) {
        throw ValueError("Born term requires positive ionic sigma_i.");
    }
    if (mode == 0) {
        return sigma_i;
    }
    if (mode == 1) {
        return sigma_i * (1.0 - 0.12);
    }
    if (mode == 2) {
        return sigma_i * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
    }
    if (mode == 3) {
        if (cppargs.d_born.size() <= static_cast<size_t>(i) || cppargs.d_born[i] <= 0.0) {
            throw ValueError("d_Born_mode=fitted_param requires positive ionic params['d_born'] values.");
        }
        return cppargs.d_born[i];
    }
    throw ValueError("Unknown d_Born_mode. Supported values are 0, 1, 2, 3.");
}

double ion_born_radius_cpp_dt(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return 0.0;
    }
    if (cppargs.d_born_mode == 2) {
        double sigma_i = cppargs.s[i];
        double expo = std::exp(-3.0 * cppargs.e[i] / t);
        return -0.36 * sigma_i * cppargs.e[i] * expo / (t * t);
    }
    return 0.0;
}

}  // namespace parameter_setup_detail
}  // namespace thermo_detail

using thermo_detail::parameter_setup_detail::AutodiffDispersionState;
using thermo_detail::parameter_setup_detail::AutodiffMixtureState;
using thermo_detail::parameter_setup_detail::dielectric_constant_rule_scalar_cpp;
using thermo_detail::parameter_setup_detail::dielectric_inputs_valid_cpp;
using thermo_detail::parameter_setup_detail::dispersion_state_autodiff_cpp;
using thermo_detail::parameter_setup_detail::ion_born_radius_cpp;
using thermo_detail::parameter_setup_detail::ion_diameter_cpp;
using thermo_detail::parameter_setup_detail::ion_diameter_cpp_dt;
using thermo_detail::parameter_setup_detail::mixed_dielectric_constant_scalar_cpp;
using thermo_detail::parameter_setup_detail::mixture_state_autodiff_cpp;
using thermo_detail::parameter_setup_detail::pair_epsilon_cpp;
using thermo_detail::parameter_setup_detail::pair_sigma_cpp;

ScalarContributionTerms make_scalar_terms(double hc, double disp, double assoc, double ion, double born, double total) {
    ScalarContributionTerms out;
    out.hc = hc;
    out.disp = disp;
    out.assoc = assoc;
    out.ion = ion;
    out.born = born;
    out.total = total;
    return out;
}

VectorContributionTerms make_vector_terms(
    const vector<double> &hc,
    const vector<double> &disp,
    const vector<double> &assoc,
    const vector<double> &ion,
    const vector<double> &born,
    const vector<double> &total
) {
    VectorContributionTerms out;
    out.hc = hc;
    out.disp = disp;
    out.assoc = assoc;
    out.ion = ion;
    out.born = born;
    out.total = total;
    return out;
}

// EqID: m_bar
// EqID: d_segment
// EqID: d_segment_dT
// EqID: d_ion_rule
// EqID: sigma_ij
// EqID: epsilon_ij_mixing
// EqID: epsilon_ij_ionic_zero
// EqID: m2epssigma3_bar
// EqID: m2eps2sigma3_bar
MixtureState mixture_state_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs, bool include_dt) {
    MixtureState state;
    int ncomp = static_cast<int>(x.size());
    state.d.assign(ncomp, 0.0);
    if (include_dt) {
        state.dd_dt.assign(ncomp, 0.0);
    }
    state.e_ij.assign(ncomp * ncomp, 0.0);
    state.s_ij.assign(ncomp * ncomp, 0.0);
    state.den = rho * N_AV / 1.0e30;

    for (int i = 0; i < ncomp; ++i) {
        state.d[i] = cppargs.s[i] * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
        if (include_dt) {
            state.dd_dt[i] = -0.36 * cppargs.s[i] * cppargs.e[i] * std::exp(-3.0 * cppargs.e[i] / t) / (t * t);
        }
        if (!cppargs.z.empty() && std::abs(cppargs.z[i]) > 1e-12) {
            state.d[i] = ion_diameter_cpp(i, t, cppargs);
            if (include_dt) {
                state.dd_dt[i] = ion_diameter_cpp_dt(i, t, cppargs);
            }
        }
    }

    for (int i = 0; i < ncomp; ++i) {
        state.m_avg += x[i] * cppargs.m[i];
    }

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            idx += 1;
            state.s_ij[idx] = pair_sigma_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.e_ij[idx] = pair_epsilon_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.m2es3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * state.e_ij[idx] / t * std::pow(state.s_ij[idx], 3);
            state.m2e2s3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * std::pow(state.e_ij[idx] / t, 2) * std::pow(state.s_ij[idx], 3);
        }
    }

    return state;
}

AutoDual dielectric_constant_rule_autodiff_cpp(int rule, const vector<AutoDual> &x, const add_args &cppargs) {
    return dielectric_constant_rule_scalar_cpp(rule, x, cppargs);
}

double dielectric_constant_rule_cpp(int rule, const vector<double> &x, const add_args &cppargs) {
    return dielectric_constant_rule_scalar_cpp(rule, x, cppargs);
}

// EqID: depsr_dxi_mole
// EqID: depsr_dxi_mass
// EqID: depsr_dxi_combo
// EqID: depsr_sf_dxi
// EqID: depsr_mix_dxi
// EqID: depsr_mix_dxi_piecewise
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
        } else {
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
            mw_bar += x[i] * cppargs.mw[i];
            eps_mix_num += x[i] * cppargs.mw[i] * cppargs.dielc[i];
        }
        if (mw_bar <= 0.0) {
            throw ValueError("Average molecular weight must be positive for dielc_rule=2.");
        }
        double eps_mix = eps_mix_num / mw_bar;
        for (int i = 0; i < ncomp; i++) {
            deps_dx[i] = (cppargs.mw[i] / mw_bar) * (cppargs.dielc[i] - eps_mix);
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
            mw_sol += x[idx] * cppargs.mw[idx];
            eps_sol_num += x[idx] * cppargs.mw[idx] * cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.");
        }
        double eps_sol_w = eps_sol_num / mw_sol;
        double x_sol = 0.0;
        for (int idx : idx_sol) x_sol += x[idx];
        for (int idx : idx_sol) {
            deps_dx[idx] = eps_sol_w + x_sol * (cppargs.mw[idx] / mw_sol) * (cppargs.dielc[idx] - eps_sol_w);
        }
        for (int idx : idx_ion) {
            deps_dx[idx] = cppargs.dielc[idx];
        }
        return deps_dx;
    }
    if (rule == 9) {
        vector<int> idx_sol;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=9 requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx] * cppargs.mw[idx];
            eps_sol_num += x[idx] * cppargs.mw[idx] * cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=9.");
        }
        double eps_sf = eps_sol_num / mw_sol;
        for (int idx : idx_sol) {
            deps_dx[idx] = (cppargs.mw[idx] / mw_sol) * (cppargs.dielc[idx] - eps_sf);
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
            mw_sol += x[idx] * cppargs.mw[idx];
            eps_sf_num += x[idx] * cppargs.mw[idx] * cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=4.");
        }
        double eps_sf = eps_sf_num / mw_sol;
        double x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        double den = 1.0 + alpha * x_ion;
        for (int idx : idx_sol) {
            deps_dx[idx] = (1.0 / den) * (cppargs.mw[idx] / mw_sol) * (cppargs.dielc[idx] - eps_sf);
        }
        for (int idx : idx_ion) {
            deps_dx[idx] = -alpha * eps_sf / (den * den);
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
    throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8, 9.");
}

vector<double> dielectric_derivative_rule_fd_cpp(int rule, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    vector<double> deps_dx(ncomp, 0.0);
    double f0 = dielectric_constant_rule_cpp(rule, x, cppargs);
    for (int i = 0; i < ncomp; i++) {
        double h = 1e-6 * std::max(1.0, std::abs(x[i]));
        vector<double> xp = x;
        xp[i] += h;
        double fp = dielectric_constant_rule_cpp(rule, xp, cppargs);
        if (x[i] - h >= 0.0) {
            vector<double> xm = x;
            xm[i] -= h;
            double fm = dielectric_constant_rule_cpp(rule, xm, cppargs);
            deps_dx[i] = (fp - fm) / (2.0 * h);
        } else {
            deps_dx[i] = (fp - f0) / h;
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
        vector<AutoDual> x_dual(ncomp, make_autodiff_scalar(0.0, 0.0));
        for (int j = 0; j < ncomp; ++j) {
            x_dual[j] = make_autodiff_scalar(x[j], (i == j) ? 1.0 : 0.0);
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
    if ((cppargs.dielc_diff_mode == 0 || cppargs.dielc_diff_mode == 3) && cppargs.dielc_rule != 8) {
        state.deps_dx = dielectric_derivative_rule_cpp(cppargs.dielc_rule, x, cppargs);
    } else if (cppargs.dielc_diff_mode == 2 || cppargs.dielc_diff_mode == 3) {
        state.deps_dx = dielectric_derivative_rule_ad_cpp(cppargs.dielc_rule, x, cppargs);
    } else {
        state.deps_dx = dielectric_derivative_rule_fd_cpp(cppargs.dielc_rule, x, cppargs);
    }
    return state;
}

vector<double> contribution_dadx_autodiff_cpp(AresContributionKind kind, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (kind == AresContributionKind::ASSOC) {
        throw ValueError("autodiff differential_mode is not implemented for this contribution yet.");
    }
    if (kind == AresContributionKind::BORN && cppargs.born_model == 2) {
        throw ValueError("autodiff differential_mode is not implemented for the SSM/DS Born composition derivative yet.");
    }

    vector<double> dadx(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        vector<AutoDual> x_dual(ncomp, make_autodiff_scalar(0.0, 0.0));
        for (int j = 0; j < ncomp; ++j) {
            x_dual[j] = make_autodiff_scalar(x[j], (i == j) ? 1.0 : 0.0);
        }

        AutoDual value = make_autodiff_scalar(0.0, 0.0);
        if (kind == AresContributionKind::HC || kind == AresContributionKind::DISP) {
            AutodiffMixtureState thermo = mixture_state_autodiff_cpp(t, rho, x_dual, cppargs);
            AutodiffHardChainState hc_state = hard_chain_state_autodiff_cpp(thermo.den, thermo.d, x_dual, cppargs);
            if (kind == AresContributionKind::HC) {
                AutoDual ares_hs = 1.0 / hc_state.zeta[0] * (
                    3.0 * hc_state.zeta[1] * hc_state.zeta[2] / (1.0 - hc_state.zeta[3])
                    + scalar_pow(hc_state.zeta[2], 3) / (hc_state.zeta[3] * scalar_pow(1.0 - hc_state.zeta[3], 2))
                    + (scalar_pow(hc_state.zeta[2], 3) / scalar_pow(hc_state.zeta[3], 2) - hc_state.zeta[0]) * scalar_log(1.0 - hc_state.zeta[3])
                );
                AutoDual log_sum = make_autodiff_scalar(0.0, 0.0);
                for (int k = 0; k < ncomp; ++k) {
                    log_sum += x_dual[k] * (cppargs.m[k] - 1.0) * scalar_log(hc_state.ghs[k * ncomp + k]);
                }
                value = thermo.m_avg * ares_hs - log_sum;
            } else {
                AutodiffDispersionState dispersion = dispersion_state_autodiff_cpp(thermo.m_avg, hc_state.eta);
                value = -2.0 * PI * thermo.den * dispersion.I1 * thermo.m2es3
                    - PI * thermo.den * thermo.m_avg * dispersion.C1 * dispersion.I2 * thermo.m2e2s3;
            }
        } else if (kind == AresContributionKind::ION) {
            AutoDual q2_sum = make_autodiff_scalar(0.0, 0.0);
            for (int k = 0; k < ncomp; ++k) {
                q2_sum += x_dual[k] * cppargs.z[k] * cppargs.z[k];
            }
            AutoDual eps = dielectric_constant_rule_autodiff_cpp(cppargs.dielc_rule, x_dual, cppargs);
            AutoDual kappa = scalar_sqrt((rho * N_AV / 1.0e30) * E_CHRG * E_CHRG / kb / t / perm_vac * q2_sum / eps);
            AutoDual chi_sum = make_autodiff_scalar(0.0, 0.0);
            for (int k = 0; k < ncomp; ++k) {
                double d_k = ion_diameter_cpp(k, t, cppargs);
                AutoDual ka = kappa * d_k;
                AutoDual chi = 3.0 / scalar_pow(ka, 3) * (1.5 + scalar_log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * scalar_pow(1.0 + ka, 2));
                chi_sum += x_dual[k] * cppargs.z[k] * cppargs.z[k] * chi;
            }
            double K0 = E_CHRG * E_CHRG / (12.0 * PI * kb * t * perm_vac);
            value = -K0 * kappa / eps * chi_sum;
        } else if (kind == AresContributionKind::BORN) {
            AutoDual eps = (cppargs.born_eps_mode == 1)
                ? reference_solvent_dielectric_constant_ad_cpp(x_dual, cppargs)
                : dielectric_constant_rule_autodiff_cpp(cppargs.dielc_rule, x_dual, cppargs);
            AutoDual charge_radius_sum = make_autodiff_scalar(0.0, 0.0);
            for (int k = 0; k < ncomp; ++k) {
                if (is_ion_species(cppargs, k)) {
                    charge_radius_sum += x_dual[k] * cppargs.z[k] * cppargs.z[k] / ion_born_radius_cpp(k, t, cppargs);
                }
            }
            const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
            value = -Kborn * (1.0 - 1.0 / eps) * charge_radius_sum;
        }

        dadx[i] = scalar_derivative(value);
        if (!std::isfinite(dadx[i])) {
            throw ValueError("Non-finite contribution autodiff derivative.");
        }
    }
    return dadx;
}

vector<double> contribution_dadx_fd_cpp(AresContributionKind kind, double t, double rho, const vector<double> &x, const add_args &cppargs, double a0) {
    int ncomp = static_cast<int>(x.size());
    vector<double> dadx(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        double h = 1e-6 * std::max(1.0, std::abs(x[i]));
        vector<double> xp = x;
        xp[i] += h;
        double fp = ares_contribution_value_cpp(ares_contributions_cpp(t, rho, xp, cppargs), kind);
        if (x[i] - h >= 0.0) {
            vector<double> xm = x;
            xm[i] -= h;
            double fm = ares_contribution_value_cpp(ares_contributions_cpp(t, rho, xm, cppargs), kind);
            dadx[i] = (fp - fm) / (2.0 * h);
        } else {
            dadx[i] = (fp - a0) / h;
        }
        if (!std::isfinite(dadx[i])) {
            throw ValueError("Non-finite contribution finite-difference derivative.");
        }
    }
    return dadx;
}
