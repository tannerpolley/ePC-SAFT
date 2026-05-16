#include "contributions/epcsaft_contrib_internal.h"

namespace {

struct BornParamPartials {
    vector<double> da;
    vector<double> ddadx_row_major;
};

static double born_prefactor_cpp(double t)
{
    return E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
}

static void assign_not_available(BornSSMDSDerivativeResult& result, int ncomp, const std::string& message)
{
    result.supported = false;
    result.backend = "not_available";
    result.message = message;
    result.ncomp = ncomp;
    result.a_born_d_d_born.assign(ncomp, 0.0);
    result.a_born_d_f_solv.assign(ncomp, 0.0);
    result.mu_res_d_d_born_row_major.assign(ncomp * ncomp, 0.0);
    result.mu_res_d_f_solv_row_major.assign(ncomp * ncomp, 0.0);
    result.lnfug_d_d_born_row_major.assign(ncomp * ncomp, 0.0);
    result.lnfug_d_f_solv_row_major.assign(ncomp * ncomp, 0.0);
    result.lngamma_d_d_born_row_major.assign(ncomp * ncomp, 0.0);
    result.lngamma_d_f_solv_row_major.assign(ncomp * ncomp, 0.0);
}

static void initialize_supported_result(BornSSMDSDerivativeResult& result, int ncomp)
{
    result.supported = true;
    result.backend = "analytic";
    result.message = "Liquid-electrolyte SSM+DS Born derivatives are available for d_born and f_solv.";
    result.ncomp = ncomp;
    result.a_born_d_d_born.assign(ncomp, 0.0);
    result.a_born_d_f_solv.assign(ncomp, 0.0);
    result.mu_res_d_d_born_row_major.assign(ncomp * ncomp, 0.0);
    result.mu_res_d_f_solv_row_major.assign(ncomp * ncomp, 0.0);
    result.lnfug_d_d_born_row_major.assign(ncomp * ncomp, 0.0);
    result.lnfug_d_f_solv_row_major.assign(ncomp * ncomp, 0.0);
    result.lngamma_d_d_born_row_major.assign(ncomp * ncomp, 0.0);
    result.lngamma_d_f_solv_row_major.assign(ncomp * ncomp, 0.0);
}

static vector<double> born_q_factors_cpp(const BornSSMDSData& shell, const add_args& args)
{
    int ncomp = static_cast<int>(shell.d_born.size());
    vector<double> q(ncomp, 1.0);
    for (int i = 0; i < ncomp; ++i) {
        if (is_ion_species(args, i) && shell.d_born[i] != 0.0) {
            q[i] = shell.D[i] / shell.d_born[i];
        }
    }
    return q;
}

static BornParamPartials born_d_born_partials_cpp(
    const BornIntermediateState& born_state,
    double t,
    const vector<double>& x,
    const add_args& args
) {
    const int ncomp = static_cast<int>(x.size());
    const double kborn = born_prefactor_cpp(t);
    const double eps = born_state.eps_value;
    const double eps_ion = 8.0;
    const double c_bulk = 1.0 - 1.0 / eps;
    const double c_ion = (args.born_dielectric_saturation != 0) ? (1.0 - 1.0 / eps_ion) : 0.0;
    const double shell_coeff = c_bulk - c_ion;
    const double inv_eps2 = 1.0 / (eps * eps);
    const bool use_deps = (args.mu_born_comp_dep_rel_perm != 0);
    const bool use_deps_sum = use_deps && (args.mu_born_include_sum_term != 0);
    const bool use_shell_chain = (args.mu_born_comp_dep_delta_d != 0);
    const vector<double> q = born_q_factors_cpp(born_state.shell, args);

    BornParamPartials out;
    out.da.assign(ncomp, 0.0);
    out.ddadx_row_major.assign(ncomp * ncomp, 0.0);

    if (args.d_born_mode != 3) {
        return out;
    }

    for (int p = 0; p < ncomp; ++p) {
        if (!is_ion_species(args, p)) {
            continue;
        }
        const double d = born_state.shell.d_born[p];
        const double absz = std::abs(args.z[p]);
        const double z2 = args.z[p] * args.z[p];
        const double qfac = q[p];
        const double dbracket = -c_ion / (d * d) - shell_coeff / (qfac * d * d);
        const double dgap = (1.0 / qfac - 1.0) / (d * d);
        const double dsum_gap = x[p] * z2 * dgap;
        const double dsum_s = x[p] * z2 * (-1.0 / (absz * d * d * qfac * qfac));

        out.da[p] = -kborn * x[p] * z2 * dbracket;

        vector<double> ddadx(ncomp, 0.0);
        for (int k = 0; k < ncomp; ++k) {
            double ddirect = 0.0;
            if (k == p) {
                ddirect = z2 * dbracket;
            }
            const double ddeps = use_deps_sum ? dsum_gap * born_state.deps_dx[k] * inv_eps2 : 0.0;
            const double ddelta = use_shell_chain ? shell_coeff * dsum_s * born_state.shell.f_k[k] : 0.0;
            ddadx[k] = -kborn * (ddirect + ddeps + ddelta);
        }
        const double dsum_x_dadx = std::inner_product(x.begin(), x.end(), ddadx.begin(), 0.0);
        for (int k = 0; k < ncomp; ++k) {
            out.ddadx_row_major[p * ncomp + k] = out.da[p] + ddadx[k] - dsum_x_dadx;
        }
    }
    return out;
}

static BornParamPartials born_f_solv_partials_cpp(
    const BornIntermediateState& born_state,
    double t,
    const vector<double>& x,
    const add_args& args
) {
    const int ncomp = static_cast<int>(x.size());
    const double kborn = born_prefactor_cpp(t);
    const double eps = born_state.eps_value;
    const double eps_ion = 8.0;
    const double c_bulk = 1.0 - 1.0 / eps;
    const double c_ion = (args.born_dielectric_saturation != 0) ? (1.0 - 1.0 / eps_ion) : 0.0;
    const double shell_coeff = c_bulk - c_ion;
    const double inv_eps2 = 1.0 / (eps * eps);
    const bool use_ssm = (args.born_solvation_shell_model != 0);
    const bool use_deps = (args.mu_born_comp_dep_rel_perm != 0);
    const bool use_deps_sum = use_deps && (args.mu_born_include_sum_term != 0);
    const bool use_shell_chain = (args.mu_born_comp_dep_delta_d != 0);
    const vector<double> q = born_q_factors_cpp(born_state.shell, args);

    BornParamPartials out;
    out.da.assign(ncomp, 0.0);
    out.ddadx_row_major.assign(ncomp * ncomp, 0.0);

    if (!use_ssm) {
        return out;
    }

    for (int p = 0; p < ncomp; ++p) {
        if (is_ion_species(args, p)) {
            continue;
        }
        const double dfmix = x[p];
        double dsum_bracket = 0.0;
        double dsum_gap = 0.0;
        double dsum_s = 0.0;
        vector<double> dbracket(ncomp, 0.0);

        for (int i = 0; i < ncomp; ++i) {
            if (!is_ion_species(args, i)) {
                continue;
            }
            const double d = born_state.shell.d_born[i];
            const double absz = std::abs(args.z[i]);
            const double z2 = args.z[i] * args.z[i];
            const double invD = 1.0 / born_state.shell.D[i];
            const double dq = dfmix / absz;
            const double dinvD = -invD * invD * d * dq;
            const double dgap = -dinvD;
            dbracket[i] = c_bulk * dinvD + c_ion * dgap;
            const double ds_term = -2.0 * dq / (absz * d * q[i] * q[i] * q[i]);

            dsum_bracket += x[i] * z2 * dbracket[i];
            dsum_gap += x[i] * z2 * dgap;
            dsum_s += x[i] * z2 * ds_term;
        }

        out.da[p] = -kborn * dsum_bracket;

        vector<double> ddadx(ncomp, 0.0);
        for (int k = 0; k < ncomp; ++k) {
            double ddirect = 0.0;
            if (is_ion_species(args, k)) {
                ddirect = args.z[k] * args.z[k] * dbracket[k];
            }
            const double ddeps = use_deps_sum ? dsum_gap * born_state.deps_dx[k] * inv_eps2 : 0.0;
            double ddelta = use_shell_chain ? shell_coeff * dsum_s * born_state.shell.f_k[k] : 0.0;
            if (use_shell_chain && k == p) {
                ddelta += shell_coeff * born_state.shell.sum_dpref_over_D2;
            }
            ddadx[k] = -kborn * (ddirect + ddeps + ddelta);
        }
        const double dsum_x_dadx = std::inner_product(x.begin(), x.end(), ddadx.begin(), 0.0);
        for (int k = 0; k < ncomp; ++k) {
            out.ddadx_row_major[p * ncomp + k] = out.da[p] + ddadx[k] - dsum_x_dadx;
        }
    }
    return out;
}

static void copy_born_partials_to_result(
    const BornParamPartials& partials,
    vector<double>& da_out,
    vector<double>& mu_out,
    vector<double>& lnfug_out,
    vector<double>& lngamma_out
) {
    da_out = partials.da;
    mu_out = partials.ddadx_row_major;
    lnfug_out = partials.ddadx_row_major;
    lngamma_out = partials.ddadx_row_major;
}

}  // namespace

BornSSMDSDerivativeResult born_ssmds_liquid_derivatives_cpp(
    double t,
    double,
    int phase,
    vector<double> x,
    const add_args &cppargs
) {
    const int ncomp = static_cast<int>(x.size());
    BornSSMDSDerivativeResult result;

    if (phase != 0) {
        assign_not_available(result, ncomp, "SSM+DS Born derivatives are liquid-electrolyte only");
        return result;
    }
    if (cppargs.z.empty() || cppargs.born_model != 2) {
        assign_not_available(result, ncomp, "not_available: SSM+DS Born derivatives require the liquid SSM/DS Born model.");
        return result;
    }
    if (cppargs.born_eps_mode == 1) {
        assign_not_available(result, ncomp, "not_available: SSM+DS Born parameter derivatives currently use mixture dielectric routing.");
        return result;
    }

    BornIntermediateState born_state = born_intermediate_state_cpp(t, x, cppargs, false);
    initialize_supported_result(result, ncomp);

    const BornParamPartials d_born = born_d_born_partials_cpp(born_state, t, x, cppargs);
    const BornParamPartials f_solv = born_f_solv_partials_cpp(born_state, t, x, cppargs);
    copy_born_partials_to_result(
        d_born,
        result.a_born_d_d_born,
        result.mu_res_d_d_born_row_major,
        result.lnfug_d_d_born_row_major,
        result.lngamma_d_d_born_row_major
    );
    copy_born_partials_to_result(
        f_solv,
        result.a_born_d_f_solv,
        result.mu_res_d_f_solv_row_major,
        result.lnfug_d_f_solv_row_major,
        result.lngamma_d_f_solv_row_major
    );
    return result;
}
