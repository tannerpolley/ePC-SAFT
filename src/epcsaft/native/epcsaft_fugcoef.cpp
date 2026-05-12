#include "epcsaft_core_internal.h"
#include "autodiff/ad_scalar.h"

using namespace thermo_detail;

namespace fugcoef_detail {

template <typename Scalar>
struct SupportedMixtureStateScalar {
    vector<double> d;
    Scalar den = scalar_constant<Scalar>(0.0);
    Scalar m_avg = scalar_constant<Scalar>(0.0);
    Scalar m2es3 = scalar_constant<Scalar>(0.0);
    Scalar m2e2s3 = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
struct SupportedHardChainStateScalar {
    vector<Scalar> zeta;
    vector<Scalar> ghs;
    Scalar eta = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
struct SupportedDispersionStateScalar {
    std::array<Scalar, 7> a{};
    std::array<Scalar, 7> b{};
    Scalar I1 = scalar_constant<Scalar>(0.0);
    Scalar I2 = scalar_constant<Scalar>(0.0);
    Scalar dEtaI1_deta = scalar_constant<Scalar>(0.0);
    Scalar dEtaI2_deta = scalar_constant<Scalar>(0.0);
    Scalar C1 = scalar_constant<Scalar>(0.0);
    Scalar C2 = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
static SupportedMixtureStateScalar<Scalar> mixture_state_scalar_cpp(
    double t,
    const Scalar& rho,
    const vector<Scalar>& x,
    const add_args& cppargs
) {
    const int ncomp = static_cast<int>(x.size());
    SupportedMixtureStateScalar<Scalar> state;
    state.d.assign(static_cast<size_t>(ncomp), 0.0);
    state.den = rho * (N_AV / 1.0e30);
    for (int i = 0; i < ncomp; ++i) {
        state.d[static_cast<size_t>(i)] = cppargs.s[static_cast<size_t>(i)]
            * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[static_cast<size_t>(i)] / t));
        if (!cppargs.z.empty() && std::abs(cppargs.z[static_cast<size_t>(i)]) > 1.0e-12) {
            state.d[static_cast<size_t>(i)] = parameter_setup_detail::ion_diameter_cpp(i, t, cppargs);
        }
        state.m_avg += x[static_cast<size_t>(i)] * cppargs.m[static_cast<size_t>(i)];
    }
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            const double sigma_ij = parameter_setup_detail::pair_sigma_cpp(static_cast<size_t>(idx), i, j, cppargs);
            const double epsilon_ij = parameter_setup_detail::pair_epsilon_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.m2es3 += x[static_cast<size_t>(i)] * x[static_cast<size_t>(j)]
                * cppargs.m[static_cast<size_t>(i)] * cppargs.m[static_cast<size_t>(j)] * epsilon_ij / t
                * scalar_pow(sigma_ij, 3.0);
            state.m2e2s3 += x[static_cast<size_t>(i)] * x[static_cast<size_t>(j)]
                * cppargs.m[static_cast<size_t>(i)] * cppargs.m[static_cast<size_t>(j)] * scalar_pow(epsilon_ij / t, 2.0)
                * scalar_pow(sigma_ij, 3.0);
        }
    }
    return state;
}

template <typename Scalar>
static SupportedHardChainStateScalar<Scalar> hard_chain_state_scalar_cpp(
    const Scalar& den,
    const vector<double>& d,
    const vector<Scalar>& x,
    const add_args& cppargs
) {
    const int ncomp = static_cast<int>(x.size());
    SupportedHardChainStateScalar<Scalar> state;
    state.zeta.assign(4, scalar_constant<Scalar>(0.0));
    state.ghs.assign(static_cast<size_t>(ncomp * ncomp), scalar_constant<Scalar>(0.0));
    for (int k = 0; k < 4; ++k) {
        Scalar summ = scalar_constant<Scalar>(0.0);
        for (int j = 0; j < ncomp; ++j) {
            summ += x[static_cast<size_t>(j)] * cppargs.m[static_cast<size_t>(j)] * scalar_pow(d[static_cast<size_t>(j)], k);
        }
        state.zeta[static_cast<size_t>(k)] = PI / 6.0 * den * summ;
    }
    state.eta = state.zeta[3];
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            const double pair_diameter = parameter_setup_detail::pair_diameter_cpp(d[static_cast<size_t>(i)], d[static_cast<size_t>(j)]);
            state.ghs[static_cast<size_t>(idx)] = 1.0 / (1.0 - state.zeta[3])
                + pair_diameter * 3.0 * state.zeta[2] / scalar_pow(1.0 - state.zeta[3], 2)
                + scalar_pow(pair_diameter, 2.0) * 2.0 * state.zeta[2] * state.zeta[2] / scalar_pow(1.0 - state.zeta[3], 3);
        }
    }
    return state;
}

template <typename Scalar>
static SupportedDispersionStateScalar<Scalar> dispersion_state_scalar_cpp(const Scalar& m_avg, const Scalar& eta) {
    SupportedDispersionStateScalar<Scalar> state;
    const Scalar c1 = (m_avg - scalar_constant<Scalar>(1.0)) / m_avg;
    const Scalar c2 = (m_avg - scalar_constant<Scalar>(2.0)) / m_avg;
    for (size_t i = 0; i < state.a.size(); ++i) {
        state.a[i] = kDispersionA0[i] + c1 * kDispersionA1[i] + c1 * c2 * kDispersionA2[i];
        state.b[i] = kDispersionB0[i] + c1 * kDispersionB1[i] + c1 * c2 * kDispersionB2[i];
        state.I1 += state.a[i] * scalar_pow(eta, static_cast<int>(i));
        state.I2 += state.b[i] * scalar_pow(eta, static_cast<int>(i));
        state.dEtaI1_deta += state.a[i] * static_cast<double>(i + 1) * scalar_pow(eta, static_cast<int>(i));
        state.dEtaI2_deta += state.b[i] * static_cast<double>(i + 1) * scalar_pow(eta, static_cast<int>(i));
    }
    state.C1 = 1.0 / (
        1.0
        + m_avg * (8.0 * eta - 2.0 * eta * eta) / scalar_pow(1.0 - eta, 4)
        + (1.0 - m_avg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * scalar_pow(eta, 3) - 2.0 * scalar_pow(eta, 4))
            / scalar_pow((1.0 - eta) * (scalar_constant<Scalar>(2.0) - eta), 2)
    );
    state.C2 = -state.C1 * state.C1 * (
        m_avg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / scalar_pow(1.0 - eta, 5)
        + (1.0 - m_avg) * (2.0 * scalar_pow(eta, 3) + 12.0 * eta * eta - 48.0 * eta + 40.0)
            / scalar_pow((1.0 - eta) * (scalar_constant<Scalar>(2.0) - eta), 3)
    );
    return state;
}

template <typename Scalar>
static Scalar dadrho_hs_scalar_cpp(const vector<Scalar>& zeta) {
    return zeta[3] / (1.0 - zeta[3])
        + 3.0 * zeta[1] * zeta[2] / zeta[0] / scalar_pow(1.0 - zeta[3], 2)
        + (3.0 * scalar_pow(zeta[2], 3) - zeta[3] * scalar_pow(zeta[2], 3)) / zeta[0] / scalar_pow(1.0 - zeta[3], 3);
}

template <typename Scalar>
static Scalar hs_contact_density_derivative_scalar_cpp(double pair_diameter, const Scalar& zeta2, const Scalar& zeta3) {
    return zeta3 / scalar_pow(1.0 - zeta3, 2)
        + pair_diameter * (3.0 * zeta2 / scalar_pow(1.0 - zeta3, 2) + 6.0 * zeta2 * zeta3 / scalar_pow(1.0 - zeta3, 3))
        + scalar_pow(pair_diameter, 2.0)
            * (4.0 * zeta2 * zeta2 / scalar_pow(1.0 - zeta3, 3) + 6.0 * zeta2 * zeta2 * zeta3 / scalar_pow(1.0 - zeta3, 4));
}

template <typename Scalar>
static Scalar supported_dielectric_constant_scalar_cpp(const vector<Scalar>& x, const add_args& cppargs) {
    if (cppargs.dielc_rule == 0) {
        return scalar_constant<Scalar>(*std::max_element(cppargs.dielc.begin(), cppargs.dielc.end()));
    }
    Scalar eps = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        eps += x[static_cast<size_t>(i)] * cppargs.dielc[static_cast<size_t>(i)];
    }
    return eps;
}

template <typename Scalar>
static Scalar supported_reference_solvent_dielectric_scalar_cpp(const vector<Scalar>& x, const add_args& cppargs) {
    if (cppargs.z.size() != x.size()) {
        return supported_dielectric_constant_scalar_cpp(x, cppargs);
    }
    Scalar x_sol = scalar_constant<Scalar>(0.0);
    Scalar eps_sol_num = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        if (std::abs(cppargs.z[static_cast<size_t>(i)]) <= 1.0e-12) {
            x_sol += x[static_cast<size_t>(i)];
            eps_sol_num += x[static_cast<size_t>(i)] * cppargs.dielc[static_cast<size_t>(i)];
        }
    }
    if (!(scalar_value(x_sol) > 0.0)) {
        return supported_dielectric_constant_scalar_cpp(x, cppargs);
    }
    return eps_sol_num / x_sol;
}

template <typename Scalar>
static Scalar supported_born_dielectric_scalar_cpp(const vector<Scalar>& x, const add_args& cppargs) {
    if (cppargs.born_eps_mode == 1) {
        return supported_reference_solvent_dielectric_scalar_cpp(x, cppargs);
    }
    return supported_dielectric_constant_scalar_cpp(x, cppargs);
}

static vector<double> supported_born_dielectric_derivative_cpp(const vector<double>& x, const add_args& cppargs) {
    if (cppargs.born_eps_mode == 1) {
        return reference_solvent_dielectric_derivative_cpp(x, cppargs);
    }
    return dielectric_diff_cpp(x, cppargs);
}

template <typename Scalar>
static Scalar active_born_radius_scalar_cpp(
    int i,
    double t,
    const add_args& cppargs,
    const std::string* parameter_kind = nullptr,
    int component_index = -1,
    const Scalar* active_value = nullptr
) {
    const bool active_born_radius =
        parameter_kind != nullptr
        && active_value != nullptr
        && (*parameter_kind == "born_radius" || *parameter_kind == "born_diameter")
        && component_index == i;
    if (!is_ion_species(cppargs, i)) {
        if (active_born_radius && cppargs.d_born.size() > static_cast<size_t>(i) && cppargs.d_born[static_cast<size_t>(i)] > 0.0) {
            return *active_value;
        }
        if (cppargs.d_born.size() > static_cast<size_t>(i) && cppargs.d_born[static_cast<size_t>(i)] > 0.0) {
            return scalar_constant<Scalar>(cppargs.d_born[static_cast<size_t>(i)]);
        }
        return scalar_constant<Scalar>(cppargs.s[static_cast<size_t>(i)]);
    }
    const int mode = cppargs.d_born_mode;
    if (mode == 0) {
        return scalar_constant<Scalar>(cppargs.s[static_cast<size_t>(i)]);
    }
    if (mode == 1) {
        return scalar_constant<Scalar>(cppargs.s[static_cast<size_t>(i)] * (1.0 - 0.12));
    }
    if (mode == 2) {
        return scalar_constant<Scalar>(
            cppargs.s[static_cast<size_t>(i)] * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[static_cast<size_t>(i)] / t))
        );
    }
    if (mode == 3) {
        if (active_born_radius) {
            return *active_value;
        }
        return scalar_constant<Scalar>(parameter_setup_detail::ion_born_radius_cpp(i, t, cppargs));
    }
    throw ValueError("Unknown d_Born_mode. Supported values are 0, 1, 2, 3.");
}

template <typename Scalar>
struct SupportedBornShellStateScalar {
    vector<Scalar> d_born;
    vector<Scalar> D;
    vector<Scalar> ddelta_prefac;
    vector<Scalar> f_k;
    vector<Scalar> bracket;
    Scalar sum_bracket = scalar_constant<Scalar>(0.0);
    Scalar sum_gap = scalar_constant<Scalar>(0.0);
    Scalar sum_dpref_over_D2 = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
static SupportedBornShellStateScalar<Scalar> supported_born_shell_state_scalar_cpp(
    const vector<Scalar>& x,
    const add_args& cppargs,
    double t,
    const Scalar& eps_r,
    double eps_r_ion,
    const std::string* parameter_kind = nullptr,
    int component_index = -1,
    const Scalar* active_value = nullptr
) {
    const bool use_ssm = (cppargs.born_solvation_shell_model != 0);
    const bool use_ds = (cppargs.born_dielectric_saturation != 0);
    SupportedBornShellStateScalar<Scalar> data;
    const int ncomp = static_cast<int>(x.size());
    data.d_born.assign(static_cast<size_t>(ncomp), scalar_constant<Scalar>(1.0));
    data.D.assign(static_cast<size_t>(ncomp), scalar_constant<Scalar>(1.0));
    data.ddelta_prefac.assign(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    data.f_k.assign(static_cast<size_t>(ncomp), scalar_constant<Scalar>(1.0));
    data.bracket.assign(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));

    Scalar f_mix = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < ncomp; ++i) {
        const bool is_ion = is_ion_species(cppargs, i);
        Scalar fi = scalar_constant<Scalar>(1.0);
        const bool active_f_solv =
            parameter_kind != nullptr
            && active_value != nullptr
            && (*parameter_kind == "solvation_factor" || *parameter_kind == "f_solv")
            && component_index == i
            && !is_ion;
        if (active_f_solv) {
            fi = *active_value;
        } else if (!is_ion && cppargs.f_solv.size() > static_cast<size_t>(i)) {
            fi = scalar_constant<Scalar>(cppargs.f_solv[static_cast<size_t>(i)]);
        }
        data.f_k[static_cast<size_t>(i)] = fi;
        f_mix += x[static_cast<size_t>(i)] * fi;
        data.d_born[static_cast<size_t>(i)] = active_born_radius_scalar_cpp(
            i,
            t,
            cppargs,
            parameter_kind,
            component_index,
            active_value
        );
        if (is_ion) {
            data.ddelta_prefac[static_cast<size_t>(i)] =
                data.d_born[static_cast<size_t>(i)] / scalar_constant<Scalar>(std::abs(cppargs.z[static_cast<size_t>(i)]));
        }
    }

    const Scalar inv_eps = scalar_constant<Scalar>(1.0) / eps_r;
    const Scalar shell_coeff = scalar_constant<Scalar>(1.0 / eps_r_ion) - inv_eps;
    for (int i = 0; i < ncomp; ++i) {
        const bool is_ion = std::abs(cppargs.z[static_cast<size_t>(i)]) > 1.0e-12;
        if (!is_ion) {
            data.D[static_cast<size_t>(i)] = data.d_born[static_cast<size_t>(i)];
            continue;
        }
        const Scalar delta_di = use_ssm
            ? ((f_mix - scalar_constant<Scalar>(1.0)) * data.ddelta_prefac[static_cast<size_t>(i)])
            : scalar_constant<Scalar>(0.0);
        data.D[static_cast<size_t>(i)] = data.d_born[static_cast<size_t>(i)] + delta_di;
        const Scalar invD = scalar_constant<Scalar>(1.0) / data.D[static_cast<size_t>(i)];
        const Scalar gap = scalar_constant<Scalar>(1.0) / data.d_born[static_cast<size_t>(i)] - invD;
        const Scalar base_term = (scalar_constant<Scalar>(1.0) - inv_eps) * invD;
        const Scalar ds_term = use_ds
            ? ((scalar_constant<Scalar>(1.0) - scalar_constant<Scalar>(1.0 / eps_r_ion)) * gap)
            : scalar_constant<Scalar>(0.0);
        data.bracket[static_cast<size_t>(i)] = base_term + ds_term;
        const double z2 = cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)];
        data.sum_bracket += x[static_cast<size_t>(i)] * scalar_constant<Scalar>(z2) * data.bracket[static_cast<size_t>(i)];
        data.sum_gap += x[static_cast<size_t>(i)] * scalar_constant<Scalar>(z2) * gap;
        if (use_ssm) {
            data.sum_dpref_over_D2 += x[static_cast<size_t>(i)] * scalar_constant<Scalar>(z2)
                * data.ddelta_prefac[static_cast<size_t>(i)] * invD * invD;
        }
    }
    return data;
}

static bool lnfug_composition_derivative_supported_cpp(const add_args& cppargs, std::string* reason) {
    if (!cppargs.assoc_num.empty() || !cppargs.assoc_matrix.empty() || !cppargs.k_hb.empty()
        || !cppargs.e_assoc.empty() || !cppargs.vol_a.empty()) {
        if (reason != nullptr) {
            *reason = "native lnphi-composition derivatives currently support nonassociating states only.";
        }
        return false;
    }
    if (cppargs.DH_model == 2) {
        if (reason != nullptr) {
            *reason = "native lnphi-composition derivatives do not support DH_model=2.";
        }
        return false;
    }
    if (cppargs.dielc_rule != 0 && cppargs.dielc_rule != 1) {
        if (reason != nullptr) {
            *reason = "native lnphi-composition derivatives currently support dielc_rule 0 or 1 only.";
        }
        return false;
    }
    return true;
}

static vector<double> exp_vector(const vector<double> &values) {
    vector<double> out(values.size(), 0.0);
    for (int i = 0; i < static_cast<int>(values.size()); ++i) {
        out[i] = std::exp(values[i]);
    }
    return out;
}

// EqID: lnphi_alpha_near_ideal
static double stable_logz_over_zminus1(double Z) {
    double dz = Z - 1.0;
    if (std::abs(dz) < 1e-8) {
        double dz2 = dz * dz;
        double dz3 = dz2 * dz;
        return 1.0 - 0.5 * dz + dz2 / 3.0 - 0.25 * dz3;
    }
    return std::log(Z) / dz;
}

template <typename Scalar>
static Scalar stable_logz_over_zminus1_scalar(const Scalar& Z) {
    const double dz_value = scalar_value(Z) - 1.0;
    if (std::abs(dz_value) < 1e-8) {
        const Scalar dz = Z - scalar_constant<Scalar>(1.0);
        const Scalar dz2 = dz * dz;
        const Scalar dz3 = dz2 * dz;
        return scalar_constant<Scalar>(1.0) - scalar_constant<Scalar>(0.5) * dz
            + dz2 / scalar_constant<Scalar>(3.0) - scalar_constant<Scalar>(0.25) * dz3;
    }
    return scalar_log(Z) / (Z - scalar_constant<Scalar>(1.0));
}

static double lnfug_correction_scale_cpp(const ScalarContributionTerms &z_raw_terms) {
    vector<double> z_terms = {
        z_raw_terms.hc,
        z_raw_terms.disp,
        z_raw_terms.assoc,
        z_raw_terms.ion,
        z_raw_terms.born,
    };
    double z_scale = z_term_scale_cpp(z_terms, z_raw_terms.total);
    double z_weight = stable_logz_over_zminus1(1.0 + z_raw_terms.total);
    return z_scale * z_weight;
}

// EqID: lnphi_alpha
static vector<double> lnfug_contribution_cpp(
    const vector<double> &mu_term,
    double z_value,
    double z_correction_scale
) {
    vector<double> out(mu_term.size(), 0.0);
    for (int i = 0; i < static_cast<int>(mu_term.size()); ++i) {
        out[i] = mu_term[i] - z_value * z_correction_scale;
    }
    return out;
}

static vector<double> lnfug_hc_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.hc, z_raw_terms.hc, z_correction_scale);
}

static vector<double> lnfug_disp_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.disp, z_raw_terms.disp, z_correction_scale);
}

static vector<double> lnfug_assoc_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.assoc, z_raw_terms.assoc, z_correction_scale);
}

static vector<double> lnfug_ion_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.ion, z_raw_terms.ion, z_correction_scale);
}

static vector<double> lnfug_born_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.born, z_raw_terms.born, z_correction_scale);
}

static vector<double> lnfug_total_cpp(
    const vector<double> &lnfug_hc,
    const vector<double> &lnfug_disp,
    const vector<double> &lnfug_assoc,
    const vector<double> &lnfug_ion,
    const vector<double> &lnfug_born
) {
    vector<double> total(lnfug_hc.size(), 0.0);
    for (int i = 0; i < static_cast<int>(total.size()); ++i) {
        total[i] = lnfug_hc[i] + lnfug_disp[i] + lnfug_assoc[i] + lnfug_ion[i] + lnfug_born[i];
    }
    return total;
}

template <typename Scalar>
static vector<Scalar> supported_lnfug_scalar_cpp(double t, const Scalar& rho, const vector<Scalar>& x, const add_args& cppargs) {
    const int ncomp = static_cast<int>(x.size());
    SupportedMixtureStateScalar<Scalar> thermo = mixture_state_scalar_cpp(t, rho, x, cppargs);
    SupportedHardChainStateScalar<Scalar> hc_state = hard_chain_state_scalar_cpp(thermo.den, thermo.d, x, cppargs);
    SupportedDispersionStateScalar<Scalar> disp_state = dispersion_state_scalar_cpp(thermo.m_avg, hc_state.eta);

    Scalar ares_hs = 1.0 / hc_state.zeta[0] * (
        3.0 * hc_state.zeta[1] * hc_state.zeta[2] / (1.0 - hc_state.zeta[3])
        + scalar_pow(hc_state.zeta[2], 3) / (hc_state.zeta[3] * scalar_pow(1.0 - hc_state.zeta[3], 2))
        + (scalar_pow(hc_state.zeta[2], 3) / scalar_pow(hc_state.zeta[3], 2) - hc_state.zeta[0]) * scalar_log(1.0 - hc_state.zeta[3])
    );
    Scalar ares_hc = thermo.m_avg * ares_hs;
    for (int i = 0; i < ncomp; ++i) {
        ares_hc -= x[static_cast<size_t>(i)] * (cppargs.m[static_cast<size_t>(i)] - 1.0)
            * scalar_log(hc_state.ghs[static_cast<size_t>(i * ncomp + i)]);
    }
    Scalar zraw_hc = thermo.m_avg * dadrho_hs_scalar_cpp(hc_state.zeta);
    for (int i = 0; i < ncomp; ++i) {
        const double pair_diameter = parameter_setup_detail::pair_diameter_cpp(
            thermo.d[static_cast<size_t>(i)],
            thermo.d[static_cast<size_t>(i)]
        );
        zraw_hc -= x[static_cast<size_t>(i)] * (cppargs.m[static_cast<size_t>(i)] - 1.0)
            / hc_state.ghs[static_cast<size_t>(i * ncomp + i)]
            * hs_contact_density_derivative_scalar_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3]);
    }
    vector<Scalar> dadx_hc(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    for (int i = 0; i < ncomp; ++i) {
        Scalar hs_base_value = 1.0 / hc_state.zeta[0] * (
            3.0 * hc_state.zeta[1] * hc_state.zeta[2] / (1.0 - hc_state.zeta[3])
            + scalar_pow(hc_state.zeta[2], 3) / (hc_state.zeta[3] * scalar_pow(1.0 - hc_state.zeta[3], 2))
            + (scalar_pow(hc_state.zeta[2], 3) / scalar_pow(hc_state.zeta[3], 2) - hc_state.zeta[0]) * scalar_log(1.0 - hc_state.zeta[3])
        );
        Scalar dzeta0_dx = PI / 6.0 * thermo.den * cppargs.m[static_cast<size_t>(i)];
        Scalar dzeta1_dx = PI / 6.0 * thermo.den * cppargs.m[static_cast<size_t>(i)] * thermo.d[static_cast<size_t>(i)];
        Scalar dzeta2_dx = PI / 6.0 * thermo.den * cppargs.m[static_cast<size_t>(i)] * scalar_pow(thermo.d[static_cast<size_t>(i)], 2);
        Scalar dzeta3_dx = PI / 6.0 * thermo.den * cppargs.m[static_cast<size_t>(i)] * scalar_pow(thermo.d[static_cast<size_t>(i)], 3);
        Scalar dadx_hs = -dzeta0_dx / hc_state.zeta[0] * hs_base_value
            + 1.0 / hc_state.zeta[0] * (
                3.0 * (dzeta1_dx * hc_state.zeta[2] + hc_state.zeta[1] * dzeta2_dx) / (1.0 - hc_state.zeta[3])
                + 3.0 * hc_state.zeta[1] * hc_state.zeta[2] * dzeta3_dx / scalar_pow(1.0 - hc_state.zeta[3], 2)
                + 3.0 * hc_state.zeta[2] * hc_state.zeta[2] * dzeta2_dx / hc_state.zeta[3] / scalar_pow(1.0 - hc_state.zeta[3], 2)
                + scalar_pow(hc_state.zeta[2], 3) * dzeta3_dx * (3.0 * hc_state.zeta[3] - 1.0)
                    / scalar_pow(hc_state.zeta[3], 2) / scalar_pow(1.0 - hc_state.zeta[3], 3)
                + scalar_log(1.0 - hc_state.zeta[3]) * (
                    (3.0 * hc_state.zeta[2] * hc_state.zeta[2] * dzeta2_dx * hc_state.zeta[3]
                        - 2.0 * scalar_pow(hc_state.zeta[2], 3) * dzeta3_dx) / scalar_pow(hc_state.zeta[3], 3)
                    - dzeta0_dx
                )
                + (hc_state.zeta[0] - scalar_pow(hc_state.zeta[2], 3) / scalar_pow(hc_state.zeta[3], 2))
                    * dzeta3_dx / (1.0 - hc_state.zeta[3])
            );
        dadx_hc[static_cast<size_t>(i)] = thermo.m_avg * dadx_hs - (cppargs.m[static_cast<size_t>(i)] - 1.0)
            * scalar_log(hc_state.ghs[static_cast<size_t>(i * ncomp + i)]);
    }

    Scalar ares_disp = -2.0 * PI * thermo.den * disp_state.I1 * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * disp_state.C1 * disp_state.I2 * thermo.m2e2s3;
    Scalar zraw_disp = -2.0 * PI * thermo.den * disp_state.dEtaI1_deta * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * (disp_state.C1 * disp_state.dEtaI2_deta + disp_state.C2 * hc_state.eta * disp_state.I2) * thermo.m2e2s3;
    vector<Scalar> dadx_disp(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    for (int i = 0; i < ncomp; ++i) {
        Scalar dzeta3_dx = PI / 6.0 * thermo.den * cppargs.m[static_cast<size_t>(i)] * scalar_pow(thermo.d[static_cast<size_t>(i)], 3);
        Scalar dI1_dx = scalar_constant<Scalar>(0.0);
        Scalar dI2_dx = scalar_constant<Scalar>(0.0);
        Scalar dm2es3_dx = scalar_constant<Scalar>(0.0);
        Scalar dm2e2s3_dx = scalar_constant<Scalar>(0.0);
        for (int l = 0; l < 7; ++l) {
            Scalar daa_dx = cppargs.m[static_cast<size_t>(i)] / thermo.m_avg / thermo.m_avg * kDispersionA1[l]
                + cppargs.m[static_cast<size_t>(i)] / thermo.m_avg / thermo.m_avg
                    * (3.0 - 4.0 / thermo.m_avg) * kDispersionA2[l];
            Scalar db_dx = cppargs.m[static_cast<size_t>(i)] / thermo.m_avg / thermo.m_avg * kDispersionB1[l]
                + cppargs.m[static_cast<size_t>(i)] / thermo.m_avg / thermo.m_avg
                    * (3.0 - 4.0 / thermo.m_avg) * kDispersionB2[l];
            dI1_dx += disp_state.a[static_cast<size_t>(l)] * l * dzeta3_dx * scalar_pow(hc_state.eta, l - 1) + daa_dx * scalar_pow(hc_state.eta, l);
            dI2_dx += disp_state.b[static_cast<size_t>(l)] * l * dzeta3_dx * scalar_pow(hc_state.eta, l - 1) + db_dx * scalar_pow(hc_state.eta, l);
        }
        for (int j = 0; j < ncomp; ++j) {
            const size_t idx = static_cast<size_t>(i * ncomp + j);
            dm2es3_dx += x[static_cast<size_t>(j)] * cppargs.m[static_cast<size_t>(j)]
                * (parameter_setup_detail::pair_epsilon_cpp(idx, i, j, cppargs) / t)
                * scalar_pow(parameter_setup_detail::pair_sigma_cpp(idx, i, j, cppargs), 3);
            dm2e2s3_dx += x[static_cast<size_t>(j)] * cppargs.m[static_cast<size_t>(j)]
                * scalar_pow(parameter_setup_detail::pair_epsilon_cpp(idx, i, j, cppargs) / t, 2.0)
                * scalar_pow(parameter_setup_detail::pair_sigma_cpp(idx, i, j, cppargs), 3);
        }
        dm2es3_dx *= 2.0 * cppargs.m[static_cast<size_t>(i)];
        dm2e2s3_dx *= 2.0 * cppargs.m[static_cast<size_t>(i)];
        Scalar dC1_dx = disp_state.C2 * dzeta3_dx - disp_state.C1 * disp_state.C1 * (
            cppargs.m[static_cast<size_t>(i)] * (8.0 * hc_state.eta - 2.0 * hc_state.eta * hc_state.eta) / scalar_pow(1.0 - hc_state.eta, 4)
            - cppargs.m[static_cast<size_t>(i)] * (20.0 * hc_state.eta - 27.0 * hc_state.eta * hc_state.eta
                + 12.0 * scalar_pow(hc_state.eta, 3) - 2.0 * scalar_pow(hc_state.eta, 4))
                / scalar_pow((1.0 - hc_state.eta) * (2.0 - hc_state.eta), 2)
        );
        dadx_disp[static_cast<size_t>(i)] = -2.0 * PI * thermo.den * (dI1_dx * thermo.m2es3 + disp_state.I1 * dm2es3_dx)
            - PI * thermo.den * ((cppargs.m[static_cast<size_t>(i)] * disp_state.C1 * disp_state.I2
                + thermo.m_avg * dC1_dx * disp_state.I2 + thermo.m_avg * disp_state.C1 * dI2_dx) * thermo.m2e2s3
                + thermo.m_avg * disp_state.C1 * disp_state.I2 * dm2e2s3_dx);
    }

    Scalar ares_ion = scalar_constant<Scalar>(0.0);
    Scalar zraw_ion = scalar_constant<Scalar>(0.0);
    vector<Scalar> dadx_ion(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    Scalar eps = scalar_constant<Scalar>(0.0);
    Scalar q2_sum = scalar_constant<Scalar>(0.0);
    Scalar kappa = scalar_constant<Scalar>(0.0);
    Scalar chi_sum = scalar_constant<Scalar>(0.0);
    Scalar sigma_sum = scalar_constant<Scalar>(0.0);
    vector<Scalar> chi(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    vector<Scalar> sigma_k(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    if (!cppargs.z.empty()) {
        vector<double> x_value(static_cast<size_t>(ncomp), 0.0);
        for (int i = 0; i < ncomp; ++i) {
            x_value[static_cast<size_t>(i)] = scalar_value(x[static_cast<size_t>(i)]);
        }
        const vector<double> deps_dx = dielectric_diff_cpp(x_value, cppargs);
        eps = supported_dielectric_constant_scalar_cpp(x, cppargs);
        for (int i = 0; i < ncomp; ++i) {
            q2_sum += x[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)];
        }
        if (scalar_value(q2_sum) > 0.0) {
            kappa = scalar_sqrt(thermo.den * E_CHRG * E_CHRG / kb / t / perm_vac * q2_sum / eps);
            for (int i = 0; i < ncomp; ++i) {
                Scalar ka = kappa * thermo.d[static_cast<size_t>(i)];
                chi[static_cast<size_t>(i)] = 3.0 / scalar_pow(ka, 3)
                    * (1.5 + scalar_log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * scalar_pow(1.0 + ka, 2));
                sigma_k[static_cast<size_t>(i)] = -2.0 * chi[static_cast<size_t>(i)] + 3.0 / (1.0 + ka);
                chi_sum += x[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * chi[static_cast<size_t>(i)];
                sigma_sum += x[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * sigma_k[static_cast<size_t>(i)];
            }
            const double K0 = E_CHRG * E_CHRG / (12.0 * PI * kb * t * perm_vac);
            ares_ion = -K0 * kappa / eps * chi_sum;
            zraw_ion = -kappa / 24.0 / PI / kb / t / (eps * perm_vac) * sigma_sum * E_CHRG * E_CHRG;
            const double a_const = scalar_value(thermo.den) * E_CHRG * E_CHRG / (kb * t * perm_vac);
            for (int i = 0; i < ncomp; ++i) {
                Scalar d_inv_eps_dx = -deps_dx[static_cast<size_t>(i)] / (eps * eps);
                Scalar dkappa_dx = a_const
                    * (cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] / eps
                        - q2_sum * deps_dx[static_cast<size_t>(i)] / (eps * eps))
                    / (2.0 * kappa);
                Scalar dchi_sum_dx = cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * chi[static_cast<size_t>(i)]
                    + dkappa_dx * (sigma_sum - chi_sum) / kappa;
                dadx_ion[static_cast<size_t>(i)] = -K0 * (
                    (dkappa_dx / eps + kappa * d_inv_eps_dx) * chi_sum
                    + kappa / eps * dchi_sum_dx
                );
            }
        }
    }

    Scalar ares_born = scalar_constant<Scalar>(0.0);
    Scalar zraw_born = scalar_constant<Scalar>(0.0);
    vector<Scalar> dadx_born(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    if (cppargs.born_model == 1) {
        vector<double> x_value(static_cast<size_t>(ncomp), 0.0);
        for (int i = 0; i < ncomp; ++i) {
            x_value[static_cast<size_t>(i)] = scalar_value(x[static_cast<size_t>(i)]);
        }
        const vector<double> deps_dx = supported_born_dielectric_derivative_cpp(x_value, cppargs);
        Scalar eps_born = supported_born_dielectric_scalar_cpp(x, cppargs);
        Scalar charge_radius_sum = scalar_constant<Scalar>(0.0);
        for (int i = 0; i < ncomp; ++i) {
            if (is_ion_species(cppargs, i)) {
                charge_radius_sum += x[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)]
                    / active_born_radius_scalar_cpp<Scalar>(i, t, cppargs);
            }
        }
        const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        ares_born = -Kborn * (1.0 - 1.0 / eps_born) * charge_radius_sum;
        for (int i = 0; i < ncomp; ++i) {
            Scalar ion_part = scalar_constant<Scalar>(0.0);
            if (is_ion_species(cppargs, i)) {
                ion_part = (1.0 - 1.0 / eps_born) * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)]
                    / active_born_radius_scalar_cpp<Scalar>(i, t, cppargs);
            }
            Scalar eps_part = charge_radius_sum * deps_dx[static_cast<size_t>(i)] / (eps_born * eps_born);
            dadx_born[static_cast<size_t>(i)] = -Kborn * (ion_part + eps_part);
        }
    } else if (cppargs.born_model == 2) {
        vector<double> x_value(static_cast<size_t>(ncomp), 0.0);
        for (int i = 0; i < ncomp; ++i) {
            x_value[static_cast<size_t>(i)] = scalar_value(x[static_cast<size_t>(i)]);
        }
        const vector<double> deps_dx = supported_born_dielectric_derivative_cpp(x_value, cppargs);
        const Scalar eps_born = supported_born_dielectric_scalar_cpp(x, cppargs);
        const double eps_r_ion = 8.0;
        const auto shell = supported_born_shell_state_scalar_cpp<Scalar>(x, cppargs, t, eps_born, eps_r_ion);
        const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        ares_born = -Kborn * shell.sum_bracket;

        const Scalar inv_eps2 = scalar_constant<Scalar>(1.0) / (eps_born * eps_born);
        const Scalar shell_coeff = scalar_constant<Scalar>(1.0 / eps_r_ion) - scalar_constant<Scalar>(1.0) / eps_born;
        const bool use_deps = (cppargs.mu_born_comp_dep_rel_perm != 0);
        const bool use_shell_chain = (cppargs.mu_born_comp_dep_delta_d != 0);
        const Scalar deps_multiplier = (cppargs.mu_born_include_sum_term != 0)
            ? shell.sum_gap
            : scalar_constant<Scalar>(1.0);
        for (int k = 0; k < ncomp; ++k) {
            Scalar direct_part = scalar_constant<Scalar>(0.0);
            if (std::abs(cppargs.z[static_cast<size_t>(k)]) > 1.0e-12) {
                direct_part = cppargs.z[static_cast<size_t>(k)] * cppargs.z[static_cast<size_t>(k)]
                    * shell.bracket[static_cast<size_t>(k)];
            }
            const Scalar deps_part = use_deps
                ? deps_multiplier * deps_dx[static_cast<size_t>(k)] * inv_eps2
                : scalar_constant<Scalar>(0.0);
            const Scalar ddelta_part = use_shell_chain
                ? shell_coeff * shell.sum_dpref_over_D2 * shell.f_k[static_cast<size_t>(k)]
                : scalar_constant<Scalar>(0.0);
            dadx_born[static_cast<size_t>(k)] = -Kborn * (direct_part + deps_part + ddelta_part);
        }
    }

    Scalar z_total = zraw_hc + zraw_disp + zraw_ion + zraw_born;
    const Scalar z_scale = stable_logz_over_zminus1_scalar(scalar_constant<Scalar>(1.0) + z_total);
    Scalar sum_x_dadx_hc = scalar_constant<Scalar>(0.0);
    Scalar sum_x_dadx_disp = scalar_constant<Scalar>(0.0);
    Scalar sum_x_dadx_ion = scalar_constant<Scalar>(0.0);
    Scalar sum_x_dadx_born = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < ncomp; ++i) {
        sum_x_dadx_hc += x[static_cast<size_t>(i)] * dadx_hc[static_cast<size_t>(i)];
        sum_x_dadx_disp += x[static_cast<size_t>(i)] * dadx_disp[static_cast<size_t>(i)];
        sum_x_dadx_ion += x[static_cast<size_t>(i)] * dadx_ion[static_cast<size_t>(i)];
        sum_x_dadx_born += x[static_cast<size_t>(i)] * dadx_born[static_cast<size_t>(i)];
    }

    vector<Scalar> lnfug(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    for (int i = 0; i < ncomp; ++i) {
        Scalar mu_hc = ares_hc + zraw_hc + dadx_hc[static_cast<size_t>(i)] - sum_x_dadx_hc;
        Scalar mu_disp = ares_disp + zraw_disp + dadx_disp[static_cast<size_t>(i)] - sum_x_dadx_disp;
        Scalar mu_ion = ares_ion + zraw_ion + dadx_ion[static_cast<size_t>(i)] - sum_x_dadx_ion;
        Scalar mu_born = ares_born + zraw_born + dadx_born[static_cast<size_t>(i)] - sum_x_dadx_born;
        lnfug[static_cast<size_t>(i)] = (mu_hc - zraw_hc * z_scale)
            + (mu_disp - zraw_disp * z_scale)
            + (mu_ion - zraw_ion * z_scale)
            + (mu_born - zraw_born * z_scale);
    }
    return lnfug;
}

}  // namespace fugcoef_detail

// EqID: lnphi_total
// EqID: lnphi_total_sum
FugacityContributionResult fugacity_coefficient_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    ResidualChemicalPotentialResult mu_result = residual_chemical_potential_result_cpp(t, rho, std::move(x), cppargs);
    double z_correction_scale = fugcoef_detail::lnfug_correction_scale_cpp(mu_result.composition.z_raw);

    vector<double> lnfug_hc = fugcoef_detail::lnfug_hc_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_disp = fugcoef_detail::lnfug_disp_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_assoc = fugcoef_detail::lnfug_assoc_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_ion = fugcoef_detail::lnfug_ion_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_born = fugcoef_detail::lnfug_born_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_total = fugcoef_detail::lnfug_total_cpp(lnfug_hc, lnfug_disp, lnfug_assoc, lnfug_ion, lnfug_born);

    FugacityContributionResult result;
    result.mu = mu_result.mu;
    result.composition = mu_result.composition;
    result.lnfugcoef = make_vector_terms(lnfug_hc, lnfug_disp, lnfug_assoc, lnfug_ion, lnfug_born, lnfug_total);
    return result;
}

vector<double> lnfug_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return fugacity_coefficient_result_cpp(t, rho, std::move(x), cppargs).lnfugcoef.total;
}

vector<double> fugcoef_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return fugcoef_detail::exp_vector(lnfug_cpp(t, rho, std::move(x), cppargs));
}

LnfugCompositionDerivativeResult lnfug_composition_derivative_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    LnfugCompositionDerivativeResult result;
    result.rows = static_cast<int>(x.size());
    result.cols = static_cast<int>(x.size());
    result.lnfug.assign(x.size(), std::numeric_limits<double>::quiet_NaN());
    result.dlnfugdx_row_major.assign(x.size() * x.size(), std::numeric_limits<double>::quiet_NaN());
    std::string unsupported_reason;
    if (!fugcoef_detail::lnfug_composition_derivative_supported_cpp(cppargs, &unsupported_reason)) {
        result.finite_difference_fallback_reason = unsupported_reason;
        return result;
    }
#ifndef EPCSAFT_HAS_CPPAD
    result.finite_difference_fallback_reason =
        "backend_unavailable: native lnphi-composition derivatives require a CppAD-enabled build.";
    return result;
#else
    using epcsaft::autodiff::CppADScalar;
    const std::size_t ncomp = x.size();
    std::vector<CppADScalar> independent(ncomp);
    for (std::size_t i = 0; i < ncomp; ++i) {
        independent[i] = x[i];
    }
    CppAD::Independent(independent);
    std::vector<CppADScalar> dependent = fugcoef_detail::supported_lnfug_scalar_cpp(t, CppADScalar(rho), independent, cppargs);
    CppAD::ADFun<double> tape(independent, dependent);
    const std::vector<double> jacobian = tape.Jacobian(x);
    const std::vector<double> canonical_lnfug = lnfug_cpp(t, rho, x, cppargs);
    result.supported = true;
    result.derivative_backend = "autodiff_composition";
    result.finite_difference_fallback_used = false;
    result.lnfug = canonical_lnfug;
    result.dlnfugdx_row_major = jacobian;
    return result;
#endif
}

LnfugDensityDerivativeResult lnfug_density_derivative_result_cpp(double t, double rho, vector<double> x, const add_args& cppargs) {
    LnfugDensityDerivativeResult result;
    result.size = static_cast<int>(x.size());
    result.lnfug.assign(x.size(), std::numeric_limits<double>::quiet_NaN());
    result.dlnfugdrho.assign(x.size(), std::numeric_limits<double>::quiet_NaN());
    std::string unsupported_reason;
    if (!fugcoef_detail::lnfug_composition_derivative_supported_cpp(cppargs, &unsupported_reason)) {
        result.finite_difference_fallback_reason = unsupported_reason;
        return result;
    }
#ifndef EPCSAFT_HAS_CPPAD
    result.finite_difference_fallback_reason =
        "backend_unavailable: native lnphi-density derivatives require a CppAD-enabled build.";
    return result;
#else
    using epcsaft::autodiff::CppADScalar;
    const std::size_t ncomp = x.size();
    std::vector<CppADScalar> independent(1);
    independent[0] = rho;
    CppAD::Independent(independent);
    std::vector<CppADScalar> x_const(ncomp);
    for (std::size_t i = 0; i < ncomp; ++i) {
        x_const[i] = x[i];
    }
    std::vector<CppADScalar> dependent = fugcoef_detail::supported_lnfug_scalar_cpp(t, independent[0], x_const, cppargs);
    CppAD::ADFun<double> tape(independent, dependent);
    const std::vector<double> jacobian = tape.Jacobian(std::vector<double>{rho});
    const std::vector<double> canonical_lnfug = lnfug_cpp(t, rho, x, cppargs);
    result.supported = true;
    result.derivative_backend = "autodiff_density";
    result.finite_difference_fallback_used = false;
    result.lnfug = canonical_lnfug;
    for (std::size_t i = 0; i < ncomp; ++i) {
        result.dlnfugdrho[i] = jacobian[i];
        if (!std::isfinite(result.dlnfugdrho[i])) {
            throw ValueError("Non-finite native lnphi-density derivative.");
        }
    }
    return result;
#endif
}

namespace fugcoef_detail {

static bool lnfug_parameter_derivative_supported_cpp(
    const add_args& cppargs,
    const std::string& parameter_kind,
    int component_index,
    std::string* reason
) {
    std::string unsupported_reason;
    if (!lnfug_composition_derivative_supported_cpp(cppargs, &unsupported_reason)) {
        if (reason != nullptr) {
            *reason = unsupported_reason;
        }
        return false;
    }
    const bool born_radius = (parameter_kind == "born_radius" || parameter_kind == "born_diameter");
    const bool f_solv = (parameter_kind == "solvation_factor" || parameter_kind == "f_solv");
    if (!born_radius && !f_solv) {
        if (reason != nullptr) {
            *reason = "native lnphi-parameter derivatives currently support born_radius and f_solv only.";
        }
        return false;
    }
    if (component_index < 0 || static_cast<size_t>(component_index) >= cppargs.s.size()) {
        if (reason != nullptr) {
            *reason = "native lnphi-parameter derivative component index is out of range.";
        }
        return false;
    }
    if (born_radius) {
        if (is_ion_species(cppargs, component_index)) {
            if (cppargs.d_born_mode != 3) {
                if (reason != nullptr) {
                    *reason = "native lnphi born_radius derivatives for ionic species require d_born_mode=3.";
                }
                return false;
            }
            return true;
        }
        if (cppargs.d_born.size() > static_cast<size_t>(component_index)
            && cppargs.d_born[static_cast<size_t>(component_index)] > 0.0) {
            return true;
        }
        if (reason != nullptr) {
            *reason = "native lnphi born_radius derivatives for neutral species require an explicit positive d_born entry.";
        }
        return false;
    }
    if (is_ion_species(cppargs, component_index)) {
        if (reason != nullptr) {
            *reason = "native lnphi f_solv derivatives require a neutral solvent component.";
        }
        return false;
    }
    if (cppargs.born_model != 2) {
        if (reason != nullptr) {
            *reason = "native lnphi f_solv derivatives currently require born_model=2.";
        }
        return false;
    }
    return true;
}

template <typename Scalar>
static vector<Scalar> born_lnfug_parameter_scalar_cpp(
    double t,
    double rho,
    const vector<double>& x,
    const add_args& cppargs,
    const std::string& parameter_kind,
    int component_index,
    const Scalar& active_value
) {
    const int ncomp = static_cast<int>(x.size());
    vector<Scalar> out(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    vector<Scalar> x_scalar(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    for (int i = 0; i < ncomp; ++i) {
        x_scalar[static_cast<size_t>(i)] = scalar_constant<Scalar>(x[static_cast<size_t>(i)]);
    }

    const vector<double> deps_dx = supported_born_dielectric_derivative_cpp(x, cppargs);
    const Scalar eps_born = supported_born_dielectric_scalar_cpp(x_scalar, cppargs);
    vector<Scalar> dadx_born(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    Scalar ares_born = scalar_constant<Scalar>(0.0);

    if (cppargs.born_model == 1) {
        Scalar charge_radius_sum = scalar_constant<Scalar>(0.0);
        for (int i = 0; i < ncomp; ++i) {
            if (is_ion_species(cppargs, i)) {
                charge_radius_sum += x_scalar[static_cast<size_t>(i)]
                    * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)]
                    / active_born_radius_scalar_cpp<Scalar>(
                        i,
                        t,
                        cppargs,
                        &parameter_kind,
                        component_index,
                        &active_value
                    );
            }
        }
        const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        ares_born = -Kborn * (scalar_constant<Scalar>(1.0) - scalar_constant<Scalar>(1.0) / eps_born) * charge_radius_sum;
        for (int i = 0; i < ncomp; ++i) {
            Scalar ion_part = scalar_constant<Scalar>(0.0);
            if (is_ion_species(cppargs, i)) {
                ion_part = (scalar_constant<Scalar>(1.0) - scalar_constant<Scalar>(1.0) / eps_born)
                    * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)]
                    / active_born_radius_scalar_cpp<Scalar>(
                        i,
                        t,
                        cppargs,
                        &parameter_kind,
                        component_index,
                        &active_value
                    );
            }
            const Scalar eps_part = charge_radius_sum * deps_dx[static_cast<size_t>(i)] / (eps_born * eps_born);
            dadx_born[static_cast<size_t>(i)] = -Kborn * (ion_part + eps_part);
        }
    } else if (cppargs.born_model == 2) {
        const double eps_r_ion = 8.0;
        const auto shell = supported_born_shell_state_scalar_cpp<Scalar>(
            x_scalar,
            cppargs,
            t,
            eps_born,
            eps_r_ion,
            &parameter_kind,
            component_index,
            &active_value
        );
        const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        ares_born = -Kborn * shell.sum_bracket;
        const Scalar inv_eps2 = scalar_constant<Scalar>(1.0) / (eps_born * eps_born);
        const Scalar shell_coeff = scalar_constant<Scalar>(1.0 / eps_r_ion) - scalar_constant<Scalar>(1.0) / eps_born;
        const bool use_deps = (cppargs.mu_born_comp_dep_rel_perm != 0);
        const bool use_shell_chain = (cppargs.mu_born_comp_dep_delta_d != 0);
        const Scalar deps_multiplier = (cppargs.mu_born_include_sum_term != 0)
            ? shell.sum_gap
            : scalar_constant<Scalar>(1.0);
        for (int k = 0; k < ncomp; ++k) {
            Scalar direct_part = scalar_constant<Scalar>(0.0);
            if (std::abs(cppargs.z[static_cast<size_t>(k)]) > 1.0e-12) {
                direct_part = cppargs.z[static_cast<size_t>(k)] * cppargs.z[static_cast<size_t>(k)]
                    * shell.bracket[static_cast<size_t>(k)];
            }
            const Scalar deps_part = use_deps
                ? deps_multiplier * deps_dx[static_cast<size_t>(k)] * inv_eps2
                : scalar_constant<Scalar>(0.0);
            const Scalar ddelta_part = use_shell_chain
                ? shell_coeff * shell.sum_dpref_over_D2 * shell.f_k[static_cast<size_t>(k)]
                : scalar_constant<Scalar>(0.0);
            dadx_born[static_cast<size_t>(k)] = -Kborn * (direct_part + deps_part + ddelta_part);
        }
    }

    Scalar sum_x_dadx_born = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < ncomp; ++i) {
        sum_x_dadx_born += x_scalar[static_cast<size_t>(i)] * dadx_born[static_cast<size_t>(i)];
    }
    for (int i = 0; i < ncomp; ++i) {
        out[static_cast<size_t>(i)] = ares_born + dadx_born[static_cast<size_t>(i)] - sum_x_dadx_born;
    }
    return out;
}

}  // namespace fugcoef_detail

LnfugParameterDerivativeResult lnfug_parameter_derivative_result_cpp(
    double t,
    double rho,
    vector<double> x,
    const add_args& cppargs,
    const std::string& parameter_kind,
    int component_index
) {
    LnfugParameterDerivativeResult result;
    result.size = static_cast<int>(x.size());
    result.dlnfugdtheta.assign(x.size(), std::numeric_limits<double>::quiet_NaN());
    std::string unsupported_reason;
    if (!fugcoef_detail::lnfug_parameter_derivative_supported_cpp(cppargs, parameter_kind, component_index, &unsupported_reason)) {
        result.finite_difference_fallback_reason = unsupported_reason;
        return result;
    }
#ifndef EPCSAFT_HAS_CPPAD
    result.finite_difference_fallback_reason =
        "backend_unavailable: native lnphi-parameter derivatives require a CppAD-enabled build.";
    return result;
#else
    using epcsaft::autodiff::CppADScalar;
    double parameter_value = 0.0;
    if (parameter_kind == "born_radius" || parameter_kind == "born_diameter") {
        parameter_value = fugcoef_detail::active_born_radius_scalar_cpp<double>(component_index, t, cppargs);
    } else {
        parameter_value = (cppargs.f_solv.size() > static_cast<size_t>(component_index))
            ? cppargs.f_solv[static_cast<size_t>(component_index)]
            : 1.0;
    }
    std::vector<CppADScalar> independent(1);
    independent[0] = parameter_value;
    CppAD::Independent(independent);
    std::vector<CppADScalar> dependent = fugcoef_detail::born_lnfug_parameter_scalar_cpp(
        t,
        rho,
        x,
        cppargs,
        parameter_kind,
        component_index,
        independent[0]
    );
    CppAD::ADFun<double> tape(independent, dependent);
    const std::vector<double> jacobian = tape.Jacobian(std::vector<double>{parameter_value});
    result.supported = true;
    result.derivative_backend = "autodiff_parameter";
    result.finite_difference_fallback_used = false;
    for (std::size_t i = 0; i < jacobian.size(); ++i) {
        result.dlnfugdtheta[i] = jacobian[i];
        if (!std::isfinite(result.dlnfugdtheta[i])) {
            throw ValueError("Non-finite native lnphi-parameter derivative.");
        }
    }
    return result;
#endif
}
