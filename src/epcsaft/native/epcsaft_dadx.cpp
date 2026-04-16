#include "epcsaft_core_internal.h"
#include "epcsaft_autodiff_internal.h"
#include "contributions/epcsaft_contrib_internal.h"

using namespace thermo_detail;

namespace {

struct AutodiffMixtureState {
    vector<double> d;
    vector<double> s_ij;
    vector<double> e_ij;
    double den = 0.0;
    AutoDual m_avg = 0.0;
    AutoDual m2es3 = 0.0;
    AutoDual m2e2s3 = 0.0;
};

struct AutodiffDispersionState {
    std::array<double, 7> a{};
    std::array<double, 7> b{};
    AutoDual I1 = 0.0;
    AutoDual I2 = 0.0;
    AutoDual C1 = 0.0;
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

    state.m_avg = 0.0;
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

vector<double> contribution_dadx_autodiff_cpp(AresContributionKind kind, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (kind == AresContributionKind::ASSOC || kind == AresContributionKind::POLAR) {
        throw ValueError("autodiff differential_mode is not implemented for this contribution yet.");
    }
    if (kind == AresContributionKind::BORN && cppargs.born_model == 2) {
        throw ValueError("autodiff differential_mode is not implemented for the SSM/DS Born composition derivative yet.");
    }

    vector<double> dadx(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        vector<AutoDual> x_dual(ncomp, AutoDual(0.0, 0.0));
        for (int j = 0; j < ncomp; ++j) {
            x_dual[j] = AutoDual(x[j], (i == j) ? 1.0 : 0.0);
        }

        AutoDual value = 0.0;
        if (kind == AresContributionKind::HC || kind == AresContributionKind::DISP) {
            AutodiffMixtureState thermo = mixture_state_autodiff_cpp(t, rho, x_dual, cppargs);
            AutodiffHardChainState hc_state = hard_chain_state_autodiff_cpp(thermo.den, thermo.d, x_dual, cppargs);
            if (kind == AresContributionKind::HC) {
                AutoDual ares_hs = 1.0 / hc_state.zeta[0] * (
                    3.0 * hc_state.zeta[1] * hc_state.zeta[2] / (1.0 - hc_state.zeta[3])
                    + scalar_pow(hc_state.zeta[2], 3) / (hc_state.zeta[3] * scalar_pow(1.0 - hc_state.zeta[3], 2))
                    + (scalar_pow(hc_state.zeta[2], 3) / scalar_pow(hc_state.zeta[3], 2) - hc_state.zeta[0]) * scalar_log(1.0 - hc_state.zeta[3])
                );
                AutoDual log_sum = 0.0;
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
            AutoDual q2_sum = 0.0;
            for (int k = 0; k < ncomp; ++k) {
                q2_sum += x_dual[k] * cppargs.z[k] * cppargs.z[k];
            }
            AutoDual eps = dielectric_constant_rule_autodiff_cpp(cppargs.dielc_rule, x_dual, cppargs);
            AutoDual kappa = scalar_sqrt((rho * N_AV / 1.0e30) * E_CHRG * E_CHRG / kb / t / perm_vac * q2_sum / eps);
            AutoDual chi_sum = 0.0;
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
            AutoDual charge_radius_sum = 0.0;
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

}  // namespace

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

// EqID: dg_hs_dxk
double hs_contact_composition_derivative_cpp(
    double pair_diameter,
    double zeta2,
    double zeta3,
    double dzeta2_dx,
    double dzeta3_dx
) {
    return dzeta3_dx / std::pow(1.0 - zeta3, 2.0)
        + pair_diameter * (
            3.0 * dzeta2_dx / std::pow(1.0 - zeta3, 2.0)
            + 6.0 * zeta2 * dzeta3_dx / std::pow(1.0 - zeta3, 3.0)
        )
        + std::pow(pair_diameter, 2.0) * (
            4.0 * zeta2 * dzeta2_dx / std::pow(1.0 - zeta3, 3.0)
            + 6.0 * zeta2 * zeta2 * dzeta3_dx / std::pow(1.0 - zeta3, 4.0)
        );
}

namespace {

// EqID: dares_hs_dxk
vector<double> dadx_hs_cpp(const MixtureState &thermo, const HardChainState &hc_state, const add_args &cppargs) {
    int ncomp = static_cast<int>(thermo.d.size());
    vector<double> result(ncomp, 0.0);
    const auto &zeta = hc_state.zeta;

    auto hs_base_value = [&]() {
        return 1.0 / zeta[0] * (
            3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3])
            + std::pow(zeta[2], 3.0) / (zeta[3] * std::pow(1.0 - zeta[3], 2.0))
            + (std::pow(zeta[2], 3.0) / std::pow(zeta[3], 2.0) - zeta[0]) * std::log(1.0 - zeta[3])
        );
    }();

    vector<double> dzeta_dx(4, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        for (int l = 0; l < 4; ++l) {
            dzeta_dx[l] = PI / 6.0 * thermo.den * cppargs.m[i] * std::pow(thermo.d[i], l);
        }
        result[i] = -dzeta_dx[0] / zeta[0] * hs_base_value
            + 1.0 / zeta[0] * (
                3.0 * (dzeta_dx[1] * zeta[2] + zeta[1] * dzeta_dx[2]) / (1.0 - zeta[3])
                + 3.0 * zeta[1] * zeta[2] * dzeta_dx[3] / std::pow(1.0 - zeta[3], 2.0)
                + 3.0 * zeta[2] * zeta[2] * dzeta_dx[2] / zeta[3] / std::pow(1.0 - zeta[3], 2.0)
                + std::pow(zeta[2], 3.0) * dzeta_dx[3] * (3.0 * zeta[3] - 1.0) / zeta[3] / zeta[3] / std::pow(1.0 - zeta[3], 3.0)
                + std::log(1.0 - zeta[3]) * (
                    (3.0 * zeta[2] * zeta[2] * dzeta_dx[2] * zeta[3] - 2.0 * std::pow(zeta[2], 3.0) * dzeta_dx[3]) / std::pow(zeta[3], 3.0)
                    - dzeta_dx[0]
                )
                + (zeta[0] - std::pow(zeta[2], 3.0) / zeta[3] / zeta[3]) * dzeta_dx[3] / (1.0 - zeta[3])
            );
    }
    return result;
}

// EqID: dares_hc_dxk
vector<double> hc_contact_composition_terms_cpp(const MixtureState &thermo, const HardChainState &hc_state, const add_args &cppargs) {
    int ncomp = static_cast<int>(thermo.d.size());
    vector<double> terms(ncomp * ncomp, 0.0);
    vector<double> dzeta_dx(4, 0.0);
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int l = 0; l < 4; ++l) {
            dzeta_dx[l] = PI / 6.0 * thermo.den * cppargs.m[i] * std::pow(thermo.d[i], l);
        }
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            double pair_diameter = thermo.d[j] * thermo.d[j] / (thermo.d[j] + thermo.d[j]);
            terms[idx] = hs_contact_composition_derivative_cpp(
                pair_diameter,
                hc_state.zeta[2],
                hc_state.zeta[3],
                dzeta_dx[2],
                dzeta_dx[3]
            );
        }
    }
    return terms;
}

ContributionDadxResult dadx_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);

    const auto &zeta = hc_state.zeta;
    double log_sum = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        log_sum += x[i] * (cppargs.m[i] - 1.0) * std::log(hc_state.ghs[i * ncomp + i]);
    }
    double ares_hs = 1.0 / zeta[0] * (
        3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3])
        + std::pow(zeta[2], 3.0) / (zeta[3] * std::pow(1.0 - zeta[3], 2.0))
        + (std::pow(zeta[2], 3.0) / std::pow(zeta[3], 2.0) - zeta[0]) * std::log(1.0 - zeta[3])
    );
    result.ares = thermo.m_avg * ares_hs - log_sum;

    vector<double> dahs_dx = dadx_hs_cpp(thermo, hc_state, cppargs);
    vector<double> contact_composition_terms = hc_contact_composition_terms_cpp(thermo, hc_state, cppargs);
    for (int i = 0; i < ncomp; ++i) {
        double correction = 0.0;
        for (int j = 0; j < ncomp; ++j) {
            correction += x[j] * (cppargs.m[j] - 1.0) / hc_state.ghs[j * ncomp + j] * contact_composition_terms[i * ncomp + j];
        }
        result.dadx[i] = cppargs.m[i] * ares_hs
            + thermo.m_avg * dahs_dx[i]
            - correction
            - (cppargs.m[i] - 1.0) * std::log(hc_state.ghs[i * ncomp + i]);
    }

    if (cppargs.hc_dadx_diff_mode == 1) {
        result.dadx = contribution_dadx_fd_cpp(AresContributionKind::HC, t, rho, x, cppargs, result.ares);
    } else if (cppargs.hc_dadx_diff_mode == 2) {
        result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::HC, t, rho, x, cppargs);
    }

    double z_correction = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        double pair_diameter = pair_diameter_cpp(thermo.d[i], thermo.d[i]);
        z_correction += x[i] * (cppargs.m[i] - 1.0) / hc_state.ghs[i * ncomp + i]
            * hs_contact_density_derivative_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3]);
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    double dadrho_hs = hc_state.zeta[3] / (1.0 - hc_state.zeta[3])
        + 3.0 * hc_state.zeta[1] * hc_state.zeta[2] / hc_state.zeta[0] / std::pow(1.0 - hc_state.zeta[3], 2.0)
        + (3.0 * std::pow(hc_state.zeta[2], 3.0) - hc_state.zeta[3] * std::pow(hc_state.zeta[2], 3.0)) / hc_state.zeta[0] / std::pow(1.0 - hc_state.zeta[3], 3.0);
    result.z = thermo.m_avg * dadrho_hs - z_correction;
    return result;
}

// EqID: dares_disp_dxk
ContributionDadxResult dadx_disp_cpp(const MixtureState &thermo, const HardChainState &hc_state, const DispersionPolynomialState &dispersion, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);

    result.ares = -2.0 * PI * thermo.den * dispersion.I1 * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * dispersion.C1 * dispersion.I2 * thermo.m2e2s3;

    for (int i = 0; i < ncomp; ++i) {
        double dzeta3_dx = PI / 6.0 * thermo.den * cppargs.m[i] * std::pow(thermo.d[i], 3.0);
        double dI1_dx = 0.0;
        double dI2_dx = 0.0;
        double dm2es3_dx = 0.0;
        double dm2e2s3_dx = 0.0;
        for (int l = 0; l < 7; ++l) {
            double daa_dx = cppargs.m[i] / thermo.m_avg / thermo.m_avg * kDispersionA1[l]
                + cppargs.m[i] / thermo.m_avg / thermo.m_avg * (3.0 - 4.0 / thermo.m_avg) * kDispersionA2[l];
            double db_dx = cppargs.m[i] / thermo.m_avg / thermo.m_avg * kDispersionB1[l]
                + cppargs.m[i] / thermo.m_avg / thermo.m_avg * (3.0 - 4.0 / thermo.m_avg) * kDispersionB2[l];
            dI1_dx += dispersion.a[l] * l * dzeta3_dx * std::pow(hc_state.eta, l - 1) + daa_dx * std::pow(hc_state.eta, l);
            dI2_dx += dispersion.b[l] * l * dzeta3_dx * std::pow(hc_state.eta, l - 1) + db_dx * std::pow(hc_state.eta, l);
        }
        for (int j = 0; j < ncomp; ++j) {
            dm2es3_dx += x[j] * cppargs.m[j] * (thermo.e_ij[i * ncomp + j] / t) * std::pow(thermo.s_ij[i * ncomp + j], 3);
            dm2e2s3_dx += x[j] * cppargs.m[j] * std::pow(thermo.e_ij[i * ncomp + j] / t, 2) * std::pow(thermo.s_ij[i * ncomp + j], 3);
        }
        dm2es3_dx *= 2.0 * cppargs.m[i];
        dm2e2s3_dx *= 2.0 * cppargs.m[i];
        double dC1_dx = dispersion.C2 * dzeta3_dx - dispersion.C1 * dispersion.C1 * (
            cppargs.m[i] * (8.0 * hc_state.eta - 2.0 * hc_state.eta * hc_state.eta) / std::pow(1.0 - hc_state.eta, 4)
            - cppargs.m[i] * (20.0 * hc_state.eta - 27.0 * hc_state.eta * hc_state.eta + 12.0 * std::pow(hc_state.eta, 3) - 2.0 * std::pow(hc_state.eta, 4))
                / std::pow((1.0 - hc_state.eta) * (2.0 - hc_state.eta), 2)
        );
        result.dadx[i] = -2.0 * PI * thermo.den * (dI1_dx * thermo.m2es3 + dispersion.I1 * dm2es3_dx)
            - PI * thermo.den * ((cppargs.m[i] * dispersion.C1 * dispersion.I2 + thermo.m_avg * dC1_dx * dispersion.I2 + thermo.m_avg * dispersion.C1 * dI2_dx) * thermo.m2e2s3
            + thermo.m_avg * dispersion.C1 * dispersion.I2 * dm2e2s3_dx);
    }

    if (cppargs.disp_dadx_diff_mode == 1) {
        result.dadx = contribution_dadx_fd_cpp(AresContributionKind::DISP, t, rho, x, cppargs, result.ares);
    } else if (cppargs.disp_dadx_diff_mode == 2) {
        result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::DISP, t, rho, x, cppargs);
    }

    result.z = -2.0 * PI * thermo.den * dispersion.dEtaI1_deta * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * (dispersion.C1 * dispersion.dEtaI2_deta + dispersion.C2 * hc_state.eta * dispersion.I2) * thermo.m2e2s3;
    for (int i = 0; i < ncomp; ++i) {
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    return result;
}

ContributionDadxResult dadx_polar_cpp(const HardChainState &hc_state, const PolarIntermediateState &polar_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);
    if (!polar_state.active) {
        return result;
    }

    if (polar_state.second_order != 0.0) {
        result.ares = polar_state.second_order / (1.0 - polar_state.third_order / polar_state.second_order);
        for (int i = 0; i < ncomp; ++i) {
            result.dadx[i] = (polar_state.second_order_composition_terms[i] * (1.0 - polar_state.third_order / polar_state.second_order)
                + (polar_state.third_order_composition_terms[i] * polar_state.second_order - polar_state.third_order * polar_state.second_order_composition_terms[i]) / polar_state.second_order)
                / std::pow(1.0 - polar_state.third_order / polar_state.second_order, 2.0);
        }
        if (cppargs.polar_dadx_diff_mode == 1) {
            result.dadx = contribution_dadx_fd_cpp(AresContributionKind::POLAR, t, rho, x, cppargs, result.ares);
        } else if (cppargs.polar_dadx_diff_mode == 2) {
            result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::POLAR, t, rho, x, cppargs);
        }
        result.z = hc_state.eta * (
            (polar_state.second_order_density_term * (1.0 - polar_state.third_order / polar_state.second_order)
                + (polar_state.third_order_density_term * polar_state.second_order - polar_state.third_order * polar_state.second_order_density_term) / polar_state.second_order)
            / std::pow(1.0 - polar_state.third_order / polar_state.second_order, 2.0)
        );
        for (int i = 0; i < ncomp; ++i) {
            result.sum_x_dadx += x[i] * result.dadx[i];
        }
    }
    return result;
}

vector<double> association_site_fraction_composition_terms_cpp(
    const vector<double> &delta_ij,
    double den,
    const vector<double> &XA,
    const vector<double> &ddelta_dx,
    const vector<int> &site_component_index,
    const vector<double> &x_assoc,
    int ncomp
) {
    int num_sites = static_cast<int>(XA.size());
    vector<double> dXA_dx(ncomp * num_sites, 0.0);

    for (int k = 0; k < ncomp; ++k) {
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_sites, 1);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites, num_sites);

        int ij = 0;
        for (int i = 0; i < num_sites; ++i) {
            double direct_sum = 0.0;
            double delta_sum = 0.0;
            for (int j = 0; j < num_sites; ++j) {
                if (site_component_index[j] == k) {
                    direct_sum += XA[j] * delta_ij[ij];
                }
                delta_sum += x_assoc[j] * XA[j] * ddelta_dx[k * num_sites * num_sites + ij];
                A(i, j) = x_assoc[j] * delta_ij[ij];
                ++ij;
            }
            B(i) = -(direct_sum + delta_sum);
            A(i, i) += 1.0 / (den * XA[i] * XA[i]);
        }

        Eigen::MatrixXd solution = A.lu().solve(B);
        for (int i = 0; i < num_sites; ++i) {
            dXA_dx[k * num_sites + i] = solution(i);
        }
    }

    return dXA_dx;
}

// EqID: dares_assoc_dxk
ContributionDadxResult dadx_assoc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const AssociationIntermediateState &assoc_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);
    if (!assoc_state.active) {
        return result;
    }

    int num_sites = static_cast<int>(assoc_state.setup.site_component_index.size());
    vector<double> ddelta_dx(num_sites * num_sites * ncomp, 0.0);
    int idx_ddelta = 0;
    for (int k = 0; k < ncomp; ++k) {
        for (int i = 0; i < num_sites; ++i) {
            int comp_i = assoc_state.setup.site_component_index[i];
            for (int j = 0; j < num_sites; ++j) {
                int comp_j = assoc_state.setup.site_component_index[j];
                if (cppargs.assoc_matrix[i * num_sites + j] != 0) {
                    double pair_diameter = pair_diameter_cpp(thermo.d[comp_i], thermo.d[comp_j]);
                    double dzeta2_dx = PI / 6.0 * thermo.den * cppargs.m[k] * std::pow(thermo.d[k], 2);
                    double dzeta3_dx = PI / 6.0 * thermo.den * cppargs.m[k] * std::pow(thermo.d[k], 3);
                    double dghsd_dx = hs_contact_composition_derivative_cpp(
                        pair_diameter,
                        hc_state.zeta[2],
                        hc_state.zeta[3],
                        dzeta2_dx,
                        dzeta3_dx
                    );
                    double eABij = 0.5 * (cppargs.e_assoc[comp_i] + cppargs.e_assoc[comp_j]);
                    double volABij = association_volume_cpp(comp_i, comp_j, ncomp, thermo.s_ij, cppargs);
                    ddelta_dx[idx_ddelta] = dghsd_dx * (std::exp(eABij / t) - 1.0)
                        * std::pow(thermo.s_ij[comp_i * ncomp + comp_j], 3.0) * volABij;
                }
                ++idx_ddelta;
            }
        }
    }

    vector<double> dXA_dx = association_site_fraction_composition_terms_cpp(
        assoc_state.setup.delta_ij,
        thermo.den,
        assoc_state.XA,
        ddelta_dx,
        assoc_state.setup.site_component_index,
        assoc_state.setup.x_assoc,
        ncomp
    );

    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < num_sites; ++j) {
            result.dadx[i] += x[assoc_state.setup.site_component_index[j]] * dXA_dx[i * num_sites + j] * (1.0 / assoc_state.XA[j] - 0.5);
        }
    }

    for (int i = 0; i < num_sites; ++i) {
        int component_index = assoc_state.setup.site_component_index[i];
        result.dadx[component_index] += std::log(assoc_state.XA[i]) - 0.5 * assoc_state.XA[i] + 0.5;
        result.ares += x[component_index] * (std::log(assoc_state.XA[i]) - 0.5 * assoc_state.XA[i] + 0.5);
    }

    if (cppargs.assoc_dadx_diff_mode == 1) {
        result.dadx = contribution_dadx_fd_cpp(AresContributionKind::ASSOC, t, rho, x, cppargs, result.ares);
    } else if (cppargs.assoc_dadx_diff_mode == 2) {
        result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::ASSOC, t, rho, x, cppargs);
    }

    for (int i = 0; i < ncomp; ++i) {
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    return result;
}

// EqID: dares_dh_dxi
ContributionDadxResult dadx_ion_cpp(const MixtureState &thermo, const IonIntermediateState &ion_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);
    if (cppargs.z.empty()) {
        return result;
    }

    int dh_model = cppargs.DH_model;
    if (dh_model == 2) {
        throw ValueError("DH_model=2 (Bjerrum treatment) is reserved and not implemented.");
    }
    if ((dh_model != 0) && (dh_model != 1)) {
        throw ValueError("Unknown DH_model. Supported values are 0, 1, and reserved 2.");
    }
    if (!ion_state.active) {
        return result;
    }

    double K0 = E_CHRG * E_CHRG / (12.0 * PI * kb * t * perm_vac);
    result.ares = -K0 * ion_state.kappa / ion_state.dielectric.eps * ion_state.chi_sum;
    result.z = -ion_state.kappa / 24.0 / PI / kb / t / (ion_state.dielectric.eps * perm_vac) * ion_state.sigma_sum * E_CHRG * E_CHRG;

    if (cppargs.mu_DH_diff_mode == 1) {
        result.dadx = contribution_dadx_fd_cpp(AresContributionKind::ION, t, rho, x, cppargs, result.ares);
    } else if (cppargs.mu_DH_diff_mode == 2) {
        result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::ION, t, rho, x, cppargs);
    } else {
        const bool use_dh_deps = (cppargs.mu_DH_comp_dep_rel_perm != 0);
        for (int i = 0; i < ncomp; ++i) {
            double d_inv_eps_dx = use_dh_deps ? -ion_state.dielectric.deps_dx[i] / (ion_state.dielectric.eps * ion_state.dielectric.eps) : 0.0;
            double term1 = (ion_state.dkappa_dx[i] / ion_state.dielectric.eps + ion_state.kappa * d_inv_eps_dx) * ion_state.chi_sum;
            double term2 = ion_state.kappa / ion_state.dielectric.eps * ion_state.dchi_sum_dx[i];
            result.dadx[i] = -K0 * (term1 + term2);
        }
    }

    for (int i = 0; i < ncomp; ++i) {
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    return result;
}

ContributionDadxResult dadx_born_cpp(const BornIntermediateState &born_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);
    if (cppargs.z.empty()) {
        return result;
    }

    if (born_state.model == 1) {
        double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        result.ares = -Kborn * (1.0 - 1.0 / born_state.eps_value) * born_state.charge_radius_sum;

        if (cppargs.born_diff_mode == 1) {
            result.dadx = contribution_dadx_fd_cpp(AresContributionKind::BORN, t, rho, x, cppargs, result.ares);
        } else if (cppargs.born_diff_mode == 4) {
            result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::BORN, t, rho, x, cppargs);
        } else {
            for (int i = 0; i < ncomp; ++i) {
                double ion_part = 0.0;
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                    ion_part = (1.0 - 1.0 / born_state.eps_value) * cppargs.z[i] * cppargs.z[i] / d_born_i;
                }
                double eps_part = 0.0;
                if (cppargs.born_diff_mode == 2) {
                    eps_part = born_state.deps_dx[i] / (born_state.eps_value * born_state.eps_value);
                } else if (cppargs.born_diff_mode == 3) {
                    eps_part = 0.0;
                } else {
                    eps_part = born_state.charge_radius_sum * born_state.deps_dx[i] / (born_state.eps_value * born_state.eps_value);
                }
                result.dadx[i] = -Kborn * (ion_part + eps_part);
            }
        }
    } else if (born_state.model == 2) {
        const double eps_r_ion = 8.0;
        const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        result.ares = -Kborn * born_state.shell.sum_bracket;

        if (cppargs.born_diff_mode == 1) {
            result.dadx = contribution_dadx_fd_cpp(AresContributionKind::BORN, t, rho, x, cppargs, result.ares);
        } else if (cppargs.born_diff_mode == 4) {
            result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::BORN, t, rho, x, cppargs);
        } else {
            const double inv_eps2 = 1.0 / (born_state.eps_value * born_state.eps_value);
            const double shell_coeff = 1.0 / eps_r_ion - 1.0 / born_state.eps_value;
            const bool use_deps = (cppargs.mu_born_comp_dep_rel_perm != 0);
            const bool use_shell_chain = (cppargs.mu_born_comp_dep_delta_d != 0);
            const double deps_multiplier = (cppargs.mu_born_include_sum_term != 0) ? born_state.shell.sum_gap : 1.0;
            for (int k = 0; k < ncomp; ++k) {
                double direct_part = 0.0;
                if (std::abs(cppargs.z[k]) > 1e-12) {
                    direct_part = cppargs.z[k] * cppargs.z[k] * born_state.shell.bracket[k];
                }
                double deps_part = use_deps ? deps_multiplier * born_state.deps_dx[k] * inv_eps2 : 0.0;
                double ddelta_part = use_shell_chain ? shell_coeff * born_state.shell.sum_dpref_over_D2 * born_state.shell.f_k[k] : 0.0;
                result.dadx[k] = -Kborn * (direct_part + deps_part + ddelta_part);
            }
        }
    } else if (born_state.model != 0) {
        throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
    }

    for (int i = 0; i < ncomp; ++i) {
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    return result;
}

}  // namespace

CompositionContributionResult composition_derivative_residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    MixtureState thermo = mixture_state_cpp(t, rho, x, cppargs, false);
    HardChainState hc_state = hard_chain_state_cpp(thermo, x, cppargs);
    DadrhoResult dadrho_result = dadrho_result_cpp(t, rho, x, cppargs);
    const ScalarContributionTerms &z_raw_terms = dadrho_result.raw;
    ScalarContributionTerms z_terms = compressibility_terms_from_dadrho_cpp(dadrho_result);

    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(thermo.m_avg, hc_state.eta);
    PolarIntermediateState polar_state = polar_intermediate_state_cpp(thermo, hc_state, 0.0, t, x, cppargs, true, false, true);
    AssociationIntermediateState assoc_state = association_intermediate_state_cpp(thermo, hc_state, t, x, cppargs, false, false);
    IonIntermediateState ion_state = ion_intermediate_state_cpp(thermo, t, x, cppargs, true);
    BornIntermediateState born_state = born_intermediate_state_cpp(t, x, cppargs, false, true);

    ContributionDadxResult hc = dadx_hc_cpp(thermo, hc_state, t, rho, x, cppargs);
    ContributionDadxResult disp = dadx_disp_cpp(thermo, hc_state, dispersion, t, rho, x, cppargs);
    ContributionDadxResult polar = dadx_polar_cpp(hc_state, polar_state, t, rho, x, cppargs);
    ContributionDadxResult assoc = dadx_assoc_cpp(thermo, hc_state, assoc_state, t, rho, x, cppargs);
    ContributionDadxResult ion = dadx_ion_cpp(thermo, ion_state, t, rho, x, cppargs);
    ContributionDadxResult born = dadx_born_cpp(born_state, t, rho, x, cppargs);

    CompositionContributionResult result;
    result.dadx = make_vector_terms(hc.dadx, disp.dadx, polar.dadx, assoc.dadx, ion.dadx, born.dadx, vector<double>());
    result.ares = make_scalar_terms(hc.ares, disp.ares, polar.ares, assoc.ares, ion.ares, born.ares,
        hc.ares + disp.ares + polar.ares + assoc.ares + ion.ares + born.ares);
    result.sum_x_dadx = make_scalar_terms(hc.sum_x_dadx, disp.sum_x_dadx, polar.sum_x_dadx, assoc.sum_x_dadx,
        ion.sum_x_dadx, born.sum_x_dadx,
        hc.sum_x_dadx + disp.sum_x_dadx + polar.sum_x_dadx + assoc.sum_x_dadx + ion.sum_x_dadx + born.sum_x_dadx);
    result.z_raw = z_raw_terms;
    result.z = z_terms;

    vector<double> total(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        total[i] = hc.dadx[i] + disp.dadx[i] + polar.dadx[i] + assoc.dadx[i] + ion.dadx[i] + born.dadx[i];
    }
    result.dadx.total = total;
    return result;
}
