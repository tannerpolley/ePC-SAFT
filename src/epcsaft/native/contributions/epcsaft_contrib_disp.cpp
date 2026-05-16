#include "epcsaft_contrib_internal.h"

using namespace thermo_detail;

// EqID: c1_disp
// EqID: c2_disp
// EqID: i1_disp
// EqID: deta_i1_deta
// EqID: i2_disp
// EqID: deta_i2_deta
// EqID: a_i_mbar
// EqID: b_i_mbar
DispersionPolynomialState dispersion_polynomials_cpp(double m_avg, double eta) {
    DispersionPolynomialState state;
    double c1 = (m_avg - 1.0) / m_avg;
    double c2 = (m_avg - 2.0) / m_avg;
    for (size_t i = 0; i < state.a.size(); ++i) {
        state.a[i] = kDispersionA0[i] + c1 * kDispersionA1[i] + c1 * c2 * kDispersionA2[i];
        state.b[i] = kDispersionB0[i] + c1 * kDispersionB1[i] + c1 * c2 * kDispersionB2[i];
        state.I1 += state.a[i] * std::pow(eta, static_cast<int>(i));
        state.I2 += state.b[i] * std::pow(eta, static_cast<int>(i));
        if (i > 0) {
            state.dI1_deta += state.a[i] * static_cast<double>(i) * std::pow(eta, static_cast<int>(i - 1));
            state.dI2_deta += state.b[i] * static_cast<double>(i) * std::pow(eta, static_cast<int>(i - 1));
        }
        state.dEtaI1_deta += state.a[i] * static_cast<double>(i + 1) * std::pow(eta, static_cast<int>(i));
        state.dEtaI2_deta += state.b[i] * static_cast<double>(i + 1) * std::pow(eta, static_cast<int>(i));
    }
    state.C1 = 1.0 / (1.0
        + m_avg * (8.0 * eta - 2.0 * eta * eta) / std::pow(1.0 - eta, 4.0)
        + (1.0 - m_avg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * std::pow(eta, 3.0) - 2.0 * std::pow(eta, 4.0))
            / std::pow((1.0 - eta) * (2.0 - eta), 2.0));
    state.C2 = -state.C1 * state.C1 * (
        m_avg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / std::pow(1.0 - eta, 5.0)
        + (1.0 - m_avg) * (2.0 * std::pow(eta, 3.0) + 12.0 * eta * eta - 48.0 * eta + 40.0)
            / std::pow((1.0 - eta) * (2.0 - eta), 3.0));
    return state;
}

// EqID: disp_ares_dadrho
double dadrho_disp_cpp(const MixtureState &thermo, const HardChainState &hc_state, const DispersionPolynomialState &dispersion) {
    return -2.0 * PI * thermo.den * dispersion.dEtaI1_deta * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * (dispersion.C1 * dispersion.dEtaI2_deta + dispersion.C2 * hc_state.eta * dispersion.I2) * thermo.m2e2s3;
}

// EqID: disp_ares_dT
double dadt_disp_cpp(const MixtureState &thermo, double deta_dt, double t, const DispersionPolynomialState &dispersion) {
    double dI1_dt = dispersion.dI1_deta * deta_dt;
    double dI2_dt = dispersion.dI2_deta * deta_dt;
    double dC1_dt = dispersion.C2 * deta_dt;
    return -2.0 * PI * thermo.den * (dI1_dt - dispersion.I1 / t) * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * (dC1_dt * dispersion.I2 + dispersion.C1 * dI2_dt - 2.0 * dispersion.C1 * dispersion.I2 / t) * thermo.m2e2s3;
}

// EqID: disp_ares_dxk
// EqID: c1_xk
// EqID: m2epssigma3_xk
// EqID: m2eps2sigma3_xk
// EqID: i1_xk
// EqID: i2_xk
// EqID: a_i_xk
// EqID: b_i_xk
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
        throw ValueError("unsupported: dispersion composition derivative backend is not enabled.");
    } else if (cppargs.disp_dadx_diff_mode == 2) {
        result.dadx = contribution_dadx_cppad_cpp(AresContributionKind::DISP, t, rho, x, cppargs);
    }

    result.z = -2.0 * PI * thermo.den * dispersion.dEtaI1_deta * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * (dispersion.C1 * dispersion.dEtaI2_deta + dispersion.C2 * hc_state.eta * dispersion.I2) * thermo.m2e2s3;
    for (int i = 0; i < ncomp; ++i) {
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    return result;
}
