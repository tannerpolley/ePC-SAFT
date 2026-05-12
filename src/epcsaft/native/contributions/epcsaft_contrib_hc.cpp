#include "epcsaft_contrib_internal.h"

using thermo_detail::AutodiffHardChainState;
using thermo_detail::ContributionDadxResult;
using thermo_detail::HardChainState;
using thermo_detail::MixtureState;
using thermo_detail::parameter_setup_detail::pair_diameter_cpp;

namespace hard_chain_detail {

template <typename Scalar>
static Scalar hs_contact_value_scalar_cpp(double pair_diameter, const Scalar &zeta2, const Scalar &zeta3) {
    const Scalar one = scalar_constant<Scalar>(1.0);
    return one / (one - zeta3)
        + pair_diameter * 3.0 * zeta2 / scalar_pow(one - zeta3, 2)
        + scalar_pow(pair_diameter, 2.0) * 2.0 * zeta2 * zeta2 / scalar_pow(one - zeta3, 3);
}

template <typename Scalar>
struct HardChainStateScalar {
    vector<Scalar> zeta;
    vector<Scalar> ghs;
    Scalar eta = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
// EqID: zeta_n
// EqID: zeta_n_xk
// EqID: zeta3_eta
static HardChainStateScalar<Scalar> hard_chain_state_scalar_cpp(double den, const vector<double> &d, const vector<Scalar> &x, const add_args &cppargs) {
    HardChainStateScalar<Scalar> state;
    int ncomp = static_cast<int>(x.size());
    state.zeta.assign(4, scalar_constant<Scalar>(0.0));
    state.ghs.assign(ncomp * ncomp, scalar_constant<Scalar>(0.0));

    for (int k = 0; k < 4; ++k) {
        Scalar summ = scalar_constant<Scalar>(0.0);
        for (int j = 0; j < ncomp; ++j) {
            summ += x[j] * cppargs.m[j] * scalar_pow(d[j], k);
        }
        state.zeta[k] = PI / 6.0 * den * summ;
    }

    state.eta = state.zeta[3];

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            double pair_diameter = pair_diameter_cpp(d[i], d[j]);
            state.ghs[idx] = hs_contact_value_scalar_cpp(pair_diameter, state.zeta[2], state.zeta[3]);
        }
    }

    return state;
}

}  // namespace hard_chain_detail

// EqID: g_hs_contact
double hs_contact_value_cpp(double pair_diameter, double zeta2, double zeta3) {
    return hard_chain_detail::hs_contact_value_scalar_cpp(pair_diameter, zeta2, zeta3);
}

HardChainState hard_chain_state_cpp(const MixtureState &thermo, const vector<double> &x, const add_args &cppargs) {
    auto scalar_state = hard_chain_detail::hard_chain_state_scalar_cpp(thermo.den, thermo.d, x, cppargs);
    HardChainState state;
    state.zeta = std::move(scalar_state.zeta);
    state.ghs = std::move(scalar_state.ghs);
    state.eta = scalar_state.eta;
    return state;
}

AutodiffHardChainState hard_chain_state_autodiff_cpp(double den, const vector<double> &d, const vector<AutoDual> &x, const add_args &cppargs) {
    auto scalar_state = hard_chain_detail::hard_chain_state_scalar_cpp(den, d, x, cppargs);
    AutodiffHardChainState state;
    state.zeta = std::move(scalar_state.zeta);
    state.ghs = std::move(scalar_state.ghs);
    state.eta = scalar_state.eta;
    return state;
}

namespace hard_chain_detail {

// EqID: hs_ares_dadrho
static double dadrho_hs_cpp(const HardChainState &hc_state) {
    const auto &zeta = hc_state.zeta;
    return zeta[3] / (1.0 - zeta[3])
        + 3.0 * zeta[1] * zeta[2] / zeta[0] / std::pow(1.0 - zeta[3], 2.0)
        + (3.0 * std::pow(zeta[2], 3.0) - zeta[3] * std::pow(zeta[2], 3.0)) / zeta[0] / std::pow(1.0 - zeta[3], 3.0);
}

// EqID: hs_ares_dT
static double dadt_hs_cpp(const HardChainState &hc_state, const vector<double> &dzeta_dt) {
    const auto &zeta = hc_state.zeta;
    return 1.0 / zeta[0] * (
        3.0 * (dzeta_dt[1] * zeta[2] + zeta[1] * dzeta_dt[2]) / (1.0 - zeta[3])
        + 3.0 * zeta[1] * zeta[2] * dzeta_dt[3] / std::pow(1.0 - zeta[3], 2.0)
        + 3.0 * std::pow(zeta[2], 2.0) * dzeta_dt[2] / zeta[3] / std::pow(1.0 - zeta[3], 2.0)
        + std::pow(zeta[2], 3.0) * dzeta_dt[3] * (3.0 * zeta[3] - 1.0) / std::pow(zeta[3], 2.0) / std::pow(1.0 - zeta[3], 3.0)
        + (3.0 * std::pow(zeta[2], 2.0) * dzeta_dt[2] * zeta[3] - 2.0 * std::pow(zeta[2], 3.0) * dzeta_dt[3]) / std::pow(zeta[3], 3.0) * std::log(1.0 - zeta[3])
        + (zeta[0] - std::pow(zeta[2], 3.0) / std::pow(zeta[3], 2.0)) * dzeta_dt[3] / (1.0 - zeta[3])
    );
}

// EqID: hs_ares_dxk
static vector<double> dadx_hs_cpp(const MixtureState &thermo, const HardChainState &hc_state, const add_args &cppargs) {
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

static vector<double> hc_contact_composition_terms_cpp(const MixtureState &thermo, const HardChainState &hc_state, const add_args &cppargs) {
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
            double pair_diameter = pair_diameter_cpp(thermo.d[j], thermo.d[j]);
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

}  // namespace hard_chain_detail

// EqID: ghs_contact_dadrho
double hs_contact_density_derivative_cpp(double pair_diameter, double zeta2, double zeta3) {
    return zeta3 / std::pow(1.0 - zeta3, 2.0)
        + pair_diameter * (3.0 * zeta2 / std::pow(1.0 - zeta3, 2.0) + 6.0 * zeta2 * zeta3 / std::pow(1.0 - zeta3, 3.0))
        + std::pow(pair_diameter, 2.0) * (4.0 * zeta2 * zeta2 / std::pow(1.0 - zeta3, 3.0) + 6.0 * zeta2 * zeta2 * zeta3 / std::pow(1.0 - zeta3, 4.0));
}

// EqID: ghs_contact_dT
double hs_contact_time_derivative_cpp(
    double pair_diameter,
    double pair_diameter_dt,
    double zeta2,
    double zeta3,
    double dzeta2_dt,
    double dzeta3_dt
) {
    return dzeta3_dt / std::pow(1.0 - zeta3, 2.0)
        + 3.0 * (pair_diameter_dt * zeta2 + pair_diameter * dzeta2_dt) / std::pow(1.0 - zeta3, 2.0)
        + 4.0 * pair_diameter * zeta2 * (1.5 * dzeta3_dt + pair_diameter_dt * zeta2 + pair_diameter * dzeta2_dt) / std::pow(1.0 - zeta3, 3.0)
        + 6.0 * std::pow(pair_diameter * zeta2, 2.0) * dzeta3_dt / std::pow(1.0 - zeta3, 4.0);
}

// EqID: ghs_contact_dxk
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

// EqID: zeta_n_dT
vector<double> dzeta_dt_cpp(const MixtureState &thermo, const vector<double> &x, const add_args &cppargs) {
    vector<double> dzeta_dt(4, 0.0);
    for (int k = 1; k < 4; ++k) {
        double summ = 0.0;
        for (int j = 0; j < static_cast<int>(x.size()); ++j) {
            summ += x[j] * cppargs.m[j] * k * thermo.dd_dt[j] * std::pow(thermo.d[j], k - 1);
        }
        dzeta_dt[k] = PI / 6.0 * thermo.den * summ;
    }
    return dzeta_dt;
}

vector<double> hc_contact_time_terms_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &dzeta_dt) {
    int ncomp = static_cast<int>(thermo.d.size());
    vector<double> terms(ncomp * ncomp, 0.0);
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            double pair_diameter = pair_diameter_cpp(thermo.d[i], thermo.d[j]);
            double pair_diameter_dt = pair_diameter * (
                thermo.dd_dt[i] / thermo.d[i] + thermo.dd_dt[j] / thermo.d[j]
                - (thermo.dd_dt[i] + thermo.dd_dt[j]) / (thermo.d[i] + thermo.d[j])
            );
            terms[idx] = hs_contact_time_derivative_cpp(
                pair_diameter,
                pair_diameter_dt,
                hc_state.zeta[2],
                hc_state.zeta[3],
                dzeta_dt[2],
                dzeta_dt[3]
            );
        }
    }
    return terms;
}

// EqID: hc_ares_dadrho
double dadrho_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    double summ = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        double pair_diameter = pair_diameter_cpp(thermo.d[i], thermo.d[i]);
        double dghs_drho = hs_contact_density_derivative_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3]);
        summ += x[i] * (cppargs.m[i] - 1.0) / hc_state.ghs[i * ncomp + i] * dghs_drho;
    }
    return thermo.m_avg * hard_chain_detail::dadrho_hs_cpp(hc_state) - summ;
}

// EqID: hc_ares_dT
double dadt_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &dzeta_dt, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    vector<double> contact_time_terms = hc_contact_time_terms_cpp(thermo, hc_state, dzeta_dt);
    double summ = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        summ += x[i] * (cppargs.m[i] - 1.0) * contact_time_terms[i * ncomp + i] / hc_state.ghs[i * ncomp + i];
    }
    return thermo.m_avg * hard_chain_detail::dadt_hs_cpp(hc_state, dzeta_dt) - summ;
}

// EqID: hc_ares_dxk
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

    vector<double> dahs_dx = hard_chain_detail::dadx_hs_cpp(thermo, hc_state, cppargs);
    vector<double> contact_composition_terms = hard_chain_detail::hc_contact_composition_terms_cpp(thermo, hc_state, cppargs);
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
        throw ValueError("backend_unavailable: hard-chain composition derivative backend is unavailable.");
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
    result.z = thermo.m_avg * hard_chain_detail::dadrho_hs_cpp(hc_state) - z_correction;
    return result;
}
