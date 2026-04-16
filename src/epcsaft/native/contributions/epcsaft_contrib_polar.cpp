#include "epcsaft_contrib_internal.h"

using namespace thermo_detail;

namespace {

vector<double> dipole_strengths_cpp(const vector<double> &x, const add_args &cppargs) {
    vector<double> dipm_sq(x.size(), 0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        dipm_sq[i] = std::pow(cppargs.dipm[i], 2.0)
            / (cppargs.m[i] * cppargs.e[i] * std::pow(cppargs.s[i], 3.0))
            * kDipoleConversion;
    }
    return dipm_sq;
}

}  // namespace

vector<double> dipole_coefficients_cpp(const std::array<double, 5> &c0, const std::array<double, 5> &c1, const std::array<double, 5> &c2, double m) {
    vector<double> coeffs(5, 0.0);
    double a = (m - 1.0) / m;
    double b = (m - 2.0) / m;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        coeffs[i] = c0[i] + a * c1[i] + a * b * c2[i];
    }
    return coeffs;
}

PolarIntermediateState polar_intermediate_state_cpp(
    const MixtureState &thermo,
    const HardChainState &hc_state,
    double deta_dt,
    double t,
    const vector<double> &x,
    const add_args &cppargs,
    bool include_density,
    bool include_dt,
    bool include_dx
) {
    PolarIntermediateState state;
    if (cppargs.dipm.empty()) {
        return state;
    }
    state.active = true;

    int ncomp = static_cast<int>(x.size());
    vector<double> dipm_sq = dipole_strengths_cpp(x, cppargs);
    if (include_dx) {
        state.second_order_composition_terms.assign(ncomp, 0.0);
        state.third_order_composition_terms.assign(ncomp, 0.0);
    }

    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            double m_ij = std::sqrt(cppargs.m[i] * cppargs.m[j]);
            if (m_ij > 2.0) {
                m_ij = 2.0;
            }
            vector<double> adip = dipole_coefficients_cpp(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
            vector<double> bdip = dipole_coefficients_cpp(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
            double J2 = 0.0;
            double dJ2_det = 0.0;
            double detJ2_det = 0.0;
            double dJ2_dt = 0.0;
            for (int l = 0; l < 5; ++l) {
                double factor = adip[l] + bdip[l] * thermo.e_ij[i * ncomp + j] / t;
                J2 += factor * std::pow(hc_state.eta, l);
                if (include_dx) {
                    dJ2_det += factor * l * std::pow(hc_state.eta, l - 1);
                }
                if (include_density) {
                    detJ2_det += factor * (l + 1) * std::pow(hc_state.eta, l);
                }
                if (include_dt) {
                    dJ2_dt += adip[l] * l * std::pow(hc_state.eta, l - 1) * deta_dt
                        + bdip[l] * thermo.e_ij[j * ncomp + j] * (
                            1.0 / t * l * std::pow(hc_state.eta, l - 1) * deta_dt
                            - 1.0 / std::pow(t, 2.0) * std::pow(hc_state.eta, l)
                        );
                }
            }

            double pair2 = x[i] * x[j] * thermo.e_ij[i * ncomp + i] / t * thermo.e_ij[j * ncomp + j] / t
                * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0)
                / std::pow(thermo.s_ij[i * ncomp + j], 3.0)
                * cppargs.dip_num[i] * cppargs.dip_num[j] * dipm_sq[i] * dipm_sq[j];
            state.second_order += pair2 * J2;
            if (include_density) {
                state.second_order_density_term += pair2 * detJ2_det;
            }
            if (include_dt) {
                state.second_order_temperature_term += x[i] * x[j] * thermo.e_ij[i * ncomp + i] * thermo.e_ij[j * ncomp + j]
                    * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0)
                    / std::pow(thermo.s_ij[i * ncomp + j], 3.0) * cppargs.dip_num[i] * cppargs.dip_num[j] * dipm_sq[i] * dipm_sq[j]
                    * (dJ2_dt / std::pow(t, 2.0) - 2.0 * J2 / std::pow(t, 3.0));
            }
            if (include_dx) {
                if (i == j) {
                    state.second_order_composition_terms[i] += thermo.e_ij[i * ncomp + i] / t * thermo.e_ij[j * ncomp + j] / t
                        * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0)
                        / std::pow(thermo.s_ij[i * ncomp + j], 3.0) * cppargs.dip_num[i] * cppargs.dip_num[j] * dipm_sq[i] * dipm_sq[j]
                        * (x[i] * x[j] * dJ2_det * PI / 6.0 * thermo.den * cppargs.m[i] * std::pow(thermo.d[i], 3.0) + 2.0 * x[j] * J2);
                } else {
                    state.second_order_composition_terms[i] += thermo.e_ij[i * ncomp + i] / t * thermo.e_ij[j * ncomp + j] / t
                        * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0)
                        / std::pow(thermo.s_ij[i * ncomp + j], 3.0) * cppargs.dip_num[i] * cppargs.dip_num[j] * dipm_sq[i] * dipm_sq[j]
                        * (x[i] * x[j] * dJ2_det * PI / 6.0 * thermo.den * cppargs.m[i] * std::pow(thermo.d[i], 3.0) + x[j] * J2);
                }
            }

            for (int k = 0; k < ncomp; ++k) {
                double m_ijk = std::pow(cppargs.m[i] * cppargs.m[j] * cppargs.m[k], 1.0 / 3.0);
                if (m_ijk > 2.0) {
                    m_ijk = 2.0;
                }
                vector<double> cdip = dipole_coefficients_cpp(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
                double J3 = 0.0;
                double detJ3_det = 0.0;
                double dJ3_dt = 0.0;
                double dJ3_det = 0.0;
                for (int l = 0; l < 5; ++l) {
                    J3 += cdip[l] * std::pow(hc_state.eta, l);
                    if (include_density) {
                        detJ3_det += cdip[l] * (l + 2) * std::pow(hc_state.eta, l + 1);
                    }
                    if (include_dt) {
                        dJ3_dt += cdip[l] * l * std::pow(hc_state.eta, l - 1) * deta_dt;
                    }
                    if (include_dx) {
                        dJ3_det += cdip[l] * l * std::pow(hc_state.eta, l - 1);
                    }
                }

                double pair3 = x[i] * x[j] * x[k] * thermo.e_ij[i * ncomp + i] / t * thermo.e_ij[j * ncomp + j] / t * thermo.e_ij[k * ncomp + k] / t
                    * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0) * std::pow(thermo.s_ij[k * ncomp + k], 3.0)
                    / thermo.s_ij[i * ncomp + j] / thermo.s_ij[i * ncomp + k] / thermo.s_ij[j * ncomp + k]
                    * cppargs.dip_num[i] * cppargs.dip_num[j] * cppargs.dip_num[k] * dipm_sq[i] * dipm_sq[j] * dipm_sq[k];
                state.third_order += pair3 * J3;
                if (include_density) {
                    state.third_order_density_term += pair3 * detJ3_det;
                }
                if (include_dt) {
                    state.third_order_temperature_term += x[i] * x[j] * x[k] * thermo.e_ij[i * ncomp + i] * thermo.e_ij[j * ncomp + j] * thermo.e_ij[k * ncomp + k]
                        * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0) * std::pow(thermo.s_ij[k * ncomp + k], 3.0)
                        / thermo.s_ij[i * ncomp + j] / thermo.s_ij[i * ncomp + k] / thermo.s_ij[j * ncomp + k]
                        * cppargs.dip_num[i] * cppargs.dip_num[j] * cppargs.dip_num[k] * dipm_sq[i] * dipm_sq[j] * dipm_sq[k]
                        * (-3.0 * J3 / std::pow(t, 4.0) + dJ3_dt / std::pow(t, 3.0));
                }
                if (include_dx) {
                    if ((i == j) && (i == k)) {
                        state.third_order_composition_terms[i] += thermo.e_ij[i * ncomp + i] / t * thermo.e_ij[j * ncomp + j] / t * thermo.e_ij[k * ncomp + k] / t
                            * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0) * std::pow(thermo.s_ij[k * ncomp + k], 3.0)
                            / thermo.s_ij[i * ncomp + j] / thermo.s_ij[i * ncomp + k] / thermo.s_ij[j * ncomp + k]
                            * cppargs.dip_num[i] * cppargs.dip_num[j] * cppargs.dip_num[k] * dipm_sq[i] * dipm_sq[j] * dipm_sq[k]
                            * (x[i] * x[j] * x[k] * dJ3_det * PI / 6.0 * thermo.den * cppargs.m[i] * std::pow(thermo.d[i], 3.0) + 3.0 * x[j] * x[k] * J3);
                    } else if ((i == j) || (i == k)) {
                        state.third_order_composition_terms[i] += thermo.e_ij[i * ncomp + i] / t * thermo.e_ij[j * ncomp + j] / t * thermo.e_ij[k * ncomp + k] / t
                            * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0) * std::pow(thermo.s_ij[k * ncomp + k], 3.0)
                            / thermo.s_ij[i * ncomp + j] / thermo.s_ij[i * ncomp + k] / thermo.s_ij[j * ncomp + k]
                            * cppargs.dip_num[i] * cppargs.dip_num[j] * cppargs.dip_num[k] * dipm_sq[i] * dipm_sq[j] * dipm_sq[k]
                            * (x[i] * x[j] * x[k] * dJ3_det * PI / 6.0 * thermo.den * cppargs.m[i] * std::pow(thermo.d[i], 3.0) + 2.0 * x[j] * x[k] * J3);
                    } else {
                        state.third_order_composition_terms[i] += thermo.e_ij[i * ncomp + i] / t * thermo.e_ij[j * ncomp + j] / t * thermo.e_ij[k * ncomp + k] / t
                            * std::pow(thermo.s_ij[i * ncomp + i], 3.0) * std::pow(thermo.s_ij[j * ncomp + j], 3.0) * std::pow(thermo.s_ij[k * ncomp + k], 3.0)
                            / thermo.s_ij[i * ncomp + j] / thermo.s_ij[i * ncomp + k] / thermo.s_ij[j * ncomp + k]
                            * cppargs.dip_num[i] * cppargs.dip_num[j] * cppargs.dip_num[k] * dipm_sq[i] * dipm_sq[j] * dipm_sq[k]
                            * (x[i] * x[j] * x[k] * dJ3_det * PI / 6.0 * thermo.den * cppargs.m[i] * std::pow(thermo.d[i], 3.0) + x[j] * x[k] * J3);
                    }
                }
            }
        }
    }

    state.second_order = -PI * thermo.den * state.second_order;
    state.third_order = -4.0 / 3.0 * PI * PI * thermo.den * thermo.den * state.third_order;
    if (include_density) {
        state.second_order_density_term = -PI * thermo.den / hc_state.eta * state.second_order_density_term;
        state.third_order_density_term = -4.0 / 3.0 * PI * PI * thermo.den / hc_state.eta * thermo.den / hc_state.eta * state.third_order_density_term;
    }
    if (include_dt) {
        state.second_order_temperature_term = -PI * thermo.den * state.second_order_temperature_term;
        state.third_order_temperature_term = -4.0 / 3.0 * PI * PI * thermo.den * thermo.den * state.third_order_temperature_term;
    }
    if (include_dx) {
        for (int i = 0; i < ncomp; ++i) {
            state.second_order_composition_terms[i] = -PI * thermo.den * state.second_order_composition_terms[i];
            state.third_order_composition_terms[i] = -4.0 / 3.0 * PI * PI * thermo.den * thermo.den * state.third_order_composition_terms[i];
        }
    }

    return state;
}
