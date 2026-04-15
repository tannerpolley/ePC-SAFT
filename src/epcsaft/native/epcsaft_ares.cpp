#include "epcsaft_core_internal.h"

using namespace thermo_detail;

namespace {

// EqID: ares_hc
double ares_hc_cpp(const ThermoCommonState &thermo, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    const auto &zeta = thermo.zeta;
    const auto &ghs = thermo.ghs;
    double ares_hs = 1.0 / zeta[0] * (
        3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3])
        + std::pow(zeta[2], 3.0) / (zeta[3] * std::pow(1.0 - zeta[3], 2.0))
        + (std::pow(zeta[2], 3.0) / std::pow(zeta[3], 2.0) - zeta[0]) * std::log(1.0 - zeta[3])
    );

    double summ = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        summ += x[i] * (cppargs.m[i] - 1.0) * std::log(ghs[i * ncomp + i]);
    }
    return thermo.m_avg * ares_hs - summ;
}

// EqID: ares_disp
double ares_disp_cpp(const ThermoCommonState &thermo) {
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(thermo.m_avg, thermo.eta);
    return -2.0 * PI * thermo.den * dispersion.I1 * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * dispersion.C1 * dispersion.I2 * thermo.m2e2s3;
}

double ares_polar_cpp(const ThermoCommonState &thermo, double t, const vector<double> &x, const add_args &cppargs) {
    if (cppargs.dipm.empty()) {
        return 0.0;
    }

    int ncomp = static_cast<int>(x.size());
    const auto &e_ij = thermo.e_ij;
    const auto &s_ij = thermo.s_ij;
    double eta = thermo.eta;
    double den = thermo.den;
    double A2 = 0.0;
    double A3 = 0.0;
    vector<double> dipmSQ(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        dipmSQ[i] = std::pow(cppargs.dipm[i], 2.0)
            / (cppargs.m[i] * cppargs.e[i] * std::pow(cppargs.s[i], 3.0))
            * kDipoleConversion;
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
            for (int l = 0; l < 5; ++l) {
                J2 += (adip[l] + bdip[l] * e_ij[j * ncomp + j] / t) * std::pow(eta, l);
            }
            A2 += x[i] * x[j] * e_ij[i * ncomp + i] / t * e_ij[j * ncomp + j] / t
                * std::pow(s_ij[i * ncomp + i], 3.0) * std::pow(s_ij[j * ncomp + j], 3.0)
                / std::pow(s_ij[i * ncomp + j], 3.0)
                * cppargs.dip_num[i] * cppargs.dip_num[j] * dipmSQ[i] * dipmSQ[j] * J2;

            for (int k = 0; k < ncomp; ++k) {
                double m_ijk = std::pow(cppargs.m[i] * cppargs.m[j] * cppargs.m[k], 1.0 / 3.0);
                if (m_ijk > 2.0) {
                    m_ijk = 2.0;
                }
                vector<double> cdip = dipole_coefficients_cpp(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
                double J3 = 0.0;
                for (int l = 0; l < 5; ++l) {
                    J3 += cdip[l] * std::pow(eta, l);
                }
                A3 += x[i] * x[j] * x[k] * e_ij[i * ncomp + i] / t * e_ij[j * ncomp + j] / t * e_ij[k * ncomp + k] / t
                    * std::pow(s_ij[i * ncomp + i], 3.0) * std::pow(s_ij[j * ncomp + j], 3.0) * std::pow(s_ij[k * ncomp + k], 3.0)
                    / s_ij[i * ncomp + j] / s_ij[i * ncomp + k] / s_ij[j * ncomp + k]
                    * cppargs.dip_num[i] * cppargs.dip_num[j] * cppargs.dip_num[k]
                    * dipmSQ[i] * dipmSQ[j] * dipmSQ[k] * J3;
            }
        }
    }

    A2 = -PI * den * A2;
    A3 = -4.0 / 3.0 * PI * PI * den * den * A3;
    if (A2 == 0.0) {
        return 0.0;
    }
    return A2 / (1.0 - A3 / A2);
}

// EqID: ares_assoc
double ares_assoc_cpp(const ThermoCommonState &thermo, double t, const vector<double> &x, const add_args &cppargs) {
    if (cppargs.e_assoc.empty()) {
        return 0.0;
    }

    AssociationSetup association = association_setup_cpp(x, cppargs, thermo.s_ij, thermo.ghs, t);
    vector<double> XA = solve_association_site_fractions_cpp(association.delta_ij, thermo.den, association.x_assoc);
    double value = 0.0;
    for (int i = 0; i < static_cast<int>(association.site_component_index.size()); ++i) {
        int component_index = association.site_component_index[i];
        value += x[component_index] * (std::log(XA[i]) - 0.5 * XA[i] + 0.5);
    }
    return value;
}

// EqID: ares_dh
double ares_ion_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs) {
    if (cppargs.z.empty()) {
        return 0.0;
    }

    int ncomp = static_cast<int>(x.size());
    vector<double> d(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        d[i] = cppargs.s[i] * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
        if (is_ion_species(cppargs, i)) {
            d[i] = ion_diameter_cpp(i, t, cppargs);
        }
    }

    double den = rho * N_AV / 1.0e30;
    double Qsum = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        Qsum += cppargs.z[i] * cppargs.z[i] * x[i];
    }
    if (Qsum == 0.0) {
        return 0.0;
    }

    DielectricState dielc_state = dielectric_state_cpp(x, cppargs);
    double eps = dielc_state.eps;
    double kappa = dh_kappa_cpp(den, t, eps, Qsum);
    if (kappa == 0.0) {
        return 0.0;
    }

    double S = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        S += x[i] * cppargs.z[i] * cppargs.z[i] * dh_chi_cpp(kappa, d[i]);
    }

    double K0 = E_CHRG * E_CHRG / (12.0 * PI * kb * t * perm_vac);
    return -K0 * kappa / eps * S;
}

// EqID: ares_born
double ares_born_cpp(double t, const vector<double> &x, const add_args &cppargs) {
    if (cppargs.z.empty() || cppargs.born_model == 0) {
        return 0.0;
    }

    double eps_mix = dielectric_constant_rule_cpp(cppargs.dielc_rule, x, cppargs);
    double eps_born = (cppargs.born_eps_mode == 1)
        ? reference_solvent_dielectric_constant_cpp(x, cppargs)
        : eps_mix;

    if (cppargs.born_model == 1) {
        double born_sum = 0.0;
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            if (is_ion_species(cppargs, i)) {
                double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                born_sum += x[i] * cppargs.z[i] * cppargs.z[i] / d_born_i;
            }
        }
        return -E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac) * (1.0 - 1.0 / eps_born) * born_sum;
    }
    if (cppargs.born_model == 2) {
        const double eps_r_ion = 8.0;
        const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        BornSSMDSData born = born_shell_data_cpp(x, cppargs, t, eps_born, eps_r_ion);
        return -Kborn * born.sum_bracket;
    }
    throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
}

}  // namespace

double ares_contribution_value_cpp(const AresContributions &terms, AresContributionKind kind) {
    switch (kind) {
        case AresContributionKind::HC:
            return terms.hc;
        case AresContributionKind::DISP:
            return terms.disp;
        case AresContributionKind::POLAR:
            return terms.polar;
        case AresContributionKind::ASSOC:
            return terms.assoc;
        case AresContributionKind::ION:
            return terms.ion;
        case AresContributionKind::BORN:
            return terms.born;
    }
    throw ValueError("Unknown AresContributionKind.");
}

// EqID: ares_total
AresContributions ares_contributions_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs) {
    AresContributions out;
    ThermoCommonState thermo = thermo_common_state_cpp(t, rho, x, cppargs, false);
    out.hc = ares_hc_cpp(thermo, x, cppargs);
    out.disp = ares_disp_cpp(thermo);
    out.polar = ares_polar_cpp(thermo, t, x, cppargs);
    out.assoc = ares_assoc_cpp(thermo, t, x, cppargs);
    out.ion = ares_ion_cpp(t, rho, x, cppargs);
    out.born = ares_born_cpp(t, x, cppargs);
    return out;
}

ScalarContributionTerms residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    AresContributions contributions = ares_contributions_cpp(t, rho, x, cppargs);
    double ares = contributions.hc + contributions.disp + contributions.polar + contributions.assoc + contributions.ion + contributions.born;
    return make_scalar_terms(
        contributions.hc,
        contributions.disp,
        contributions.polar,
        contributions.assoc,
        contributions.ion,
        contributions.born,
        ares
    );
}

double ares_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    AresContributions contributions = ares_contributions_cpp(t, rho, x, cppargs);
    return contributions.hc + contributions.disp + contributions.polar + contributions.assoc + contributions.ion + contributions.born;
}
