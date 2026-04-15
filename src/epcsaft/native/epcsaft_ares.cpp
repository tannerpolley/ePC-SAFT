#include "epcsaft_core_internal.h"

using namespace thermo_detail;

// EqID: ares_born
double born_ares_only_cpp(double t, const vector<double> &x, const add_args &cppargs) {
    if (cppargs.born_model == 0) {
        return 0.0;
    }
    double eps_mix = dielectric_constant_rule_cpp(cppargs.dielc_rule, x, cppargs);
    double eps_born = (cppargs.born_eps_mode == 1) ? reference_solvent_dielectric_constant_cpp(x, cppargs) : eps_mix;
    if (cppargs.born_model == 1) {
        double born_sum = 0.0;
        for (int i = 0; i < static_cast<int>(x.size()); i++) {
            if (is_ion_species(cppargs, i)) {
                double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
            }
        }
        return -E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac)*(1.0 - 1.0/eps_born)*born_sum;
    }
    if (cppargs.born_model == 2) {
        const double eps_r_ion = 8.0;
        const double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
        BornSSMDSData born = born_shell_data_cpp(x, cppargs, t, eps_born, eps_r_ion);
        return -Kborn*born.sum_bracket;
    }
    throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
}

// EqID: ares_dh
double dh_ares_only_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs) {
    if (cppargs.z.empty()) {
        return 0.0;
    }
    int ncomp = static_cast<int>(x.size());
    vector<double> d(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++) {
        d[i] = cppargs.s[i]*(1.0 - 0.12*std::exp(-3.0*cppargs.e[i]/t));
        if (is_ion_species(cppargs, i)) {
            d[i] = ion_diameter_cpp(i, t, cppargs);
        }
    }

    double den = rho*N_AV/1.0e30;
    double Qsum = 0.0;
    for (int i = 0; i < ncomp; i++) {
        Qsum += cppargs.z[i]*cppargs.z[i]*x[i];
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
    for (int i = 0; i < ncomp; i++) {
        S += x[i]*cppargs.z[i]*cppargs.z[i]*dh_chi_cpp(kappa, d[i]);
    }

    double K0 = E_CHRG*E_CHRG/(12.0*PI*kb*t*perm_vac);
    return -K0*kappa/eps*S;
}

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
// EqID: ares_hc
// EqID: ares_disp
// EqID: ares_assoc
// EqID: ares_dh
// EqID: ares_born
AresContributions ares_contributions_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs) {
    AresContributions out;
    int ncomp = static_cast<int>(x.size());
    ThermoCommonState thermo = thermo_common_state_cpp(t, rho, x, cppargs, false);
    auto &zeta = thermo.zeta;
    auto &e_ij = thermo.e_ij;
    auto &s_ij = thermo.s_ij;
    auto &ghs = thermo.ghs;
    double den = thermo.den;
    double eta = thermo.eta;
    double m_avg = thermo.m_avg;
    double m2es3 = thermo.m2es3;
    double m2e2s3 = thermo.m2e2s3;
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(m_avg, eta);
    double ares_hs = 1.0 / zeta[0] * (3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3])
            + std::pow(zeta[2], 3.0) / (zeta[3] * std::pow(1.0 - zeta[3], 2.0))
            + (std::pow(zeta[2], 3.0) / std::pow(zeta[3], 2.0) - zeta[0]) * std::log(1.0 - zeta[3]));

    double summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i]*(cppargs.m[i]-1)*log(ghs[i*ncomp+i]);
    }
    out.hc = m_avg*ares_hs - summ;
    out.disp = -2*PI*den*dispersion.I1*m2es3 - PI*den*m_avg*dispersion.C1*dispersion.I2*m2e2s3;

    if (!cppargs.dipm.empty()) {
        double A2 = 0.0;
        double A3 = 0.0;
        vector<double> dipmSQ(ncomp, 0.0);
        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(cppargs.dipm[i], 2.)/(cppargs.m[i]*cppargs.e[i]*pow(cppargs.s[i],3.))*kDipoleConversion;
        }

        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                double m_ij = sqrt(cppargs.m[i]*cppargs.m[j]);
                if (m_ij > 2) {
                    m_ij = 2;
                }
                vector<double> adip = dipole_coefficients_cpp(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = dipole_coefficients_cpp(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
                double J2 = 0.0;
                for (int l = 0; l < 5; l++) {
                    J2 += (adip[l] + bdip[l]*e_ij[j*ncomp+j]/t)*pow(eta, l);
                }
                A2 += x[i]*x[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                    pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*J2;

                for (int k = 0; k < ncomp; k++) {
                    double m_ijk = pow((cppargs.m[i]*cppargs.m[j]*cppargs.m[k]),1/3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    vector<double> cdip = dipole_coefficients_cpp(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
                    double J3 = 0.0;
                    for (int l = 0; l < 5; l++) {
                        J3 += cdip[l]*pow(eta, l);
                    }
                    A3 += x[i]*x[j]*x[k]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*J3;
                }
            }
        }

        A2 = -PI*den*A2;
        A3 = -4/3.*PI*PI*den*den*A3;
        if (A2 != 0) {
            out.polar = A2/(1-A3/A2);
        }
    }

    if (!cppargs.e_assoc.empty()) {
        AssociationSetup association = association_setup_cpp(x, cppargs, s_ij, ghs, t);
        vector<double> XA = solve_association_site_fractions_cpp(association.delta_ij, den, association.x_assoc);
        for (int i = 0; i < static_cast<int>(association.site_component_index.size()); i++) {
            int component_index = association.site_component_index[i];
            out.assoc += x[component_index]*(log(XA[i])-0.5*XA[i] + 0.5);
        }
    }

    if (!cppargs.z.empty()) {
        DielectricState dielc_state = dielectric_state_cpp(x, cppargs);
        double eps = dielc_state.eps;
        double eps_born = (cppargs.born_eps_mode == 1) ? reference_solvent_dielectric_constant_cpp(x, cppargs) : eps;
        out.ion = dh_ares_only_cpp(t, rho, x, cppargs);

        if (cppargs.born_model == 1) {
            double born_sum = 0.0;
            for (int i = 0; i < ncomp; i++) {
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                    born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
                }
            }
            out.born = -E_CHRG*E_CHRG/(4.*PI*kb*t*perm_vac)*(1.-1./eps_born)*born_sum;
        }
        else if (cppargs.born_model == 2) {
            const double eps_r_ion = 8.0;
            const double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
            BornSSMDSData born = born_shell_data_cpp(x, cppargs, t, eps_born, eps_r_ion);
            out.born = -Kborn*born.sum_bracket;
        }
        else if (cppargs.born_model != 0) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }

    return out;
}

// EqID: ares_total
ScalarContributionTerms residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size()); // number of components
    double summ = 0.0;
    ThermoCommonState thermo = thermo_common_state_cpp(t, rho, x, cppargs, false);
    auto &d = thermo.d;
    auto &zeta = thermo.zeta;
    auto &e_ij = thermo.e_ij;
    auto &s_ij = thermo.s_ij;
    auto &ghs = thermo.ghs;
    double den = thermo.den;
    double eta = thermo.eta;
    double m_avg = thermo.m_avg;
    double m2es3 = thermo.m2es3;
    double m2e2s3 = thermo.m2e2s3;
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(m_avg, eta);
    double ares_hs = 1.0 / zeta[0] * (3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3])
        + std::pow(zeta[2], 3.0) / (zeta[3] * std::pow(1.0 - zeta[3], 2.0))
        + (std::pow(zeta[2], 3.0) / std::pow(zeta[3], 2.0) - zeta[0]) * std::log(1.0 - zeta[3]));
    double I1 = dispersion.I1;
    double I2 = dispersion.I2;
    double C1 = dispersion.C1;

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i]*(cppargs.m[i]-1)*log(ghs[i*ncomp+i]);
    }

    double ares_hc = m_avg*ares_hs - summ;
    double ares_disp = -2*PI*den*I1*m2es3 - PI*den*m_avg*C1*I2*m2e2s3;

    // Dipole term (Gross and Vrabec term) --------------------------------------
    double ares_polar = 0.;
    if (!cppargs.dipm.empty()) {
        double A2 = 0.;
        double A3 = 0.;
        vector<double> dipmSQ (ncomp, 0);
        const static double conv = 7242.702976750923; // conversion factor, see the note below Table 2 in Gross and Vrabec 2006

        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(cppargs.dipm[i], 2.)/(cppargs.m[i]*cppargs.e[i]*pow(cppargs.s[i],3.))*conv;
        }

        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                double m_ij = sqrt(cppargs.m[i]*cppargs.m[j]);
                if (m_ij > 2) {
                    m_ij = 2;
                }
                vector<double> adip = dipole_coefficients_cpp(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = dipole_coefficients_cpp(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
                double J2 = 0.;
                for (int l = 0; l < 5; l++) {
                    J2 += (adip[l] + bdip[l]*e_ij[j*ncomp+j]/t)*pow(eta, l); // j*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                }
                A2 += x[i]*x[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                    pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*J2;

                for (int k = 0; k < ncomp; k++) {
                    double m_ijk = pow((cppargs.m[i]*cppargs.m[j]*cppargs.m[k]),1/3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    vector<double> cdip = dipole_coefficients_cpp(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
                    double J3 = 0.;
                    for (int l = 0; l < 5; l++) {
                        J3 += cdip[l]*pow(eta, l);
                    }
                    A3 += x[i]*x[j]*x[k]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*J3;
                }
            }
        }

        A2 = -PI*den*A2;
        A3 = -4/3.*PI*PI*den*den*A3;

        if (A2 != 0) { // when the mole fraction of the polar compounds is 0 then A2 = 0 and division by 0 occurs
            ares_polar = A2/(1-A3/A2);
        }
    }

    // Association term -------------------------------------------------------
    double ares_assoc = 0.;
    if (!cppargs.e_assoc.empty()) {
        AssociationSetup assoc = association_setup_cpp(x, cppargs, s_ij, ghs, t);
        const vector<int> &iA = assoc.site_component_index;
        const vector<double> &x_assoc = assoc.x_assoc;
        const vector<double> &delta_ij = assoc.delta_ij;
        int num_sites = static_cast<int>(iA.size());

        vector<double> XA = solve_association_site_fractions_cpp(delta_ij, den, x_assoc);

        ares_assoc = 0.;
        for (int i = 0; i < num_sites; i++) {
            ares_assoc += x[iA[i]]*(log(XA[i])-0.5*XA[i] + 0.5);
        }
    }

    // Ion term ---------------------------------------------------------------
    double ares_ion = 0.;
    double ares_born = 0.;
    if (!cppargs.z.empty()) {
        DielectricState dielc_state = dielectric_state_cpp(x, cppargs);
        double eps = dielc_state.eps;
        double eps_born = (cppargs.born_eps_mode == 1) ? reference_solvent_dielectric_constant_cpp(x, cppargs) : eps;
        ares_ion = dh_ares_only_cpp(t, rho, x, cppargs);

        if (cppargs.born_model == 1) {
            // Born term (Bulow 2021a, non-SSM+DS): use d_born,i in denominator
            double born_sum = 0.;
            for (int i = 0; i < ncomp; i++) {
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                    born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
                }
            }
            ares_born = -E_CHRG*E_CHRG/(4.*PI*kb*t*perm_vac)*(1.-1./eps_born)*born_sum;
        }
        else if (cppargs.born_model == 2) {
            const double eps_r_ion = 8.0;
            const double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
            BornSSMDSData born = born_shell_data_cpp(x, cppargs, t, eps_born, eps_r_ion);
            ares_born = -Kborn*born.sum_bracket;
        }
        else if (cppargs.born_model != 0) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }

    double ares = ares_hc + ares_disp + ares_polar + ares_assoc + ares_ion + ares_born;
    return make_scalar_terms(ares_hc, ares_disp, ares_polar, ares_assoc, ares_ion, ares_born, ares);
}


double ares_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return residual_helmholtz_result_cpp(t, rho, x, cppargs).total;
}
