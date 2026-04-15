#include "epcsaft_core_internal.h"

using namespace thermo_detail;

vector<double> contribution_dadx_fd_cpp(AresContributionKind kind, double t, double rho, const vector<double> &x, const add_args &cppargs, double a0) {
    int ncomp = static_cast<int>(x.size());
    vector<double> dadx(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++) {
        double h = 1e-6*std::max(1.0, std::abs(x[i]));
        vector<double> xp = x;
        xp[i] += h;
        double fp = ares_contribution_value_cpp(ares_contributions_cpp(t, rho, xp, cppargs), kind);
        if (x[i] - h >= 0.0) {
            vector<double> xm = x;
            xm[i] -= h;
            double fm = ares_contribution_value_cpp(ares_contributions_cpp(t, rho, xm, cppargs), kind);
            dadx[i] = (fp - fm)/(2.0*h);
        }
        else {
            dadx[i] = (fp - a0)/h;
        }
        if (!std::isfinite(dadx[i])) {
            throw ValueError("Non-finite contribution finite-difference derivative.");
        }
    }
    return dadx;
}

CompositionContributionResult composition_derivative_residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size()); // number of components
    vector<double> d (ncomp);
    for (int i = 0; i < ncomp; i++) {
        d[i] = cppargs.s[i]*(1-0.12*exp(-3*cppargs.e[i]/t));
        if (!cppargs.z.empty() && is_ion_species(cppargs, i)) {
            d[i] = ion_diameter_cpp(i, t, cppargs);
        }
    }

    double den = rho*N_AV/1.0e30;

    vector<double> zeta (4, 0);
    double summ;
    for (int i = 0; i < 4; i++) {
        summ = 0;
        for (int j = 0; j < ncomp; j++) {
            summ += x[j]*cppargs.m[j]*pow(d[j], i);
        }
        zeta[i] = PI/6*den*summ;
    }

    double eta = zeta[3];
    double m_avg = 0;
    for (int i = 0; i < ncomp; i++) {
        m_avg += x[i]*cppargs.m[i];
    }

    vector<double> ghs(ncomp*ncomp, 0);
    vector<double> denghs(ncomp*ncomp, 0);
    vector<double> e_ij(ncomp*ncomp, 0);
    vector<double> s_ij(ncomp*ncomp, 0);
    double m2es3 = 0.;
    double m2e2s3 = 0.;
    int idx = -1;
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < ncomp; j++) {
            idx += 1;
            if (cppargs.l_ij.empty()) {
                s_ij[idx] = (cppargs.s[i] + cppargs.s[j])/2.;
            }
            else {
                s_ij[idx] = (cppargs.s[i] + cppargs.s[j])/2.*(1-cppargs.l_ij[idx]);
            }
            if (!cppargs.z.empty()) {
                if (cppargs.z[i]*cppargs.z[j] <= 0) { // for two cations or two anions e_ij is kept at zero to avoid dispersion between like ions (see Held et al. 2014)
                    if (cppargs.k_ij.empty()) {
                        e_ij[idx] = sqrt(cppargs.e[i]*cppargs.e[j]);
                    }
                    else {
                        e_ij[idx] = sqrt(cppargs.e[i]*cppargs.e[j])*(1-cppargs.k_ij[idx]);
                    }
                }
            } else {
                if (cppargs.k_ij.empty()) {
                    e_ij[idx] = sqrt(cppargs.e[i]*cppargs.e[j]);
                }
                else {
                    e_ij[idx] = sqrt(cppargs.e[i]*cppargs.e[j])*(1-cppargs.k_ij[idx]);
                }
            }
            m2es3 = m2es3 + x[i]*x[j]*cppargs.m[i]*cppargs.m[j]*e_ij[idx]/t*pow(s_ij[idx], 3);
            m2e2s3 = m2e2s3 + x[i]*x[j]*cppargs.m[i]*cppargs.m[j]*pow(e_ij[idx]/t,2)*pow(s_ij[idx], 3);
            ghs[idx] = 1/(1-zeta[3]) + (d[i]*d[j]/(d[i]+d[j]))*3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
                    pow(d[i]*d[j]/(d[i]+d[j]), 2)*2*zeta[2]*zeta[2]/pow(1-zeta[3], 3);
            denghs[idx] = zeta[3]/(1-zeta[3])/(1-zeta[3]) +
                (d[i]*d[j]/(d[i]+d[j]))*(3*zeta[2]/(1-zeta[3])/(1-zeta[3]) +
                6*zeta[2]*zeta[3]/pow(1-zeta[3], 3)) +
                pow(d[i]*d[j]/(d[i]+d[j]), 2)*(4*zeta[2]*zeta[2]/pow(1-zeta[3], 3) +
                6*zeta[2]*zeta[2]*zeta[3]/pow(1-zeta[3], 4));
        }
    }

    double ares_hs = 1/zeta[0]*(3*zeta[1]*zeta[2]/(1-zeta[3]) + pow(zeta[2], 3.)/(zeta[3]*pow(1-zeta[3],2))
            + (pow(zeta[2], 3.)/pow(zeta[3], 2.) - zeta[0])*log(1-zeta[3]));

    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(m_avg, eta);
    const auto &a = dispersion.a;
    const auto &b = dispersion.b;
    double I1 = dispersion.I1;
    double I2 = dispersion.I2;
    double C1 = dispersion.C1;
    double C2 = dispersion.C2;

    DadrhoResult dadrho_result = dadrho_result_cpp(t, rho, x, cppargs);
    const ScalarContributionTerms &z_raw_terms = dadrho_result.raw;
    ScalarContributionTerms z_terms = compressibility_terms_from_dadrho_cpp(dadrho_result);
    const double Zhc = z_raw_terms.hc;
    const double Zdisp = z_raw_terms.disp;
    const double Zpolar = z_raw_terms.polar;
    const double Zassoc = z_raw_terms.assoc;
    const double Zion = z_raw_terms.ion;
    const double Zborn = z_raw_terms.born;
    const double Z = z_terms.total;

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i]*(cppargs.m[i]-1)*log(ghs[i*ncomp+i]);
    }

    double ares_hc = m_avg*ares_hs - summ;
    double ares_disp = -2*PI*den*I1*m2es3 - PI*den*m_avg*C1*I2*m2e2s3;

    vector<double> dghsii_dx(ncomp*ncomp, 0);
    vector<double> dahs_dx(ncomp, 0);
    vector<double> dzeta_dx(4, 0);
    idx = -1;
    for (int i = 0; i < ncomp; i++) {
        for (int l = 0; l < 4; l++) {
            dzeta_dx[l] = PI/6.*den*cppargs.m[i]*pow(d[i],l);
        }
        for (int j = 0; j < ncomp; j++) {
            idx += 1;
            dghsii_dx[idx] = dzeta_dx[3]/(1-zeta[3])/(1-zeta[3]) + (d[j]*d[j]/(d[j]+d[j]))*
                    (3*dzeta_dx[2]/(1-zeta[3])/(1-zeta[3]) + 6*zeta[2]*dzeta_dx[3]/pow(1-zeta[3],3))
                    + pow(d[j]*d[j]/(d[j]+d[j]),2)*(4*zeta[2]*dzeta_dx[2]/pow(1-zeta[3],3)
                    + 6*zeta[2]*zeta[2]*dzeta_dx[3]/pow(1-zeta[3],4));
        }
        dahs_dx[i] = -dzeta_dx[0]/zeta[0]*ares_hs + 1/zeta[0]*(3*(dzeta_dx[1]*zeta[2]
                + zeta[1]*dzeta_dx[2])/(1-zeta[3]) + 3*zeta[1]*zeta[2]*dzeta_dx[3]
                /(1-zeta[3])/(1-zeta[3]) + 3*zeta[2]*zeta[2]*dzeta_dx[2]/zeta[3]/(1-zeta[3])/(1-zeta[3])
                + pow(zeta[2],3)*dzeta_dx[3]*(3*zeta[3]-1)/zeta[3]/zeta[3]/pow(1-zeta[3],3)
                + log(1-zeta[3])*((3*zeta[2]*zeta[2]*dzeta_dx[2]*zeta[3] -
                2*pow(zeta[2],3)*dzeta_dx[3])/pow(zeta[3],3) - dzeta_dx[0]) +
                (zeta[0]-pow(zeta[2],3)/zeta[3]/zeta[3])*dzeta_dx[3]/(1-zeta[3]));
    }

    vector<double> dadisp_dx(ncomp, 0);
    vector<double> dahc_dx(ncomp, 0);
    double dzeta3_dx, daa_dx, db_dx, dI1_dx, dI2_dx, dm2es3_dx, dm2e2s3_dx, dC1_dx;
    for (int i = 0; i < ncomp; i++) {
        dzeta3_dx = PI/6.*den*cppargs.m[i]*pow(d[i],3);
        dI1_dx = 0.0;
        dI2_dx = 0.0;
        dm2es3_dx = 0.0;
        dm2e2s3_dx = 0.0;
        for (int l = 0; l < 7; l++) {
            daa_dx = cppargs.m[i]/m_avg/m_avg*kDispersionA1[l] + cppargs.m[i]/m_avg/m_avg*(3-4/m_avg)*kDispersionA2[l];
            db_dx = cppargs.m[i]/m_avg/m_avg*kDispersionB1[l] + cppargs.m[i]/m_avg/m_avg*(3-4/m_avg)*kDispersionB2[l];
            dI1_dx += a[l]*l*dzeta3_dx*pow(eta,l-1) + daa_dx*pow(eta,l);
            dI2_dx += b[l]*l*dzeta3_dx*pow(eta,l-1) + db_dx*pow(eta,l);
        }
        for (int j = 0; j < ncomp; j++) {
            dm2es3_dx += x[j]*cppargs.m[j]*(e_ij[i*ncomp+j]/t)*pow(s_ij[i*ncomp+j],3);
            dm2e2s3_dx += x[j]*cppargs.m[j]*pow(e_ij[i*ncomp+j]/t,2)*pow(s_ij[i*ncomp+j],3);
            dahc_dx[i] += x[j]*(cppargs.m[j]-1)/ghs[j*ncomp+j]*dghsii_dx[i*ncomp+j];
        }
        dm2es3_dx = dm2es3_dx*2*cppargs.m[i];
        dm2e2s3_dx = dm2e2s3_dx*2*cppargs.m[i];
        dahc_dx[i] = cppargs.m[i]*ares_hs + m_avg*dahs_dx[i] - dahc_dx[i] - (cppargs.m[i]-1)*log(ghs[i*ncomp+i]);
        dC1_dx = C2*dzeta3_dx - C1*C1*(cppargs.m[i]*(8*eta-2*eta*eta)/pow(1-eta,4) -
            cppargs.m[i]*(20*eta-27*eta*eta+12*pow(eta,3)-2*pow(eta,4))/pow((1-eta)*(2-eta),2));

        dadisp_dx[i] = -2*PI*den*(dI1_dx*m2es3 + I1*dm2es3_dx) - PI*den
            *((cppargs.m[i]*C1*I2 + m_avg*dC1_dx*I2 + m_avg*C1*dI2_dx)*m2e2s3
            + m_avg*C1*I2*dm2e2s3_dx);
    }

    if (cppargs.hc_dadx_diff_mode == 1) {
        dahc_dx = contribution_dadx_fd_cpp(AresContributionKind::HC, t, rho, x, cppargs, ares_hc);
    }
    if (cppargs.disp_dadx_diff_mode == 1) {
        dadisp_dx = contribution_dadx_fd_cpp(AresContributionKind::DISP, t, rho, x, cppargs, ares_disp);
    }

    vector<double> mu_hc(ncomp, 0);
    vector<double> mu_disp(ncomp, 0);
    double sum_x_dahc_dx = 0.0;
    double sum_x_dadisp_dx = 0.0;
    for (int i = 0; i < ncomp; i++) {
        sum_x_dahc_dx += x[i]*dahc_dx[i];
        sum_x_dadisp_dx += x[i]*dadisp_dx[i];
    }
    for (int i = 0; i < ncomp; i++) {
        mu_hc[i] = ares_hc + Zhc + dahc_dx[i] - sum_x_dahc_dx;
        mu_disp[i] = ares_disp + Zdisp + dadisp_dx[i] - sum_x_dadisp_dx;
    }

    // Dipole term (Gross and Vrabec term) --------------------------------------
    vector<double> mu_polar(ncomp, 0);
    vector<double> dapolar_dx(ncomp, 0.0);
    double ares_polar = 0.0;
    double sum_x_dapolar_dx = 0.0;
    if (!cppargs.dipm.empty()) {
        double A2 = 0.;
        double A3 = 0.;
        double dA2_det = 0.;
        double dA3_det = 0.;
        vector<double> dA2_dx(ncomp, 0);
        vector<double> dA3_dx(ncomp, 0);

        const static double conv = 7242.702976750923; // conversion factor, see the note below Table 2 in Gross and Vrabec 2006

        vector<double> dipmSQ (ncomp, 0);
        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(cppargs.dipm[i], 2.)/(cppargs.m[i]*cppargs.e[i]*pow(cppargs.s[i],3.))*conv;
        }

        vector<double> adip (5, 0);
        vector<double> bdip (5, 0);
        vector<double> cdip (5, 0);
        double J2, dJ2_det, detJ2_det, J3, dJ3_det, detJ3_det;
        double m_ij;
        double m_ijk;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                m_ij = sqrt(cppargs.m[i]*cppargs.m[j]);
                if (m_ij > 2) {
                    m_ij = 2;
                }
                vector<double> adip = dipole_coefficients_cpp(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = dipole_coefficients_cpp(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
                J2 = 0.;
                dJ2_det = 0.;
                detJ2_det = 0;
                for (int l = 0; l < 5; l++) {
                    J2 += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/t)*pow(eta, l); // i*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    dJ2_det += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/t)*l*pow(eta, l-1);
                    detJ2_det += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/t)*(l+1)*pow(eta, l);
                }
                A2 += x[i]*x[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                    pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*J2;
                dA2_det += x[i]*x[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*
                    pow(s_ij[j*ncomp+j],3)/pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*detJ2_det;
                if (i == j) {
                    dA2_dx[i] += e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)
                        /pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*
                        (x[i]*x[j]*dJ2_det*PI/6.*den*cppargs.m[i]*pow(d[i],3) + 2*x[j]*J2);
                }
                else {
                    dA2_dx[i] += e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)
                        /pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*
                        (x[i]*x[j]*dJ2_det*PI/6.*den*cppargs.m[i]*pow(d[i],3) + x[j]*J2);
                }

                for (int k = 0; k < ncomp; k++) {
                    m_ijk = pow((cppargs.m[i]*cppargs.m[j]*cppargs.m[k]),1/3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    vector<double> cdip = dipole_coefficients_cpp(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
                    J3 = 0.;
                    dJ3_det = 0.;
                    detJ3_det = 0.;
                    for (int l = 0; l < 5; l++) {
                        J3 += cdip[l]*pow(eta, l);
                        dJ3_det += cdip[l]*l*pow(eta, (l-1));
                        detJ3_det += cdip[l]*(l+2)*pow(eta, (l+1));
                    }
                    A3 += x[i]*x[j]*x[k]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*J3;
                    dA3_det += x[i]*x[j]*x[k]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*detJ3_det;
                    if ((i == j) && (i == k)) {
                        dA3_dx[i] += e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*pow(s_ij[i*ncomp+i],3)
                            *pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/s_ij[j*ncomp+k]
                            *cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]*dipmSQ[j]
                            *dipmSQ[k]*(x[i]*x[j]*x[k]*dJ3_det*PI/6.*den*cppargs.m[i]*pow(d[i],3)
                            + 3*x[j]*x[k]*J3);
                    }
                    else if ((i == j) || (i == k)) {
                        dA3_dx[i] += e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*pow(s_ij[i*ncomp+i],3)
                            *pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/s_ij[j*ncomp+k]
                            *cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]*dipmSQ[j]
                            *dipmSQ[k]*(x[i]*x[j]*x[k]*dJ3_det*PI/6.*den*cppargs.m[i]*pow(d[i],3)
                            + 2*x[j]*x[k]*J3);
                    }
                    else {
                        dA3_dx[i] += e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*pow(s_ij[i*ncomp+i],3)
                            *pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/s_ij[j*ncomp+k]
                            *cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]*dipmSQ[j]
                            *dipmSQ[k]*(x[i]*x[j]*x[k]*dJ3_det*PI/6.*den*cppargs.m[i]*pow(d[i],3)
                            + x[j]*x[k]*J3);
                    }
                }
            }
        }

        A2 = -PI*den*A2;
        A3 = -4/3.*PI*PI*den*den*A3;
        dA2_det = -PI*den/eta*dA2_det;
        dA3_det = -4/3.*PI*PI*den/eta*den/eta*dA3_det;
        for (int i = 0; i < ncomp; i++) {
            dA2_dx[i] = -PI*den*dA2_dx[i];
            dA3_dx[i] = -4/3.*PI*PI*den*den*dA3_dx[i];
        }

        if (A2 != 0) { // when the mole fraction of the polar compounds is 0 then A2 = 0 and division by 0 occurs
            ares_polar = A2/(1-A3/A2);
            for (int i = 0; i < ncomp; i++) {
                dapolar_dx[i] = (dA2_dx[i]*(1-A3/A2) + (dA3_dx[i]*A2 - A3*dA2_dx[i])/A2)/pow(1-A3/A2,2);
            }
            if (cppargs.polar_dadx_diff_mode == 1) {
                dapolar_dx = contribution_dadx_fd_cpp(AresContributionKind::POLAR, t, rho, x, cppargs, ares_polar);
            }
            for (int i = 0; i < ncomp; i++) {
                sum_x_dapolar_dx += x[i]*dapolar_dx[i];
            }
            for (int i = 0; i < ncomp; i++) {
                mu_polar[i] = ares_polar + Zpolar + dapolar_dx[i] - sum_x_dapolar_dx;
            }
        }
    }

    // Association term -------------------------------------------------------
    vector<double> mu_assoc(ncomp, 0);
    vector<double> daassoc_dx(ncomp, 0.0);
    double ares_assoc = 0.0;
    double sum_x_daassoc_dx = 0.0;
    if (!cppargs.e_assoc.empty()) {
        AssociationSetup assoc = association_setup_cpp(x, cppargs, s_ij, ghs, t);
        const vector<int> &iA = assoc.site_component_index;
        const vector<double> &x_assoc = assoc.x_assoc;
        const vector<double> &delta_ij = assoc.delta_ij;
        int num_sites = static_cast<int>(iA.size());

        vector<double> ddelta_dx(num_sites * num_sites * ncomp, 0);
        int idx_ddelta = 0;
        for (int k = 0; k < ncomp; k++) {
            for (int i = 0; i < num_sites; i++) {
                int idxi = iA[i]*ncomp+iA[i];
                for (int j = 0; j < num_sites; j++) {
                    int idxj = iA[j]*ncomp+iA[j];
                    if (cppargs.assoc_matrix[i*num_sites+j] != 0) {
                        double eABij = (cppargs.e_assoc[iA[i]]+cppargs.e_assoc[iA[j]])/2.;
                        double volABij = HUGE_DBL;
                        if (cppargs.k_hb.empty()) {
                            volABij = sqrt(cppargs.vol_a[iA[i]]*cppargs.vol_a[iA[j]])*pow(sqrt(s_ij[idxi]*
                                s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                        }
                        else {
                            volABij = sqrt(cppargs.vol_a[iA[i]]*cppargs.vol_a[iA[j]])*pow(sqrt(s_ij[idxi]*
                                s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3)*(1-cppargs.k_hb[iA[i]*ncomp+iA[j]]);
                        }
                        double dghsd_dx = PI/6.*cppargs.m[k]*(pow(d[k], 3)/(1-zeta[3])/(1-zeta[3]) + 3*d[iA[i]]*d[iA[j]]/
                            (d[iA[i]]+d[iA[j]])*(d[k]*d[k]/(1-zeta[3])/(1-zeta[3])+2*pow(d[k], 3)*
                            zeta[2]/pow(1-zeta[3], 3)) + 2*pow((d[iA[i]]*d[iA[j]]/(d[iA[i]]+d[iA[j]])), 2)*
                            (2*d[k]*d[k]*zeta[2]/pow(1-zeta[3], 3)+3*(pow(d[k], 3)*zeta[2]*zeta[2]
                            /pow(1-zeta[3], 4))));
                        ddelta_dx[idx_ddelta] = dghsd_dx*(exp(eABij/t)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij;
                    }
                    idx_ddelta += 1;
                }
            }
        }

        vector<double> XA = solve_association_site_fractions_cpp(delta_ij, den, x_assoc);

        vector<double> dXA_dx(num_sites*ncomp, 0);
        dXA_dx = ::association_site_fraction_dx_cpp(cppargs.assoc_num, delta_ij, den, XA, ddelta_dx, x_assoc);

        int ij = 0;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < num_sites; j++) {
                daassoc_dx[i] += x[iA[j]]*den*dXA_dx[ij]*(1/XA[j]-0.5);
                ij += 1;
            }
        }

        for (int i = 0; i < num_sites; i++) {
            daassoc_dx[iA[i]] += log(XA[i]) - 0.5*XA[i] + 0.5;
            ares_assoc += x[iA[i]]*(log(XA[i]) - 0.5*XA[i] + 0.5);
        }
        if (cppargs.assoc_dadx_diff_mode == 1) {
            daassoc_dx = contribution_dadx_fd_cpp(AresContributionKind::ASSOC, t, rho, x, cppargs, ares_assoc);
        }
        for (int i = 0; i < ncomp; i++) {
            sum_x_daassoc_dx += x[i]*daassoc_dx[i];
        }
        for (int i = 0; i < ncomp; i++) {
            mu_assoc[i] = ares_assoc + Zassoc + daassoc_dx[i] - sum_x_daassoc_dx;
        }
    }

    // Ion terms --------------------------------------------------------------
    vector<double> mu_ion(ncomp, 0);
    vector<double> mu_born(ncomp, 0);
    vector<double> dadx_ion(ncomp, 0.0);
    vector<double> dadx_born_diag(ncomp, 0.0);
    double a_ion = 0.0;
    double a_born_diag = 0.0;
    double sum_x_dadx_ion = 0.0;
    double sum_x_dadx_born_diag = 0.0;
    if (!cppargs.z.empty()) {
        int dh_model = cppargs.DH_model;
        if (dh_model == 2) {
            throw ValueError("DH_model=2 (Bjerrum treatment) is reserved and not implemented.");
        }
        if ((dh_model != 0) && (dh_model != 1)) {
            throw ValueError("Unknown DH_model. Supported values are 0, 1, and reserved 2.");
        }

        // Debye-Huckel term
        double Qsum = 0.;
        for (int i = 0; i < ncomp; i++) {
            Qsum += cppargs.z[i]*cppargs.z[i]*x[i];
        }
        DielectricState dielc_state = dielectric_state_cpp(x, cppargs);
        double eps = dielc_state.eps; // mixed dielectric constant (relative)
        vector<double> deps_dx = dielc_state.deps_dx; // d(eps_r)/dx_i
        double eps_born = eps;
        vector<double> deps_dx_born = deps_dx;
        if ((cppargs.born_model >= 1) && (cppargs.born_eps_mode == 1)) {
            eps_born = reference_solvent_dielectric_constant_cpp(x, cppargs);
            deps_dx_born = reference_solvent_dielectric_derivative_cpp(x, cppargs);
        }

        double kappa = dh_kappa_cpp(den, t, eps, Qsum); // inverse Debye screening length
        if ((kappa != 0) && (Qsum != 0)) {
            vector<double> chi(ncomp, 0.0);
            vector<double> sigma_k(ncomp, 0.0);
            double S = 0.;
            double Tsum = 0.;

            for (int i = 0; i < ncomp; i++) {
                double ka = kappa*d[i];
                chi[i] = dh_chi_cpp(kappa, d[i]);
                sigma_k[i] = -2*chi[i] + 3/(1+ka);

                S += x[i]*cppargs.z[i]*cppargs.z[i]*chi[i];
                Tsum += x[i]*cppargs.z[i]*cppargs.z[i]*sigma_k[i];
            }

            double K0 = E_CHRG*E_CHRG/(12.0*PI*kb*t*perm_vac); // without epsilon
            double a_DH = -K0*kappa/(eps)*S;
            double Z_DH = Zion;
            a_ion = a_DH;

            vector<double> dadx(ncomp, 0.0);
            vector<double> dkappa_dx(ncomp, 0.0);
            vector<double> dS_dx(ncomp, 0.0);
            const bool use_dh_deps = (cppargs.mu_DH_comp_dep_rel_perm != 0);
            const double dh_deps_multiplier = (cppargs.mu_DH_include_sum_term != 0) ? Qsum : 1.0;
            if (cppargs.mu_DH_diff_mode == 1) {
                dadx = contribution_dadx_fd_cpp(AresContributionKind::ION, t, rho, x, cppargs, a_DH);
            }
            else {
                double Aconst = den*E_CHRG*E_CHRG/(kb*t*perm_vac);
                for (int i = 0; i < ncomp; i++) {
                    double deps_term = use_dh_deps ? dh_deps_multiplier*deps_dx[i]/(eps*eps) : 0.0;
                    dkappa_dx[i] = Aconst*(cppargs.z[i]*cppargs.z[i]/eps - deps_term)/(2.0*kappa);
                }

                for (int i = 0; i < ncomp; i++) {
                    dS_dx[i] = cppargs.z[i]*cppargs.z[i]*chi[i] + dkappa_dx[i]*(Tsum - S)/kappa;
                }

                for (int i = 0; i < ncomp; i++) {
                    double d_inv_eps_dx = use_dh_deps ? -deps_dx[i]/(eps*eps) : 0.0;
                    double term1 = (dkappa_dx[i]/eps + kappa*d_inv_eps_dx)*S;
                    double term2 = kappa/eps*dS_dx[i];
                    dadx[i] = -K0*(term1 + term2);
                }
            }

            double sum_x_dadx = 0.0;
            for (int i = 0; i < ncomp; i++) {
                sum_x_dadx += x[i]*dadx[i];
            }
            dadx_ion = dadx;
            sum_x_dadx_ion = sum_x_dadx;

            for (int i = 0; i < ncomp; i++) {
                mu_ion[i] = a_DH + Z_DH + dadx[i] - sum_x_dadx;
            }
        }

        if (cppargs.born_model == 1) {
            // Born term (Bulow 2021a, non-SSM+DS), using user-selected born_radius_model for D_born
            double born_sum = 0.;
            for (int i = 0; i < ncomp; i++) {
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                    born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
                }
            }
            double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
            double a_born = -Kborn*(1.0 - 1.0/eps_born)*born_sum;

            vector<double> dadx_born(ncomp, 0.0);
            vector<double> ion_part_vec(ncomp, 0.0);
            vector<double> eps_part_vec(ncomp, 0.0);
            if (cppargs.born_diff_mode == 1) {
                dadx_born = contribution_dadx_fd_cpp(AresContributionKind::BORN, t, rho, x, cppargs, a_born);
            }
            else {
                for (int i = 0; i < ncomp; i++) {
                    double ion_part = 0.0;
                    if (is_ion_species(cppargs, i)) {
                        double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                        ion_part = (1.0 - 1.0/eps_born)*cppargs.z[i]*cppargs.z[i]/d_born_i;
                    }
                    // born_diff_mode=2 follows Eq.133-style: remove the born_sum multiplier on the dielectric term.
                    // born_diff_mode=3 disables dielectric-concentration coupling in Born model 1.
                    double eps_part = 0.0;
                    if (cppargs.born_diff_mode == 2) {
                        eps_part = deps_dx_born[i]/(eps_born*eps_born);
                    }
                    else if (cppargs.born_diff_mode == 3) {
                        eps_part = 0.0;
                    }
                    else {
                        eps_part = born_sum*deps_dx_born[i]/(eps_born*eps_born);
                    }
                    ion_part_vec[i] = ion_part;
                    eps_part_vec[i] = eps_part;
                    dadx_born[i] = -Kborn*(ion_part + eps_part);
                }
            }

            double sum_x_dadx_born = 0.0;
            for (int i = 0; i < ncomp; i++) {
                sum_x_dadx_born += x[i]*dadx_born[i];
            }
            a_born_diag = a_born;
            dadx_born_diag = dadx_born;
            sum_x_dadx_born_diag = sum_x_dadx_born;
            for (int i = 0; i < ncomp; i++) {
                mu_born[i] = a_born + Zborn + dadx_born[i] - sum_x_dadx_born;
            }
        }
        else if (cppargs.born_model == 2) {
            const double eps_r_ion = 8.0;
            const double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
            BornSSMDSData born = born_shell_data_cpp(x, cppargs, t, eps_born, eps_r_ion);
            double a_born = -Kborn*born.sum_bracket;
            vector<double> dadx_born(ncomp, 0.0);
            vector<double> direct_part_vec(ncomp, 0.0);
            vector<double> deps_part_vec(ncomp, 0.0);
            vector<double> ddelta_part_vec(ncomp, 0.0);
            if (cppargs.born_diff_mode == 1) {
                dadx_born = contribution_dadx_fd_cpp(AresContributionKind::BORN, t, rho, x, cppargs, a_born);
            }
            else {
                const double inv_eps2 = 1.0/(eps_born*eps_born);
                const double shell_coeff = 1.0/eps_r_ion - 1.0/eps_born;
                const bool use_deps = (cppargs.mu_born_comp_dep_rel_perm != 0);
                const bool use_shell_chain = (cppargs.mu_born_comp_dep_delta_d != 0);
                const double deps_multiplier = (cppargs.mu_born_include_sum_term != 0) ? born.sum_gap : 1.0;
                for (int k = 0; k < ncomp; k++) {
                    double direct_part = 0.0;
                    if (std::abs(cppargs.z[k]) > 1e-12) {
                        direct_part = cppargs.z[k]*cppargs.z[k]*born.bracket[k];
                    }
                    double deps_part = use_deps ? deps_multiplier*deps_dx_born[k]*inv_eps2 : 0.0;
                    double ddelta_part = use_shell_chain ? shell_coeff*born.sum_dpref_over_D2*born.f_k[k] : 0.0;
                    direct_part_vec[k] = direct_part;
                    deps_part_vec[k] = deps_part;
                    ddelta_part_vec[k] = ddelta_part;
                    dadx_born[k] = -Kborn*(direct_part + deps_part + ddelta_part);
                }
            }

            double sum_x_dadx_born = 0.0;
            for (int i = 0; i < ncomp; i++) {
                sum_x_dadx_born += x[i]*dadx_born[i];
            }
            a_born_diag = a_born;
            dadx_born_diag = dadx_born;
            sum_x_dadx_born_diag = sum_x_dadx_born;
            for (int i = 0; i < ncomp; i++) {
                mu_born[i] = a_born + Zborn + dadx_born[i] - sum_x_dadx_born;
            }
        }
        else if (cppargs.born_model != 0) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }
    vector<double> mu(ncomp, 0);
    for (int i = 0; i < ncomp; i++) {
        mu[i] = mu_hc[i] + mu_disp[i] + mu_polar[i] + mu_assoc[i] + mu_ion[i] + mu_born[i];
    }
    CompositionContributionResult result;
    result.dadx = make_vector_terms(dahc_dx, dadisp_dx, dapolar_dx, daassoc_dx, dadx_ion, dadx_born_diag,
        vector<double>());
    result.ares = make_scalar_terms(ares_hc, ares_disp, ares_polar, ares_assoc, a_ion, a_born_diag,
        ares_hc + ares_disp + ares_polar + ares_assoc + a_ion + a_born_diag);
    result.sum_x_dadx = make_scalar_terms(sum_x_dahc_dx, sum_x_dadisp_dx, sum_x_dapolar_dx, sum_x_daassoc_dx,
        sum_x_dadx_ion, sum_x_dadx_born_diag,
        sum_x_dahc_dx + sum_x_dadisp_dx + sum_x_dapolar_dx + sum_x_daassoc_dx + sum_x_dadx_ion + sum_x_dadx_born_diag);
    result.z_raw = z_raw_terms;
    result.z = z_terms;
    vector<double> dadx_total(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        dadx_total[i] = dahc_dx[i] + dadisp_dx[i] + dapolar_dx[i] + daassoc_dx[i] + dadx_ion[i] + dadx_born_diag[i];
    }
    result.dadx.total = dadx_total;
    return result;
}
