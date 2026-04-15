#include "epcsaft_core_internal.h"

using namespace thermo_detail;

// EqID: z_eta_rho_identity
DadrhoResult dadrho_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ThermoCommonState thermo = thermo_common_state_cpp(t, rho, x, cppargs, false);
    auto &d = thermo.d;
    auto &zeta = thermo.zeta;
    auto &e_ij = thermo.e_ij;
    auto &s_ij = thermo.s_ij;
    auto &ghs = thermo.ghs;
    auto &denghs = thermo.denghs;
    double den = thermo.den;
    double eta = thermo.eta;
    double m_avg = thermo.m_avg;
    double m2es3 = thermo.m2es3;
    double m2e2s3 = thermo.m2e2s3;
    double summ = 0.0;
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(m_avg, eta);
    double detI1_det = dispersion.dEtaI1_deta;
    double detI2_det = dispersion.dEtaI2_deta;
    double I2 = dispersion.I2;
    double C1 = dispersion.C1;
    double C2 = dispersion.C2;
    double Zhs = zeta[3]/(1-zeta[3]) + 3.*zeta[1]*zeta[2]/zeta[0]/(1.-zeta[3])/(1.-zeta[3]) +
        (3.*pow(zeta[2], 3.) - zeta[3]*pow(zeta[2], 3.))/zeta[0]/pow(1.-zeta[3], 3.);

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i]*(cppargs.m[i]-1)/ghs[i*ncomp+i]*denghs[i*ncomp+i];
    }

    double Zhc = m_avg*Zhs - summ;
    double Zdisp = -2*PI*den*detI1_det*m2es3 - PI*den*m_avg*(C1*detI2_det + C2*eta*I2)*m2e2s3;

    // Dipole term (Gross and Vrabec term) --------------------------------------
    double Zpolar = 0;
    if (!cppargs.dipm.empty()) {
        double A2 = 0.;
        double A3 = 0.;
        double dA2_det = 0.;
        double dA3_det = 0.;
        vector<double> adip (5, 0);
        vector<double> bdip (5, 0);
        vector<double> cdip (5, 0);
        vector<double> dipmSQ (ncomp, 0);
        double J2, detJ2_det, J3, detJ3_det;

        const static double conv = 7242.702976750923; // conversion factor, see the note below Table 2 in Gross and Vrabec 2006

        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(cppargs.dipm[i], 2.)/(cppargs.m[i]*cppargs.e[i]*pow(cppargs.s[i],3.))*conv;
        }

        double m_ij;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                m_ij = sqrt(cppargs.m[i]*cppargs.m[j]);
                if (m_ij > 2) {
                    m_ij = 2;
                }
                vector<double> adip = dipole_coefficients_cpp(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = dipole_coefficients_cpp(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
                J2 = 0.;
                detJ2_det = 0.;
                for (int l = 0; l < 5; l++) {
                    J2 += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/t)*pow(eta, l); // i*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    detJ2_det += (adip[l] + bdip[l]*e_ij[i*ncomp+j]/t)*(l+1)*pow(eta, l);
                }
                A2 += x[i]*x[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                    pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*J2;
                dA2_det += x[i]*x[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*
                    pow(s_ij[j*ncomp+j],3)/pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*detJ2_det;
            }
        }

        double m_ijk;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                for (int k = 0; k < ncomp; k++) {
                    m_ijk = pow((cppargs.m[i]*cppargs.m[j]*cppargs.m[k]),1/3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    vector<double> cdip = dipole_coefficients_cpp(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
                    J3 = 0.;
                    detJ3_det = 0.;
                    for (int l = 0; l < 5; l++) {
                        J3 += cdip[l]*pow(eta, l);
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
                }
            }
        }

        A2 = -PI*den*A2;
        A3 = -4/3.*PI*PI*den*den*A3;
        dA2_det = -PI*den/eta*dA2_det;
        dA3_det = -4/3.*PI*PI*den/eta*den/eta*dA3_det;

        if (A2 != 0) { // when the mole fraction of the polar compounds is 0 then A2 = 0 and division by 0 occurs
            Zpolar = eta*((dA2_det*(1-A3/A2)+(dA3_det*A2-A3*dA2_det)/A2)/(1-A3/A2)/(1-A3/A2));
        }
    }

    // Association term -------------------------------------------------------
    double Zassoc = 0;
    if (!cppargs.e_assoc.empty()) {
        AssociationSetup assoc = association_setup_cpp(x, cppargs, s_ij, ghs, t);
        const vector<int> &iA = assoc.site_component_index;
        const vector<double> &x_assoc = assoc.x_assoc;
        const vector<double> &delta_ij = assoc.delta_ij;
        int num_sites = static_cast<int>(iA.size());

        vector<double> XA(num_sites, 0.0);

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

        XA = solve_association_site_fractions_cpp(delta_ij, den, x_assoc);

        vector<double> dXA_dx(num_sites*ncomp, 0);
        dXA_dx = ::association_site_fraction_dx_cpp(cppargs.assoc_num, delta_ij, den, XA, ddelta_dx, x_assoc);

        summ = 0.;
        int ij = 0;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < num_sites; j++) {
                summ += x[i]*den*x[iA[j]]*(1/XA[j]-0.5)*dXA_dx[ij];
                ij += 1;
            }
        }

        Zassoc = summ;
    }

    // Ion term ---------------------------------------------------------------
    double Zion = 0;
    double Zborn = 0;
    if (!cppargs.z.empty()) {
        DielectricState dielc_state = dielectric_state_cpp(x, cppargs);
        double eps = dielc_state.eps;
        double eps_born = (cppargs.born_eps_mode == 1) ? reference_solvent_dielectric_constant_cpp(x, cppargs) : eps;
        vector<double> q(cppargs.z.begin(), cppargs.z.end());
        for (int i = 0; i < ncomp; i++) {
            q[i] = q[i]*E_CHRG;
        }

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            summ += cppargs.z[i]*cppargs.z[i]*x[i];
        }

        double kappa = dh_kappa_cpp(den, t, eps, summ); // the inverse Debye screening length. Equation 4 in Held et al. 2008.

        if (kappa != 0) {
            double chi, sigma_k;
            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi = dh_chi_cpp(kappa, d[i]);
                sigma_k = -2*chi+3/(1+kappa*d[i]);
                summ += q[i]*q[i]*x[i]*sigma_k;
            }
            Zion = -1*kappa/24./PI/kb/t/(eps*perm_vac)*summ;
        }
    }

    // Born term (Bulow 2021a, non-SSM+DS) has no explicit density dependence.
    double z_increment = Zhc + Zdisp + Zpolar + Zassoc + Zion + Zborn;
    ScalarContributionTerms raw_terms = make_scalar_terms(Zhc, Zdisp, Zpolar, Zassoc, Zion, Zborn, z_increment);
    DadrhoResult result;
    result.raw = raw_terms;
    result.terms = normalized_dadrho_terms_cpp(raw_terms);
    return result;
}
