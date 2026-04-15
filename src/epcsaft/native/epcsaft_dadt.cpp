#include "epcsaft_core_internal.h"

using namespace thermo_detail;

// EqID: dares_dT
ScalarContributionTerms temperature_derivative_residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size()); // number of components
    double summ = 0.0;
    ThermoCommonState thermo = thermo_common_state_cpp(t, rho, x, cppargs, true);
    auto &d = thermo.d;
    auto &dd_dt = thermo.dd_dt;
    auto &zeta = thermo.zeta;
    auto &dzeta_dt = thermo.dzeta_dt;
    auto &e_ij = thermo.e_ij;
    auto &s_ij = thermo.s_ij;
    auto &ghs = thermo.ghs;
    auto &denghs = thermo.denghs;
    double den = thermo.den;
    double eta = thermo.eta;
    double m_avg = thermo.m_avg;
    double m2es3 = thermo.m2es3;
    double m2e2s3 = thermo.m2e2s3;
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(m_avg, eta);
    vector<double> dghs_dt(ncomp * ncomp, 0.0);
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            double pair_diameter = pair_diameter_cpp(d[i], d[j]);
            double pair_diameter_dt = pair_diameter * (
                dd_dt[i] / d[i] + dd_dt[j] / d[j] - (dd_dt[i] + dd_dt[j]) / (d[i] + d[j])
            );
            dghs_dt[idx] = hs_contact_time_derivative_cpp(
                pair_diameter,
                pair_diameter_dt,
                zeta[2],
                zeta[3],
                dzeta_dt[2],
                dzeta_dt[3]
            );
        }
    }

    double dadt_hs = 1/zeta[0]*(3*(dzeta_dt[1]*zeta[2] + zeta[1]*dzeta_dt[2])/(1-zeta[3])
        + 3*zeta[1]*zeta[2]*dzeta_dt[3]/pow(1-zeta[3], 2.)
        + 3*pow(zeta[2], 2.)*dzeta_dt[2]/zeta[3]/pow(1-zeta[3], 2.)
        + pow(zeta[2],3.)*dzeta_dt[3]*(3*zeta[3]-1)/pow(zeta[3], 2.)/pow(1-zeta[3], 3.)
        + (3*pow(zeta[2], 2.)*dzeta_dt[2]*zeta[3] - 2*pow(zeta[2], 3.)*dzeta_dt[3])/pow(zeta[3], 3.)
        * log(1-zeta[3])
        + (zeta[0]-pow(zeta[2],3)/pow(zeta[3],2.))*dzeta_dt[3]/(1-zeta[3]));
    double I1 = dispersion.I1;
    double I2 = dispersion.I2;
    double dI1_dt = dispersion.dI1_deta * dzeta_dt[3];
    double dI2_dt = dispersion.dI2_deta * dzeta_dt[3];
    double C1 = dispersion.C1;
    double C2 = dispersion.C2;
    double dC1_dt = C2 * dzeta_dt[3];

    summ = 0.;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i]*(cppargs.m[i]-1)*dghs_dt[i*ncomp+i]/ghs[i*ncomp+i];
    }

    double dadt_hc = m_avg*dadt_hs - summ;
    double dadt_disp = -2*PI*den*(dI1_dt-I1/t)*m2es3 - PI*den*m_avg*(dC1_dt*I2+C1*dI2_dt-2*C1*I2/t)*m2e2s3;

    // Dipole term (Gross and Vrabec term) --------------------------------------
    double dadt_polar = 0.;
    if (!cppargs.dipm.empty()) {
        double A2 = 0.;
        double A3 = 0.;
        double dA2_dt = 0.;
        double dA3_dt = 0.;
        vector<double> dipmSQ (ncomp, 0);

        const static double conv = 7242.702976750923; // conversion factor, see the note below Table 2 in Gross and Vrabec 2006

        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(cppargs.dipm[i], 2.)/(cppargs.m[i]*cppargs.e[i]*pow(cppargs.s[i],3.))*conv;
        }

        double J2, J3, dJ2_dt, dJ3_dt;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                double m_ij = sqrt(cppargs.m[i]*cppargs.m[j]);
                if (m_ij > 2) {
                    m_ij = 2;
                }
                vector<double> adip = dipole_coefficients_cpp(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = dipole_coefficients_cpp(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
                J2 = 0.;
                dJ2_dt = 0.;
                for (int l = 0; l < 5; l++) {
                    J2 += (adip[l] + bdip[l]*e_ij[j*ncomp+j]/t)*pow(eta, l); // j*ncomp+j needs to be used for e_ij because it is formatted as a 1D vector
                    dJ2_dt += adip[l]*l*pow(eta, l-1)*dzeta_dt[3]
                        + bdip[l]*e_ij[j*ncomp+j]*(1/t*l*pow(eta, l-1)*dzeta_dt[3]
                        - 1/pow(t,2.)*pow(eta,l));
                }
                A2 += x[i]*x[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                    pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*J2;
                dA2_dt += x[i]*x[j]*e_ij[i*ncomp+i]*e_ij[j*ncomp+j]*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)
                    /pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*
                    (dJ2_dt/pow(t,2)-2*J2/pow(t,3));

                for (int k = 0; k < ncomp; k++) {
                    double m_ijk = pow((cppargs.m[i]*cppargs.m[j]*cppargs.m[k]),1/3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    vector<double> cdip = dipole_coefficients_cpp(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
                    J3 = 0.;
                    dJ3_dt = 0.;
                    for (int l = 0; l < 5; l++) {
                        J3 += cdip[l]*pow(eta, l);
                        dJ3_dt += cdip[l]*l*pow(eta, l-1)*dzeta_dt[3];
                    }
                    A3 += x[i]*x[j]*x[k]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*e_ij[k*ncomp+k]/t*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]/
                        s_ij[j*ncomp+k]*cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]*
                        dipmSQ[j]*dipmSQ[k]*J3;
                    dA3_dt += x[i]*x[j]*x[k]*e_ij[i*ncomp+i]*e_ij[j*ncomp+j]*e_ij[k*ncomp+k]*
                        pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)*pow(s_ij[k*ncomp+k],3)/s_ij[i*ncomp+j]/s_ij[i*ncomp+k]
                        /s_ij[j*ncomp+k]*cppargs.dip_num[i]*cppargs.dip_num[j]*cppargs.dip_num[k]*dipmSQ[i]
                        *dipmSQ[j]*dipmSQ[k]*(-3*J3/pow(t,4) + dJ3_dt/pow(t,3));
                }
            }
        }

        A2 = -PI*den*A2;
        A3 = -4/3.*PI*PI*den*den*A3;
        dA2_dt = -PI*den*dA2_dt;
        dA3_dt = -4/3.*PI*PI*den*den*dA3_dt;

        if (A2 != 0) { // when the mole fraction of the polar compounds is 0 then A2 = 0 and division by 0 occurs
            dadt_polar = (dA2_dt-2*A3/A2*dA2_dt+dA3_dt)/pow(1-A3/A2, 2.);
        }
    }

    // Association term -------------------------------------------------------
    // only the 2B association type is currently implemented
    double dadt_assoc = 0.;
    if (!cppargs.e_assoc.empty()) {
        AssociationSetup assoc = association_setup_cpp(x, cppargs, s_ij, ghs, t);
        const vector<int> &iA = assoc.site_component_index;
        const vector<double> &x_assoc = assoc.x_assoc;
        const vector<double> &delta_ij = assoc.delta_ij;
        int num_sites = static_cast<int>(iA.size());

        vector<double> XA(num_sites, 0.0);
        vector<double> ddelta_dt(num_sites * num_sites, 0.0);
        for (int i = 0; i < num_sites; i++) {
            int idxi = iA[i]*ncomp+iA[i];
            for (int j = 0; j < num_sites; j++) {
                int idxj = iA[j]*ncomp+iA[j];
                if (cppargs.assoc_matrix[i*num_sites+j] != 0) {
                    double eABij = (cppargs.e_assoc[iA[i]]+cppargs.e_assoc[iA[j]])/2.;
                    double volABij = _HUGE;
                    if (cppargs.k_hb.empty()) {
                        volABij = sqrt(cppargs.vol_a[iA[i]]*cppargs.vol_a[iA[j]])*pow(sqrt(s_ij[idxi]*
                            s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3);
                    }
                    else {
                        volABij = sqrt(cppargs.vol_a[iA[i]]*cppargs.vol_a[iA[j]])*pow(sqrt(s_ij[idxi]*
                            s_ij[idxj])/(0.5*(s_ij[idxi]+s_ij[idxj])), 3)*(1-cppargs.k_hb[iA[i]*ncomp+iA[j]]);
                    }
                    ddelta_dt[i*num_sites+j] = pow(s_ij[idxj],3)*volABij*(-eABij/pow(t,2)
                        *exp(eABij/t)*ghs[iA[i]*ncomp+iA[j]] + dghs_dt[iA[i]*ncomp+iA[j]]
                        *(exp(eABij/t)-1));
                }
            }
        }
        XA = solve_association_site_fractions_cpp(delta_ij, den, x_assoc);

        vector<double> dXA_dt(num_sites, 0);
        dXA_dt = ::association_site_fraction_dt_cpp(delta_ij, den, XA, ddelta_dt, x_assoc);

        for (int i = 0; i < num_sites; i++) {
            dadt_assoc += x[iA[i]]*(1/XA[i]-0.5)*dXA_dt[i];
        }
    }

    // Ion term ---------------------------------------------------------------
    double dadt_ion = 0.;
    double dadt_born = 0.;
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

        double dkappa_dt;
        if (kappa != 0) {
            vector<double> chi(ncomp);
            vector<double> dchikap_dk(ncomp);
            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi[i] = dh_chi_cpp(kappa, d[i]);
                dchikap_dk[i] = -2*chi[i]+3/(1+kappa*d[i]);
                summ += x[i]*cppargs.z[i]*cppargs.z[i];
            }
            dkappa_dt = -0.5*den*E_CHRG*E_CHRG/kb/t/t/(eps*perm_vac)*summ/kappa;

            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                summ += x[i]*q[i]*q[i]*(dchikap_dk[i]*dkappa_dt/t-kappa*chi[i]/t/t);
            }
            dadt_ion = -1/12./PI/kb/(eps*perm_vac)*summ;
        }

        if (cppargs.born_model == 1) {
            // Born term temperature derivative (non-SSM+DS) with user-selected born_radius_model
            double born_sum = 0.;
            double born_sum_dt = 0.;
            for (int i = 0; i < ncomp; i++) {
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                    double d_born_dt = ion_born_radius_cpp_dt(i, t, cppargs);
                    born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
                    born_sum_dt += x[i]*cppargs.z[i]*cppargs.z[i]*(-d_born_dt)/(d_born_i*d_born_i);
                }
            }
            double born_factor = (1.-1./eps_born);
            double prefactor = E_CHRG*E_CHRG/(4.*PI*kb*perm_vac);
            dadt_born = prefactor*born_factor*(born_sum/(t*t) - born_sum_dt/t);
        }
        else if (cppargs.born_model == 2) {
            const double eps_r_ion = 8.0;
            BornSSMDSData born = born_shell_data_cpp(x, cppargs, t, eps_born, eps_r_ion);
            dadt_born = E_CHRG*E_CHRG/(4.*PI*kb*perm_vac*t*t)*born.sum_bracket;
        }
        else if (cppargs.born_model != 0) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }

    double dadt = dadt_hc + dadt_disp + dadt_assoc + dadt_polar + dadt_ion + dadt_born;
    return make_scalar_terms(dadt_hc, dadt_disp, dadt_polar, dadt_assoc, dadt_ion, dadt_born, dadt);
}

double dadt_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return temperature_derivative_residual_helmholtz_result_cpp(t, rho, std::move(x), cppargs).total;
}

// EqID: h_res
double hres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double Z = Z_cpp(t, rho, x, cppargs);
    double dares_dt = dadt_cpp(t, rho, std::move(x), cppargs);
    return (-t * dares_dt + (Z - 1.0)) * kb * N_AV * t;
}

// EqID: g_res_from_ares
double gres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double ares = ares_cpp(t, rho, x, cppargs);
    double Z = Z_cpp(t, rho, std::move(x), cppargs);
    return (ares + (Z - 1.0) - std::log(Z)) * kb * N_AV * t;
}

// EqID: s_res
double sres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double gres = gres_cpp(t, rho, x, cppargs);
    double hres = hres_cpp(t, rho, std::move(x), cppargs);
    return (hres - gres) / t;
}
