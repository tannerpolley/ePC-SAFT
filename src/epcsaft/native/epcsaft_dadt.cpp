#include "epcsaft_core_internal.h"
#include "contributions/epcsaft_contrib_internal.h"

using namespace thermo_detail;

namespace {

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

// EqID: dares_hs_dT
double dadt_hs_cpp(const HardChainState &hc_state, const vector<double> &dzeta_dt) {
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

// EqID: dg_hs_dT
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

// EqID: dares_hc_dT
double dadt_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &dzeta_dt, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    vector<double> contact_time_terms = hc_contact_time_terms_cpp(thermo, hc_state, dzeta_dt);
    double summ = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        summ += x[i] * (cppargs.m[i] - 1.0) * contact_time_terms[i * ncomp + i] / hc_state.ghs[i * ncomp + i];
    }
    return thermo.m_avg * dadt_hs_cpp(hc_state, dzeta_dt) - summ;
}

// EqID: dares_disp_dT
double dadt_disp_cpp(const MixtureState &thermo, double deta_dt, double t, const DispersionPolynomialState &dispersion) {
    double dI1_dt = dispersion.dI1_deta * deta_dt;
    double dI2_dt = dispersion.dI2_deta * deta_dt;
    double dC1_dt = dispersion.C2 * deta_dt;
    return -2.0 * PI * thermo.den * (dI1_dt - dispersion.I1 / t) * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * (dC1_dt * dispersion.I2 + dispersion.C1 * dI2_dt - 2.0 * dispersion.C1 * dispersion.I2 / t) * thermo.m2e2s3;
}

// EqID: dares_assoc_dT
double dadt_assoc_cpp(const AssociationIntermediateState &assoc_state, const vector<double> &x) {
    if (!assoc_state.active) {
        return 0.0;
    }
    double value = 0.0;
    for (int i = 0; i < static_cast<int>(assoc_state.setup.site_component_index.size()); ++i) {
        value += x[assoc_state.setup.site_component_index[i]] * (1.0 / assoc_state.XA[i] - 0.5) * assoc_state.dXA_dt[i];
    }
    return value;
}

// EqID: dares_dh_dT
double dadt_ion_cpp(const IonIntermediateState &ion_state, double t, const vector<double> &x, const add_args &cppargs) {
    if (!ion_state.active) {
        return 0.0;
    }
    vector<double> q(cppargs.z.begin(), cppargs.z.end());
    for (double &qi : q) {
        qi *= E_CHRG;
    }
    double dkappa_dt = -0.5 * ion_state.kappa / t;
    double dadt_ion = 0.0;
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        dadt_ion += x[i] * q[i] * q[i] * (ion_state.sigma_k[i] * dkappa_dt / t - ion_state.chi[i] * ion_state.kappa / (t * t));
    }
    return -1.0 / 12.0 / PI / kb / (ion_state.dielectric.eps * perm_vac) * dadt_ion;
}

// EqID: dares_born_dT
double dadt_born_cpp(double t, const BornIntermediateState &born_state) {
    if (born_state.model == 0) {
        return 0.0;
    }
    if (born_state.model == 1) {
        double born_factor = 1.0 - 1.0 / born_state.eps_value;
        double prefactor = E_CHRG * E_CHRG / (4.0 * PI * kb * perm_vac);
        return prefactor * born_factor * (born_state.charge_radius_sum / (t * t) - born_state.charge_radius_sum_dt / t);
    }
    if (born_state.model == 2) {
        return E_CHRG * E_CHRG / (4.0 * PI * kb * perm_vac * t * t) * born_state.shell.sum_bracket;
    }
    throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
}

}  // namespace

// EqID: dares_dT
ScalarContributionTerms temperature_derivative_residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    MixtureState thermo = mixture_state_cpp(t, rho, x, cppargs, true);
    HardChainState hc_state = hard_chain_state_cpp(thermo, x, cppargs);
    vector<double> dzeta_dt = dzeta_dt_cpp(thermo, x, cppargs);
    vector<double> dghs_dt = hc_contact_time_terms_cpp(thermo, hc_state, dzeta_dt);
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(thermo.m_avg, hc_state.eta);
    AssociationIntermediateState assoc_state = association_intermediate_state_cpp(thermo, hc_state, t, x, cppargs, true, false, &dghs_dt);
    IonIntermediateState ion_state = ion_intermediate_state_cpp(thermo, t, x, cppargs, false);
    BornIntermediateState born_state = born_intermediate_state_cpp(t, x, cppargs, true, false);

    double hc = dadt_hc_cpp(thermo, hc_state, dzeta_dt, x, cppargs);
    double disp = dadt_disp_cpp(thermo, dzeta_dt[3], t, dispersion);
    double assoc = dadt_assoc_cpp(assoc_state, x);
    double ion = dadt_ion_cpp(ion_state, t, x, cppargs);
    double born = dadt_born_cpp(t, born_state);
    double total = hc + disp + assoc + ion + born;
    return make_scalar_terms(hc, disp, assoc, ion, born, total);
}

double dadt_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return temperature_derivative_residual_helmholtz_result_cpp(t, rho, std::move(x), cppargs).total;
}
