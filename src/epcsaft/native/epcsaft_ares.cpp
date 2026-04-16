#include "epcsaft_core_internal.h"
#include "contributions/epcsaft_contrib_internal.h"

using namespace thermo_detail;

namespace {

// EqID: ares_hs
double ares_hs_cpp(const HardChainState &hc_state) {
    const auto &zeta = hc_state.zeta;
    return 1.0 / zeta[0] * (
        3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3])
        + std::pow(zeta[2], 3.0) / (zeta[3] * std::pow(1.0 - zeta[3], 2.0))
        + (std::pow(zeta[2], 3.0) / std::pow(zeta[3], 2.0) - zeta[0]) * std::log(1.0 - zeta[3])
    );
}

// EqID: ares_hc
double ares_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    double summ = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        summ += x[i] * (cppargs.m[i] - 1.0) * std::log(hc_state.ghs[i * ncomp + i]);
    }
    return thermo.m_avg * ares_hs_cpp(hc_state) - summ;
}

// EqID: ares_disp
double ares_disp_cpp(const MixtureState &thermo, const DispersionPolynomialState &dispersion) {
    return -2.0 * PI * thermo.den * dispersion.I1 * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * dispersion.C1 * dispersion.I2 * thermo.m2e2s3;
}

// EqID: ares_assoc
double ares_assoc_cpp(const AssociationIntermediateState &assoc_state, const vector<double> &x) {
    if (!assoc_state.active) {
        return 0.0;
    }
    double value = 0.0;
    for (int i = 0; i < static_cast<int>(assoc_state.setup.site_component_index.size()); ++i) {
        int component_index = assoc_state.setup.site_component_index[i];
        value += x[component_index] * (std::log(assoc_state.XA[i]) - 0.5 * assoc_state.XA[i] + 0.5);
    }
    return value;
}

// EqID: ares_dh
double ares_ion_cpp(double t, const IonIntermediateState &ion_state) {
    if (!ion_state.active) {
        return 0.0;
    }
    double K0 = E_CHRG * E_CHRG / (12.0 * PI * kb * t * perm_vac);
    return -K0 * ion_state.kappa / ion_state.dielectric.eps * ion_state.chi_sum;
}

// EqID: ares_born
double ares_born_cpp(double t, const BornIntermediateState &born_state) {
    if (born_state.model == 0) {
        return 0.0;
    }
    if (born_state.model == 1) {
        return -E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac) * (1.0 - 1.0 / born_state.eps_value) * born_state.charge_radius_sum;
    }
    if (born_state.model == 2) {
        const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        return -Kborn * born_state.shell.sum_bracket;
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
    MixtureState thermo = mixture_state_cpp(t, rho, x, cppargs, false);
    HardChainState hc_state = hard_chain_state_cpp(thermo, x, cppargs);
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(thermo.m_avg, hc_state.eta);
    AssociationIntermediateState assoc_state = association_intermediate_state_cpp(thermo, hc_state, t, x, cppargs, false, false);
    IonIntermediateState ion_state = ion_intermediate_state_cpp(thermo, t, x, cppargs, false);
    BornIntermediateState born_state = born_intermediate_state_cpp(t, x, cppargs, false, false);

    out.hc = ares_hc_cpp(thermo, hc_state, x, cppargs);
    out.disp = ares_disp_cpp(thermo, dispersion);
    out.assoc = ares_assoc_cpp(assoc_state, x);
    out.ion = ares_ion_cpp(t, ion_state);
    out.born = ares_born_cpp(t, born_state);
    return out;
}

ScalarContributionTerms residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    AresContributions contributions = ares_contributions_cpp(t, rho, x, cppargs);
    double ares = contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born;
    return make_scalar_terms(
        contributions.hc,
        contributions.disp,
        contributions.assoc,
        contributions.ion,
        contributions.born,
        ares
    );
}

double ares_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    AresContributions contributions = ares_contributions_cpp(t, rho, x, cppargs);
    return contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born;
}
