#include "epcsaft_contrib_internal.h"

using namespace thermo_detail;

double ion_diameter_cpp(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return cppargs.s[i];
    }
    int mode = cppargs.d_ion_mode;
    double sigma_i = cppargs.s[i];
    if (sigma_i <= 0.0) {
        throw ValueError("DH/ion diameter requires positive ionic sigma_i.");
    }
    if (mode == 0) {
        return sigma_i;
    }
    if (mode == 1) {
        return sigma_i * (1.0 - 0.12);
    }
    if (mode == 2) {
        return sigma_i * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
    }
    throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
}

double ion_diameter_cpp_dt(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return 0.0;
    }
    if (cppargs.d_ion_mode == 2) {
        double sigma_i = cppargs.s[i];
        double expo = std::exp(-3.0 * cppargs.e[i] / t);
        return -0.36 * sigma_i * cppargs.e[i] * expo / (t * t);
    }
    return 0.0;
}

double dh_kappa_cpp(double den, double t, double eps, double q2_sum) {
    return std::sqrt(den * E_CHRG * E_CHRG / kb / t / (eps * perm_vac) * q2_sum);
}

double dh_chi_cpp(double kappa, double diameter) {
    double ka = kappa * diameter;
    return 3.0 / std::pow(ka, 3.0) * (1.5 + std::log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * std::pow(1.0 + ka, 2.0));
}

IonIntermediateState ion_intermediate_state_cpp(
    const MixtureState &thermo,
    double t,
    const vector<double> &x,
    const add_args &cppargs,
    bool include_dx
) {
    IonIntermediateState state;
    if (cppargs.z.empty()) {
        return state;
    }

    state.dielectric = dielectric_state_cpp(x, cppargs);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        state.charge_square_sum += cppargs.z[i] * cppargs.z[i] * x[i];
    }
    if (state.charge_square_sum == 0.0) {
        return state;
    }

    state.kappa = dh_kappa_cpp(thermo.den, t, state.dielectric.eps, state.charge_square_sum);
    if (state.kappa == 0.0) {
        return state;
    }

    state.active = true;
    int ncomp = static_cast<int>(x.size());
    state.chi.assign(ncomp, 0.0);
    state.sigma_k.assign(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        state.chi[i] = dh_chi_cpp(state.kappa, thermo.d[i]);
        state.sigma_k[i] = -2.0 * state.chi[i] + 3.0 / (1.0 + state.kappa * thermo.d[i]);
        state.chi_sum += x[i] * cppargs.z[i] * cppargs.z[i] * state.chi[i];
        state.sigma_sum += x[i] * cppargs.z[i] * cppargs.z[i] * state.sigma_k[i];
    }

    if (include_dx) {
        state.dkappa_dx.assign(ncomp, 0.0);
        state.dchi_sum_dx.assign(ncomp, 0.0);
        const bool use_dh_deps = (cppargs.mu_DH_comp_dep_rel_perm != 0);
        const double dh_deps_multiplier = (cppargs.mu_DH_include_sum_term != 0) ? state.charge_square_sum : 1.0;
        double a_const = thermo.den * E_CHRG * E_CHRG / (kb * t * perm_vac);
        for (int i = 0; i < ncomp; ++i) {
            double deps_term = use_dh_deps ? dh_deps_multiplier * state.dielectric.deps_dx[i] / (state.dielectric.eps * state.dielectric.eps) : 0.0;
            state.dkappa_dx[i] = a_const * (cppargs.z[i] * cppargs.z[i] / state.dielectric.eps - deps_term) / (2.0 * state.kappa);
        }
        for (int i = 0; i < ncomp; ++i) {
            state.dchi_sum_dx[i] = cppargs.z[i] * cppargs.z[i] * state.chi[i] + state.dkappa_dx[i] * (state.sigma_sum - state.chi_sum) / state.kappa;
        }
    }

    return state;
}

// EqID: dadrho_dh_explicit
double dadrho_ion_cpp(double t, const IonIntermediateState &ion_state) {
    if (!ion_state.active) {
        return 0.0;
    }
    return -ion_state.kappa / 24.0 / PI / kb / t / (ion_state.dielectric.eps * perm_vac) * ion_state.sigma_sum * E_CHRG * E_CHRG;
}

// EqID: dh_ares_dT
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

// EqID: dh_ares_dxi
ContributionDadxResult dadx_ion_cpp(const MixtureState &thermo, const IonIntermediateState &ion_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);
    if (cppargs.z.empty()) {
        return result;
    }

    int dh_model = cppargs.DH_model;
    if (dh_model == 2) {
        throw ValueError("DH_model=2 (Bjerrum treatment) is reserved and not implemented.");
    }
    if ((dh_model != 0) && (dh_model != 1)) {
        throw ValueError("Unknown DH_model. Supported values are 0, 1, and reserved 2.");
    }
    if (!ion_state.active) {
        return result;
    }

    double K0 = E_CHRG * E_CHRG / (12.0 * PI * kb * t * perm_vac);
    result.ares = -K0 * ion_state.kappa / ion_state.dielectric.eps * ion_state.chi_sum;
    result.z = -ion_state.kappa / 24.0 / PI / kb / t / (ion_state.dielectric.eps * perm_vac) * ion_state.sigma_sum * E_CHRG * E_CHRG;

    if (cppargs.mu_DH_diff_mode == 1) {
        result.dadx = contribution_dadx_fd_cpp(AresContributionKind::ION, t, rho, x, cppargs, result.ares);
    } else if (cppargs.mu_DH_diff_mode == 2) {
        result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::ION, t, rho, x, cppargs);
    } else {
        const bool use_dh_deps = (cppargs.mu_DH_comp_dep_rel_perm != 0);
        for (int i = 0; i < ncomp; ++i) {
            double d_inv_eps_dx = use_dh_deps ? -ion_state.dielectric.deps_dx[i] / (ion_state.dielectric.eps * ion_state.dielectric.eps) : 0.0;
            double term1 = (ion_state.dkappa_dx[i] / ion_state.dielectric.eps + ion_state.kappa * d_inv_eps_dx) * ion_state.chi_sum;
            double term2 = ion_state.kappa / ion_state.dielectric.eps * ion_state.dchi_sum_dx[i];
            result.dadx[i] = -K0 * (term1 + term2);
        }
    }

    for (int i = 0; i < ncomp; ++i) {
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    return result;
}
