#include "epcsaft_core_internal.h"
#include "contributions/epcsaft_contrib_internal.h"

#ifdef EPCSAFT_HAS_CPPAD
#include <cppad/cppad.hpp>
#endif

#include <map>
#include <string>

using thermo_detail::AresContributionKind;
using thermo_detail::AresContributions;
using thermo_detail::AssociationIntermediateState;
using thermo_detail::BornIntermediateState;
using thermo_detail::DadrhoResult;
using thermo_detail::DispersionPolynomialState;
using thermo_detail::HardChainState;
using thermo_detail::IonIntermediateState;
using thermo_detail::MixtureState;

namespace ares_detail {

template <typename Scalar>
struct MixtureStateScalar {
    vector<double> d;
    vector<double> e_ij;
    vector<double> s_ij;
    Scalar den = scalar_constant<Scalar>(0.0);
    Scalar m_avg = scalar_constant<Scalar>(0.0);
    Scalar m2es3 = scalar_constant<Scalar>(0.0);
    Scalar m2e2s3 = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
struct HardChainStateScalar {
    vector<Scalar> zeta;
    vector<Scalar> ghs;
    Scalar eta = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
struct DispersionPolynomialStateScalar {
    std::array<Scalar, 7> a{};
    std::array<Scalar, 7> b{};
    Scalar I1 = scalar_constant<Scalar>(0.0);
    Scalar I2 = scalar_constant<Scalar>(0.0);
    Scalar dEtaI1_deta = scalar_constant<Scalar>(0.0);
    Scalar dEtaI2_deta = scalar_constant<Scalar>(0.0);
    Scalar C1 = scalar_constant<Scalar>(0.0);
    Scalar C2 = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
struct AresContributionsScalar {
    Scalar hc = scalar_constant<Scalar>(0.0);
    Scalar disp = scalar_constant<Scalar>(0.0);
    Scalar assoc = scalar_constant<Scalar>(0.0);
    Scalar ion = scalar_constant<Scalar>(0.0);
    Scalar born = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
static MixtureStateScalar<Scalar> mixture_state_scalar_cpp(double t, const Scalar &rho, const vector<Scalar> &x, const add_args &cppargs) {
    MixtureStateScalar<Scalar> state;
    int ncomp = static_cast<int>(x.size());
    state.d.assign(ncomp, 0.0);
    state.e_ij.assign(ncomp * ncomp, 0.0);
    state.s_ij.assign(ncomp * ncomp, 0.0);
    state.den = rho * N_AV / 1.0e30;

    for (int i = 0; i < ncomp; ++i) {
        state.d[i] = cppargs.s[i] * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
        if (!cppargs.z.empty() && std::abs(cppargs.z[i]) > 1e-12) {
            state.d[i] = thermo_detail::parameter_setup_detail::ion_diameter_cpp(i, t, cppargs);
        }
        state.m_avg += x[i] * cppargs.m[i];
    }

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            state.s_ij[idx] = thermo_detail::parameter_setup_detail::pair_sigma_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.e_ij[idx] = thermo_detail::parameter_setup_detail::pair_epsilon_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.m2es3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * state.e_ij[idx] / t * std::pow(state.s_ij[idx], 3);
            state.m2e2s3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * std::pow(state.e_ij[idx] / t, 2) * std::pow(state.s_ij[idx], 3);
        }
    }
    return state;
}

template <typename Scalar>
static Scalar hs_contact_value_scalar_cpp(double pair_diameter, const Scalar &zeta2, const Scalar &zeta3) {
    const Scalar one = scalar_constant<Scalar>(1.0);
    return one / (one - zeta3)
        + pair_diameter * 3.0 * zeta2 / scalar_pow(one - zeta3, 2)
        + scalar_pow(pair_diameter, 2.0) * 2.0 * zeta2 * zeta2 / scalar_pow(one - zeta3, 3);
}

template <typename Scalar>
static HardChainStateScalar<Scalar> hard_chain_state_scalar_cpp(const MixtureStateScalar<Scalar> &thermo, const vector<Scalar> &x, const add_args &cppargs) {
    HardChainStateScalar<Scalar> state;
    int ncomp = static_cast<int>(x.size());
    state.zeta.assign(4, scalar_constant<Scalar>(0.0));
    state.ghs.assign(ncomp * ncomp, scalar_constant<Scalar>(0.0));
    for (int k = 0; k < 4; ++k) {
        Scalar summ = scalar_constant<Scalar>(0.0);
        for (int j = 0; j < ncomp; ++j) {
            summ += x[j] * cppargs.m[j] * scalar_pow(thermo.d[j], k);
        }
        state.zeta[k] = PI / 6.0 * thermo.den * summ;
    }
    state.eta = state.zeta[3];

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            double pair_diameter = thermo_detail::parameter_setup_detail::pair_diameter_cpp(thermo.d[i], thermo.d[j]);
            state.ghs[idx] = hs_contact_value_scalar_cpp(pair_diameter, state.zeta[2], state.zeta[3]);
        }
    }
    return state;
}

template <typename Scalar>
static DispersionPolynomialStateScalar<Scalar> dispersion_polynomials_scalar_cpp(const Scalar &m_avg, const Scalar &eta) {
    DispersionPolynomialStateScalar<Scalar> state;
    Scalar c1 = (m_avg - 1.0) / m_avg;
    Scalar c2 = (m_avg - 2.0) / m_avg;
    for (size_t i = 0; i < state.a.size(); ++i) {
        state.a[i] = thermo_detail::kDispersionA0[i] + c1 * thermo_detail::kDispersionA1[i] + c1 * c2 * thermo_detail::kDispersionA2[i];
        state.b[i] = thermo_detail::kDispersionB0[i] + c1 * thermo_detail::kDispersionB1[i] + c1 * c2 * thermo_detail::kDispersionB2[i];
        state.I1 += state.a[i] * scalar_pow(eta, static_cast<int>(i));
        state.I2 += state.b[i] * scalar_pow(eta, static_cast<int>(i));
        state.dEtaI1_deta += state.a[i] * static_cast<double>(i + 1) * scalar_pow(eta, static_cast<int>(i));
        state.dEtaI2_deta += state.b[i] * static_cast<double>(i + 1) * scalar_pow(eta, static_cast<int>(i));
    }
    const Scalar one = scalar_constant<Scalar>(1.0);
    state.C1 = one / (one
        + m_avg * (8.0 * eta - 2.0 * eta * eta) / scalar_pow(one - eta, 4)
        + (one - m_avg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * scalar_pow(eta, 3) - 2.0 * scalar_pow(eta, 4))
            / scalar_pow((one - eta) * (2.0 - eta), 2));
    state.C2 = -state.C1 * state.C1 * (
        m_avg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / scalar_pow(one - eta, 5)
        + (one - m_avg) * (2.0 * scalar_pow(eta, 3) + 12.0 * eta * eta - 48.0 * eta + 40.0)
            / scalar_pow((one - eta) * (2.0 - eta), 3));
    return state;
}

// EqID: ares_hs
template <typename Scalar>
static Scalar ares_hs_scalar_cpp(const HardChainStateScalar<Scalar> &hc_state) {
    const auto &zeta = hc_state.zeta;
    return 1.0 / zeta[0] * (
        3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3])
        + scalar_pow(zeta[2], 3) / (zeta[3] * scalar_pow(1.0 - zeta[3], 2))
        + (scalar_pow(zeta[2], 3) / scalar_pow(zeta[3], 2) - zeta[0]) * scalar_log(1.0 - zeta[3])
    );
}

// EqID: ares_hc
template <typename Scalar>
static Scalar ares_hc_scalar_cpp(const MixtureStateScalar<Scalar> &thermo, const HardChainStateScalar<Scalar> &hc_state, const vector<Scalar> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    Scalar summ = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < ncomp; ++i) {
        summ += x[i] * (cppargs.m[i] - 1.0) * scalar_log(hc_state.ghs[i * ncomp + i]);
    }
    return thermo.m_avg * ares_hs_scalar_cpp(hc_state) - summ;
}

// EqID: ares_disp
template <typename Scalar>
static Scalar ares_disp_scalar_cpp(const MixtureStateScalar<Scalar> &thermo, const DispersionPolynomialStateScalar<Scalar> &dispersion) {
    return -2.0 * PI * thermo.den * dispersion.I1 * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * dispersion.C1 * dispersion.I2 * thermo.m2e2s3;
}

static double ares_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &x, const add_args &cppargs) {
    MixtureStateScalar<double> scalar_thermo;
    scalar_thermo.d = thermo.d;
    scalar_thermo.e_ij = thermo.e_ij;
    scalar_thermo.s_ij = thermo.s_ij;
    scalar_thermo.den = thermo.den;
    scalar_thermo.m_avg = thermo.m_avg;
    scalar_thermo.m2es3 = thermo.m2es3;
    scalar_thermo.m2e2s3 = thermo.m2e2s3;

    HardChainStateScalar<double> scalar_hc;
    scalar_hc.zeta = hc_state.zeta;
    scalar_hc.ghs = hc_state.ghs;
    scalar_hc.eta = hc_state.eta;
    return ares_hc_scalar_cpp(scalar_thermo, scalar_hc, x, cppargs);
}

static double ares_disp_cpp(const MixtureState &thermo, const DispersionPolynomialState &dispersion) {
    MixtureStateScalar<double> scalar_thermo;
    scalar_thermo.d = thermo.d;
    scalar_thermo.e_ij = thermo.e_ij;
    scalar_thermo.s_ij = thermo.s_ij;
    scalar_thermo.den = thermo.den;
    scalar_thermo.m_avg = thermo.m_avg;
    scalar_thermo.m2es3 = thermo.m2es3;
    scalar_thermo.m2e2s3 = thermo.m2e2s3;

    DispersionPolynomialStateScalar<double> scalar_dispersion;
    for (size_t i = 0; i < scalar_dispersion.a.size(); ++i) {
        scalar_dispersion.a[i] = dispersion.a[i];
        scalar_dispersion.b[i] = dispersion.b[i];
    }
    scalar_dispersion.I1 = dispersion.I1;
    scalar_dispersion.I2 = dispersion.I2;
    scalar_dispersion.dEtaI1_deta = dispersion.dEtaI1_deta;
    scalar_dispersion.dEtaI2_deta = dispersion.dEtaI2_deta;
    scalar_dispersion.C1 = dispersion.C1;
    scalar_dispersion.C2 = dispersion.C2;
    return ares_disp_scalar_cpp(scalar_thermo, scalar_dispersion);
}

// EqID: ares_assoc
template <typename Scalar>
static Scalar ares_assoc_scalar_cpp(const vector<Scalar> &x, const add_args &cppargs) {
    (void)x;
    if (!cppargs.assoc_num.empty()) {
        throw ValueError("backend_unavailable: CppAD association recording is unavailable because site fractions require an iterative solve.");
    }
    return scalar_constant<Scalar>(0.0);
}

static double ares_assoc_cpp(const AssociationIntermediateState &assoc_state, const vector<double> &x) {
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
template <typename Scalar>
static Scalar ares_ion_scalar_cpp(double t, const MixtureStateScalar<Scalar> &thermo, const vector<Scalar> &x, const add_args &cppargs) {
    if (cppargs.z.empty()) {
        return scalar_constant<Scalar>(0.0);
    }
    Scalar q2_sum = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        q2_sum += x[i] * cppargs.z[i] * cppargs.z[i];
    }
    Scalar eps = scalar_constant<Scalar>(0.0);
    if (cppargs.dielc_rule == 0) {
        eps = scalar_constant<Scalar>(*std::max_element(cppargs.dielc.begin(), cppargs.dielc.end()));
    } else if (cppargs.dielc_rule == 1) {
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            eps += x[i] * cppargs.dielc[i];
        }
    } else {
        throw ValueError("backend_unavailable: CppAD ionic recording currently supports dielc_rule 0 or 1.");
    }
    Scalar kappa = scalar_sqrt(thermo.den * E_CHRG * E_CHRG / kb / t / (eps * perm_vac) * q2_sum);
    Scalar chi_sum = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        Scalar ka = kappa * thermo.d[i];
        Scalar chi = 3.0 / scalar_pow(ka, 3) * (1.5 + scalar_log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * scalar_pow(1.0 + ka, 2));
        chi_sum += x[i] * cppargs.z[i] * cppargs.z[i] * chi;
    }
    double K0 = E_CHRG * E_CHRG / (12.0 * PI * kb * t * perm_vac);
    return -K0 * kappa / eps * chi_sum;
}

static double ares_ion_cpp(double t, const IonIntermediateState &ion_state) {
    if (!ion_state.active) {
        return 0.0;
    }
    double K0 = E_CHRG * E_CHRG / (12.0 * PI * kb * t * perm_vac);
    return -K0 * ion_state.kappa / ion_state.dielectric.eps * ion_state.chi_sum;
}

// EqID: ares_born
template <typename Scalar>
static Scalar ares_born_scalar_cpp(double t, const vector<Scalar> &x, const add_args &cppargs) {
    if (cppargs.z.empty() || cppargs.born_model == 0) {
        return scalar_constant<Scalar>(0.0);
    }
    if (cppargs.born_model != 1) {
        throw ValueError("backend_unavailable: CppAD Born recording supports direct Born model=1 formulas only.");
    }
    Scalar eps = scalar_constant<Scalar>(0.0);
    if (cppargs.dielc_rule == 0) {
        eps = scalar_constant<Scalar>(*std::max_element(cppargs.dielc.begin(), cppargs.dielc.end()));
    } else if (cppargs.dielc_rule == 1) {
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            eps += x[i] * cppargs.dielc[i];
        }
    } else {
        throw ValueError("backend_unavailable: CppAD Born recording currently supports dielc_rule 0 or 1.");
    }
    if (cppargs.born_eps_mode == 1) {
        throw ValueError("backend_unavailable: CppAD Born reference-solvent dielectric routing is not implemented.");
    }
    Scalar charge_radius_sum = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        if (is_ion_species(cppargs, i)) {
            double d_born_i = thermo_detail::parameter_setup_detail::ion_born_radius_cpp(i, t, cppargs);
            charge_radius_sum += x[i] * cppargs.z[i] * cppargs.z[i] / d_born_i;
        }
    }
    return -E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac) * (1.0 - 1.0 / eps) * charge_radius_sum;
}

static double ares_born_cpp(double t, const BornIntermediateState &born_state) {
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

template <typename Scalar>
static AresContributionsScalar<Scalar> ares_contributions_scalar_cpp(double t, const Scalar &rho, const vector<Scalar> &x, const add_args &cppargs) {
    AresContributionsScalar<Scalar> out;
    MixtureStateScalar<Scalar> thermo = mixture_state_scalar_cpp(t, rho, x, cppargs);
    HardChainStateScalar<Scalar> hc_state = hard_chain_state_scalar_cpp(thermo, x, cppargs);
    DispersionPolynomialStateScalar<Scalar> dispersion = dispersion_polynomials_scalar_cpp(thermo.m_avg, hc_state.eta);
    out.hc = ares_hc_scalar_cpp(thermo, hc_state, x, cppargs);
    out.disp = ares_disp_scalar_cpp(thermo, dispersion);
    out.assoc = ares_assoc_scalar_cpp(x, cppargs);
    out.ion = ares_ion_scalar_cpp(t, thermo, x, cppargs);
    out.born = ares_born_scalar_cpp(t, x, cppargs);
    return out;
}

}  // namespace ares_detail

namespace {

std::string contribution_backend_name(int mode) {
    if (mode == 0) return "analytic";
    if (mode == 1) return "backend_unavailable";
    if (mode == 2) return "autodiff";
    if (mode == 3 || mode == 5) {
        return "analytic";
    }
    if (mode == 4) return "autodiff";
    return "unknown";
}

std::map<std::string, std::string> composition_derivative_backend_map(const add_args &cppargs) {
    std::map<std::string, std::string> backends;
    backends["hc"] = contribution_backend_name(cppargs.hc_dadx_diff_mode);
    backends["disp"] = contribution_backend_name(cppargs.disp_dadx_diff_mode);
    backends["assoc"] = contribution_backend_name(cppargs.assoc_dadx_diff_mode);
    backends["ion"] = contribution_backend_name(cppargs.mu_DH_diff_mode);
    backends["born"] = contribution_backend_name(cppargs.born_diff_mode);
    return backends;
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

    out.hc = ares_detail::ares_hc_cpp(thermo, hc_state, x, cppargs);
    out.disp = ares_detail::ares_disp_cpp(thermo, dispersion);
    out.assoc = ares_detail::ares_assoc_cpp(assoc_state, x);
    out.ion = ares_detail::ares_ion_cpp(t, ion_state);
    out.born = ares_detail::ares_born_cpp(t, born_state);
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

epcsaft::native::autodiff::ADDerivativeResult cppad_eos_contribution_derivatives_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs
) {
#ifdef EPCSAFT_HAS_CPPAD
    using CppADScalar = CppAD::AD<double>;
    if (!cppargs.assoc_num.empty()) {
        throw ValueError("backend_unavailable: CppAD association recording is unavailable because site fractions require an iterative solve.");
    }
    if (!cppargs.z.empty() && cppargs.born_model > 1) {
        throw ValueError("backend_unavailable: CppAD Born recording supports direct Born model=1 formulas only.");
    }
    int ncomp = static_cast<int>(x.size());
    std::vector<CppADScalar> ax(ncomp);
    for (int i = 0; i < ncomp; ++i) {
        ax[i] = x[i];
    }
    CppAD::Independent(ax);

    CppADScalar arho = rho;
    auto contributions = ares_detail::ares_contributions_scalar_cpp(t, arho, ax, cppargs);
    std::vector<CppADScalar> ay(6);
    ay[0] = contributions.hc;
    ay[1] = contributions.disp;
    ay[2] = contributions.assoc;
    ay[3] = contributions.ion;
    ay[4] = contributions.born;
    ay[5] = contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born;

    CppAD::ADFun<double> function(ax, ay);
    std::vector<double> point(x.begin(), x.end());
    auto value = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);

    epcsaft::native::autodiff::ADDerivativeResult result;
    result.supported = true;
    result.backend = "cppad";
    result.message = "CppAD EOS contribution composition derivatives available";
    result.value = std::move(value);
    result.jacobian_row_major = std::move(jacobian);
    result.rows = 6;
    result.cols = ncomp;
    return result;
#else
    (void)t;
    (void)rho;
    (void)x;
    (void)cppargs;
    epcsaft::native::autodiff::ADDerivativeResult result;
    result.supported = false;
    result.backend = "backend_unavailable";
    result.message = "CppAD support is disabled in this native build";
    return result;
#endif
}

DadrhoResult dadrho_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    MixtureState thermo = mixture_state_cpp(t, rho, x, cppargs, false);
    HardChainState hc_state = hard_chain_state_cpp(thermo, x, cppargs);
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(thermo.m_avg, hc_state.eta);
    AssociationIntermediateState assoc_state = association_intermediate_state_cpp(thermo, hc_state, t, x, cppargs, false, false);
    IonIntermediateState ion_state = ion_intermediate_state_cpp(thermo, t, x, cppargs, false);

    double hc = dadrho_hc_cpp(thermo, hc_state, x, cppargs);
    double disp = dadrho_disp_cpp(thermo, hc_state, dispersion);
    double assoc = dadrho_assoc_cpp(thermo, hc_state, assoc_state, x, cppargs, t);
    double ion = dadrho_ion_cpp(t, ion_state);
    double born = dadrho_born_cpp();
    double total = hc + disp + assoc + ion + born;

    ScalarContributionTerms raw_terms = make_scalar_terms(hc, disp, assoc, ion, born, total);
    DadrhoResult result;
    result.raw = raw_terms;
    result.terms = normalized_dadrho_terms_cpp(raw_terms);
    return result;
}

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

CompositionContributionResult composition_derivative_residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    MixtureState thermo = mixture_state_cpp(t, rho, x, cppargs, false);
    HardChainState hc_state = hard_chain_state_cpp(thermo, x, cppargs);
    DadrhoResult dadrho_result = dadrho_result_cpp(t, rho, x, cppargs);
    const ScalarContributionTerms &z_raw_terms = dadrho_result.raw;
    ScalarContributionTerms z_terms = compressibility_terms_from_dadrho_cpp(dadrho_result);

    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(thermo.m_avg, hc_state.eta);
    AssociationIntermediateState assoc_state = association_intermediate_state_cpp(thermo, hc_state, t, x, cppargs, false, false);
    IonIntermediateState ion_state = ion_intermediate_state_cpp(thermo, t, x, cppargs, true);
    BornIntermediateState born_state = born_intermediate_state_cpp(t, x, cppargs, false, true);

    ContributionDadxResult hc = dadx_hc_cpp(thermo, hc_state, t, rho, x, cppargs);
    ContributionDadxResult disp = dadx_disp_cpp(thermo, hc_state, dispersion, t, rho, x, cppargs);
    ContributionDadxResult assoc = dadx_assoc_cpp(thermo, hc_state, assoc_state, t, rho, x, cppargs);
    ContributionDadxResult ion = dadx_ion_cpp(thermo, ion_state, t, rho, x, cppargs);
    ContributionDadxResult born = dadx_born_cpp(born_state, t, rho, x, cppargs);

    CompositionContributionResult result;
    result.dadx = make_vector_terms(hc.dadx, disp.dadx, assoc.dadx, ion.dadx, born.dadx, vector<double>());
    result.ares = make_scalar_terms(hc.ares, disp.ares, assoc.ares, ion.ares, born.ares,
        hc.ares + disp.ares + assoc.ares + ion.ares + born.ares);
    result.sum_x_dadx = make_scalar_terms(hc.sum_x_dadx, disp.sum_x_dadx, assoc.sum_x_dadx,
        ion.sum_x_dadx, born.sum_x_dadx,
        hc.sum_x_dadx + disp.sum_x_dadx + assoc.sum_x_dadx + ion.sum_x_dadx + born.sum_x_dadx);
    result.z_raw = z_raw_terms;
    result.z = z_terms;

    vector<double> total(ncomp, 0.0);
    for (int i = 0; i < ncomp; ++i) {
        total[i] = hc.dadx[i] + disp.dadx[i] + assoc.dadx[i] + ion.dadx[i] + born.dadx[i];
    }
    result.dadx.total = total;
    result.derivative_backend = composition_derivative_backend_map(cppargs);
    result.derivative_available = true;
    for (const auto& item : result.derivative_backend) {
        if (item.second == "backend_unavailable") {
            result.derivative_available = false;
            result.backend_unavailable_reason = "backend_unavailable";
            break;
        }
    }
    return result;
}
