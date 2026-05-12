#include "epcsaft_core_internal.h"
#include "contributions/epcsaft_contrib_internal.h"
#include "autodiff/ad_scalar.h"

using namespace thermo_detail;

namespace {

template <typename Scalar>
struct PressureMixtureStateScalar {
    vector<double> d;
    Scalar den = scalar_constant<Scalar>(0.0);
    Scalar m_avg = scalar_constant<Scalar>(0.0);
    Scalar m2es3 = scalar_constant<Scalar>(0.0);
    Scalar m2e2s3 = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
struct PressureHardChainStateScalar {
    vector<Scalar> zeta;
    vector<Scalar> ghs;
    Scalar eta = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
PressureHardChainStateScalar<Scalar> pressure_hard_chain_state_scalar_cpp(
    const Scalar &den,
    const vector<double> &d,
    const vector<Scalar> &x,
    const add_args &cppargs
) {
    PressureHardChainStateScalar<Scalar> state;
    const int ncomp = static_cast<int>(x.size());
    state.zeta.assign(4, scalar_constant<Scalar>(0.0));
    state.ghs.assign(static_cast<size_t>(ncomp * ncomp), scalar_constant<Scalar>(0.0));
    for (int k = 0; k < 4; ++k) {
        Scalar summ = scalar_constant<Scalar>(0.0);
        for (int j = 0; j < ncomp; ++j) {
            summ += x[static_cast<size_t>(j)] * cppargs.m[static_cast<size_t>(j)] * scalar_pow(d[static_cast<size_t>(j)], k);
        }
        state.zeta[static_cast<size_t>(k)] = PI / 6.0 * den * summ;
    }
    state.eta = state.zeta[3];
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            const double pair_diameter = parameter_setup_detail::pair_diameter_cpp(d[static_cast<size_t>(i)], d[static_cast<size_t>(j)]);
            state.ghs[static_cast<size_t>(idx)] = 1.0 / (1.0 - state.zeta[3])
                + pair_diameter * 3.0 * state.zeta[2] / scalar_pow(1.0 - state.zeta[3], 2)
                + scalar_pow(pair_diameter, 2.0) * 2.0 * state.zeta[2] * state.zeta[2] / scalar_pow(1.0 - state.zeta[3], 3);
        }
    }
    return state;
}

template <typename Scalar>
PressureMixtureStateScalar<Scalar> pressure_mixture_state_scalar_cpp(
    double t,
    const Scalar &rho,
    const vector<Scalar> &x,
    const add_args &cppargs
) {
    const int ncomp = static_cast<int>(x.size());
    PressureMixtureStateScalar<Scalar> state;
    state.d.assign(static_cast<size_t>(ncomp), 0.0);
    state.den = rho * (N_AV / 1.0e30);

    for (int i = 0; i < ncomp; ++i) {
        state.d[static_cast<size_t>(i)] = cppargs.s[static_cast<size_t>(i)]
            * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[static_cast<size_t>(i)] / t));
        if (!cppargs.z.empty() && std::abs(cppargs.z[static_cast<size_t>(i)]) > 1.0e-12) {
            state.d[static_cast<size_t>(i)] = parameter_setup_detail::ion_diameter_cpp(i, t, cppargs);
        }
        state.m_avg += x[static_cast<size_t>(i)] * cppargs.m[static_cast<size_t>(i)];
    }

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            const double sigma_ij = parameter_setup_detail::pair_sigma_cpp(static_cast<size_t>(idx), i, j, cppargs);
            const double epsilon_ij = parameter_setup_detail::pair_epsilon_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.m2es3 += x[static_cast<size_t>(i)] * x[static_cast<size_t>(j)]
                * cppargs.m[static_cast<size_t>(i)] * cppargs.m[static_cast<size_t>(j)] * epsilon_ij / t
                * scalar_pow(sigma_ij, 3.0);
            state.m2e2s3 += x[static_cast<size_t>(i)] * x[static_cast<size_t>(j)]
                * cppargs.m[static_cast<size_t>(i)] * cppargs.m[static_cast<size_t>(j)] * scalar_pow(epsilon_ij / t, 2.0)
                * scalar_pow(sigma_ij, 3.0);
        }
    }
    return state;
}

template <typename Scalar>
Scalar dadrho_hs_scalar_cpp(const vector<Scalar> &zeta) {
    const Scalar one = scalar_constant<Scalar>(1.0);
    return zeta[3] / (one - zeta[3])
        + 3.0 * zeta[1] * zeta[2] / zeta[0] / scalar_pow(one - zeta[3], 2)
        + (3.0 * scalar_pow(zeta[2], 3) - zeta[3] * scalar_pow(zeta[2], 3)) / zeta[0] / scalar_pow(one - zeta[3], 3);
}

template <typename Scalar>
Scalar hs_contact_density_derivative_scalar_cpp(double pair_diameter, const Scalar &zeta2, const Scalar &zeta3) {
    const Scalar one = scalar_constant<Scalar>(1.0);
    return zeta3 / scalar_pow(one - zeta3, 2)
        + pair_diameter * (3.0 * zeta2 / scalar_pow(one - zeta3, 2) + 6.0 * zeta2 * zeta3 / scalar_pow(one - zeta3, 3))
        + scalar_pow(pair_diameter, 2.0)
            * (4.0 * zeta2 * zeta2 / scalar_pow(one - zeta3, 3) + 6.0 * zeta2 * zeta2 * zeta3 / scalar_pow(one - zeta3, 4));
}

template <typename Scalar>
struct DispersionAutodiffState {
    std::array<Scalar, 7> a{};
    std::array<Scalar, 7> b{};
    Scalar dEtaI1_deta = scalar_constant<Scalar>(0.0);
    Scalar dEtaI2_deta = scalar_constant<Scalar>(0.0);
    Scalar I2 = scalar_constant<Scalar>(0.0);
    Scalar C1 = scalar_constant<Scalar>(0.0);
    Scalar C2 = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
DispersionAutodiffState<Scalar> dispersion_autodiff_state_cpp(const Scalar &m_avg, const Scalar &eta) {
    DispersionAutodiffState<Scalar> state;
    const Scalar one = scalar_constant<Scalar>(1.0);
    const Scalar c1 = (m_avg - one) / m_avg;
    const Scalar c2 = (m_avg - scalar_constant<Scalar>(2.0)) / m_avg;
    for (size_t i = 0; i < state.a.size(); ++i) {
        state.a[i] = kDispersionA0[i] + c1 * kDispersionA1[i] + c1 * c2 * kDispersionA2[i];
        state.b[i] = kDispersionB0[i] + c1 * kDispersionB1[i] + c1 * c2 * kDispersionB2[i];
        state.I2 += state.b[i] * scalar_pow(eta, static_cast<int>(i));
        state.dEtaI1_deta += state.a[i] * static_cast<double>(i + 1) * scalar_pow(eta, static_cast<int>(i));
        state.dEtaI2_deta += state.b[i] * static_cast<double>(i + 1) * scalar_pow(eta, static_cast<int>(i));
    }
    state.C1 = one / (
        one
        + m_avg * (8.0 * eta - 2.0 * eta * eta) / scalar_pow(one - eta, 4)
        + (one - m_avg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * scalar_pow(eta, 3) - 2.0 * scalar_pow(eta, 4))
            / scalar_pow((one - eta) * (scalar_constant<Scalar>(2.0) - eta), 2)
    );
    state.C2 = -state.C1 * state.C1 * (
        m_avg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / scalar_pow(one - eta, 5)
        + (one - m_avg) * (2.0 * scalar_pow(eta, 3) + 12.0 * eta * eta - 48.0 * eta + 40.0)
            / scalar_pow((one - eta) * (scalar_constant<Scalar>(2.0) - eta), 3)
    );
    return state;
}

template <typename Scalar>
Scalar pressure_dielectric_constant_supported_scalar_cpp(const vector<Scalar>& x, const add_args& cppargs) {
    if (cppargs.dielc_rule == 0) {
        return scalar_constant<Scalar>(*std::max_element(cppargs.dielc.begin(), cppargs.dielc.end()));
    }
    if (cppargs.dielc_rule == 1) {
        Scalar eps = scalar_constant<Scalar>(0.0);
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            eps += x[static_cast<size_t>(i)] * cppargs.dielc[static_cast<size_t>(i)];
        }
        return eps;
    }
    throw ValueError("native pressure derivatives currently support dielc_rule 0 or 1 only.");
}

bool pressure_composition_derivative_supported_cpp(const add_args &cppargs, std::string *reason) {
    if (!cppargs.assoc_num.empty() || !cppargs.assoc_matrix.empty() || !cppargs.k_hb.empty()
        || !cppargs.e_assoc.empty() || !cppargs.vol_a.empty()) {
        if (reason != nullptr) {
            *reason = "native pressure-composition derivatives currently support nonassociating states only.";
        }
        return false;
    }
    if (cppargs.DH_model == 2) {
        if (reason != nullptr) {
            *reason = "native pressure-composition derivatives do not support DH_model=2.";
        }
        return false;
    }
    if (cppargs.dielc_rule != 0 && cppargs.dielc_rule != 1) {
        if (reason != nullptr) {
            *reason = "native pressure-composition derivatives currently support dielc_rule 0 or 1 only.";
        }
        return false;
    }
    return true;
}

template <typename Scalar>
Scalar pressure_scalar_supported_cpp(double t, const Scalar &rho, const vector<Scalar> &x, const add_args &cppargs) {
    const int ncomp = static_cast<int>(x.size());
    PressureMixtureStateScalar<Scalar> thermo = pressure_mixture_state_scalar_cpp(t, rho, x, cppargs);
    PressureHardChainStateScalar<Scalar> hc_state = pressure_hard_chain_state_scalar_cpp(thermo.den, thermo.d, x, cppargs);

    Scalar z_hc = thermo.m_avg * dadrho_hs_scalar_cpp(hc_state.zeta);
    for (int i = 0; i < ncomp; ++i) {
        const double pair_diameter = parameter_setup_detail::pair_diameter_cpp(thermo.d[static_cast<size_t>(i)], thermo.d[static_cast<size_t>(i)]);
        z_hc -= x[static_cast<size_t>(i)] * (cppargs.m[static_cast<size_t>(i)] - 1.0) / hc_state.ghs[static_cast<size_t>(i * ncomp + i)]
            * hs_contact_density_derivative_scalar_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3]);
    }

    DispersionAutodiffState<Scalar> dispersion = dispersion_autodiff_state_cpp(thermo.m_avg, hc_state.eta);
    Scalar z_disp = -2.0 * PI * thermo.den * dispersion.dEtaI1_deta * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * (dispersion.C1 * dispersion.dEtaI2_deta + dispersion.C2 * hc_state.eta * dispersion.I2) * thermo.m2e2s3;

    Scalar z_ion = scalar_constant<Scalar>(0.0);
    if (!cppargs.z.empty()) {
        Scalar q2_sum = scalar_constant<Scalar>(0.0);
        for (int i = 0; i < ncomp; ++i) {
            q2_sum += x[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)];
        }
        if (scalar_value(q2_sum) > 0.0) {
            Scalar eps = pressure_dielectric_constant_supported_scalar_cpp(x, cppargs);
            Scalar kappa = scalar_sqrt(thermo.den * E_CHRG * E_CHRG / kb / t / perm_vac * q2_sum / eps);
            Scalar sigma_sum = scalar_constant<Scalar>(0.0);
            for (int i = 0; i < ncomp; ++i) {
                const double d_i = parameter_setup_detail::ion_diameter_cpp(i, t, cppargs);
                Scalar ka = kappa * d_i;
                Scalar chi = 3.0 / scalar_pow(ka, 3)
                    * (1.5 + scalar_log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * scalar_pow(1.0 + ka, 2));
                Scalar sigma_k = -2.0 * chi + 3.0 / (1.0 + ka);
                sigma_sum += x[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * cppargs.z[static_cast<size_t>(i)] * sigma_k;
            }
            z_ion = -kappa / 24.0 / PI / kb / t / (eps * perm_vac) * sigma_sum * E_CHRG * E_CHRG;
        }
    }

    const Scalar z_total = scalar_constant<Scalar>(1.0) + z_hc + z_disp + z_ion;
    return z_total * kb * t * thermo.den * 1.0e30;
}

#ifdef EPCSAFT_HAS_CPPAD
std::vector<double> pressure_composition_jacobian_cppad_cpp(
    double t,
    double rho,
    const vector<double>& x,
    const add_args& cppargs
) {
    using epcsaft::autodiff::CppADScalar;
    const std::size_t ncomp = x.size();
    std::vector<CppADScalar> independent(ncomp);
    for (std::size_t i = 0; i < ncomp; ++i) {
        independent[i] = x[i];
    }
    CppAD::Independent(independent);
    std::vector<CppADScalar> dependent(1);
    dependent[0] = pressure_scalar_supported_cpp(t, CppADScalar(rho), independent, cppargs);
    CppAD::ADFun<double> tape(independent, dependent);
    return tape.Jacobian(x);
}

double pressure_density_jacobian_cppad_cpp(
    double t,
    double rho,
    const vector<double>& x,
    const add_args& cppargs
) {
    using epcsaft::autodiff::CppADScalar;
    std::vector<CppADScalar> independent(1);
    independent[0] = rho;
    CppAD::Independent(independent);
    std::vector<CppADScalar> x_const(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        x_const[i] = x[i];
    }
    std::vector<CppADScalar> dependent(1);
    dependent[0] = pressure_scalar_supported_cpp(t, independent[0], x_const, cppargs);
    CppAD::ADFun<double> tape(independent, dependent);
    const std::vector<double> jacobian = tape.Jacobian(std::vector<double>{rho});
    return jacobian[0];
}
#endif

}  // namespace

double z_term_scale_cpp(const vector<double> &z_term, double increment_total) {
    double raw_sum = 0.0;
    for (double value : z_term) {
        raw_sum += value;
    }
    if (std::abs(raw_sum) <= 1e-14) {
        if (std::abs(increment_total) <= 1e-12) {
            return 0.0;
        }
        throw ValueError("Could not normalize density-derivative terms because their sum is ~0 while the compressibility increment is non-zero.");
    }
    return increment_total / raw_sum;
}

double normalized_dadrho_scale_cpp(const ScalarContributionTerms &raw_terms) {
    vector<double> raw = {
        raw_terms.hc,
        raw_terms.disp,
        raw_terms.assoc,
        raw_terms.ion,
        raw_terms.born
    };
    return z_term_scale_cpp(raw, raw_terms.total);
}

double normalized_dadrho_term_cpp(double raw_term, double scale) {
    return raw_term * scale;
}

double z_total_cpp(double dadrho_total) {
    return 1.0 + dadrho_total;
}

ScalarContributionTerms normalized_dadrho_terms_cpp(const ScalarContributionTerms &raw_terms) {
    double scale = normalized_dadrho_scale_cpp(raw_terms);
    return make_scalar_terms(
        normalized_dadrho_term_cpp(raw_terms.hc, scale),
        normalized_dadrho_term_cpp(raw_terms.disp, scale),
        normalized_dadrho_term_cpp(raw_terms.assoc, scale),
        normalized_dadrho_term_cpp(raw_terms.ion, scale),
        normalized_dadrho_term_cpp(raw_terms.born, scale),
        raw_terms.total
    );
}

// EqID: z_alpha
ScalarContributionTerms compressibility_terms_from_dadrho_cpp(const DadrhoResult &result) {
    return make_scalar_terms(
        result.terms.hc,
        result.terms.disp,
        result.terms.assoc,
        result.terms.ion,
        result.terms.born,
        z_total_cpp(result.terms.total)
    );
}

// EqID: z_from_rho
// EqID: z_total
CompressibilityFactorResult compressibility_factor_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    DadrhoResult dadrho_result = dadrho_result_cpp(t, rho, std::move(x), cppargs);
    CompressibilityFactorResult result;
    result.raw = dadrho_result.raw;
    result.terms = compressibility_terms_from_dadrho_cpp(dadrho_result);
    return result;
}

double Z_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return compressibility_factor_result_cpp(t, rho, std::move(x), cppargs).terms.total;
}

// EqID: pressure_from_z
double p_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    double den = rho * N_AV / 1.0e30;
    double Z = Z_cpp(t, rho, std::move(x), cppargs);
    return Z * kb * t * den * 1.0e30;
}

PressureCompositionDerivativeResult pressure_composition_derivative_result_cpp(
    double t,
    double rho,
    vector<double> x,
    const add_args &cppargs
) {
    PressureCompositionDerivativeResult result;
    result.pressure = p_cpp(t, rho, x, cppargs);
    std::string unsupported_reason;
    if (!pressure_composition_derivative_supported_cpp(cppargs, &unsupported_reason)) {
        result.supported = false;
        result.derivative_backend = "unsupported";
        result.unsupported_derivative_fallback_reason = unsupported_reason;
        result.dpdx.assign(x.size(), std::numeric_limits<double>::quiet_NaN());
        return result;
    }

    const int ncomp = static_cast<int>(x.size());
    result.dpdx.assign(ncomp, 0.0);
#ifdef EPCSAFT_HAS_CPPAD
    const std::vector<double> jacobian = pressure_composition_jacobian_cppad_cpp(t, rho, x, cppargs);
    for (int i = 0; i < ncomp; ++i) {
        result.dpdx[static_cast<size_t>(i)] = jacobian[static_cast<size_t>(i)];
        if (!std::isfinite(result.dpdx[static_cast<size_t>(i)])) {
            throw ValueError("Non-finite native pressure-composition derivative.");
        }
    }
    result.derivative_backend = "cppad_composition";
#else
    for (int i = 0; i < ncomp; ++i) {
        vector<AutoDual> x_dual(static_cast<size_t>(ncomp), make_autodiff_scalar(0.0, 0.0));
        for (int j = 0; j < ncomp; ++j) {
            x_dual[static_cast<size_t>(j)] = make_autodiff_scalar(x[static_cast<size_t>(j)], (i == j) ? 1.0 : 0.0);
        }
        AutoDual pressure = pressure_autodiff_supported_cpp(t, rho, x_dual, cppargs);
        result.dpdx[static_cast<size_t>(i)] = scalar_derivative(pressure);
        if (!std::isfinite(result.dpdx[static_cast<size_t>(i)])) {
            throw ValueError("Non-finite native pressure-composition derivative.");
        }
    }
    result.derivative_backend = "autodiff_composition";
#endif
    result.supported = true;
    return result;
}

PressureDensityDerivativeResult pressure_density_derivative_result_cpp(
    double t,
    double rho,
    vector<double> x,
    const add_args &cppargs
) {
    PressureDensityDerivativeResult result;
    result.pressure = p_cpp(t, rho, x, cppargs);
    std::string unsupported_reason;
    if (!pressure_composition_derivative_supported_cpp(cppargs, &unsupported_reason)) {
        result.supported = false;
        result.derivative_backend = "unsupported";
        result.unsupported_derivative_fallback_reason = unsupported_reason;
        result.dpdrho = std::numeric_limits<double>::quiet_NaN();
        return result;
    }

#ifdef EPCSAFT_HAS_CPPAD
    result.dpdrho = pressure_density_jacobian_cppad_cpp(t, rho, x, cppargs);
    result.derivative_backend = "cppad_density";
#else
    const int ncomp = static_cast<int>(x.size());
    vector<AutoDual> x_const(static_cast<size_t>(ncomp), make_autodiff_scalar(0.0, 0.0));
    for (int i = 0; i < ncomp; ++i) {
        x_const[static_cast<size_t>(i)] = make_autodiff_scalar(x[static_cast<size_t>(i)], 0.0);
    }
    AutoDual pressure = pressure_scalar_supported_cpp(t, make_autodiff_scalar(rho, 1.0), x_const, cppargs);
    result.dpdrho = scalar_derivative(pressure);
    result.derivative_backend = "autodiff_density";
#endif
    if (!std::isfinite(result.dpdrho)) {
        throw ValueError("Non-finite native pressure-density derivative.");
    }
    result.supported = true;
    return result;
}



