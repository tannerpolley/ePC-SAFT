#include "epcsaft_core_internal.h"
#include "contributions/epcsaft_contrib_internal.h"

#include <cppad/cppad.hpp>

#include <map>
#include <numeric>
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

constexpr int kGenericTargetSLocal = 1;
constexpr int kGenericTargetELocal = 2;
constexpr int kGenericTargetDBornLocal = 5;
constexpr int kGenericTargetFSolvLocal = 9;
constexpr int kGenericTargetDielcLocal = 10;

template <typename Scalar>
struct MixtureStateScalar {
    vector<Scalar> d;
    vector<Scalar> e_ij;
    vector<Scalar> s_ij;
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
struct AssociationImplicitTermsScalar {
    vector<int> site_component_index;
    vector<Scalar> x_assoc;
    vector<Scalar> delta_ij;
    vector<Scalar> residuals;
    Scalar ares = scalar_constant<Scalar>(0.0);
    Scalar zraw = scalar_constant<Scalar>(0.0);
    vector<Scalar> dadx;
    Scalar sum_x_dadx = scalar_constant<Scalar>(0.0);
    vector<Scalar> mu;
};

struct AssociationDensityResponse {
    double zraw = 0.0;
    vector<double> mu;
    double dzraw_drho = 0.0;
    vector<double> dmu_drho;
};

struct AssociationPhaseStateResponse {
    bool active = false;
    int ncomp = 0;
    int base_var_count = 0;
    double zraw = 0.0;
    vector<double> mu;
    vector<double> dzraw_dvar;
    vector<double> dmu_dvar_row_major;
    vector<double> d2zraw_dvar2_row_major;
    vector<double> d2mu_dvar2_tensor_row_major;
};

std::vector<double> vector_output_component_hessian_cppad(
    CppAD::ADFun<double>& function,
    const std::vector<double>& point,
    std::size_t output_index
) {
    return function.Hessian(point, output_index);
}

std::vector<double> scalar_function_third_derivative_tensor_cppad(
    const CppAD::ADFun<double>& function,
    const std::vector<double>& point
) {
    using ADScalar = CppAD::AD<double>;
    const std::size_t nvar = point.size();
    std::vector<ADScalar> ax(nvar);
    for (std::size_t i = 0; i < nvar; ++i) {
        ax[i] = point[i];
    }
    CppAD::Independent(ax);
    CppAD::ADFun<ADScalar, double> af(function.base2ad());
    std::vector<ADScalar> ay = af.Forward(0, ax);
    if (ay.size() != 1) {
        throw ValueError("Third-derivative tensor helper requires a scalar recorded function.");
    }
    std::vector<ADScalar> aw(1);
    aw[0] = ADScalar(1.0);
    std::vector<ADScalar> agrad = af.Reverse(1, aw);
    CppAD::ADFun<double> grad_fun(ax, agrad);
    std::vector<double> tensor(nvar * nvar * nvar, 0.0);
    for (std::size_t component = 0; component < nvar; ++component) {
        const std::vector<double> hessian = grad_fun.Hessian(point, component);
        for (std::size_t row = 0; row < nvar; ++row) {
            for (std::size_t col = 0; col < nvar; ++col) {
                tensor[component * nvar * nvar + row * nvar + col] =
                    hessian[row * nvar + col];
            }
        }
    }
    return tensor;
}

template <typename Scalar>
static Scalar scalar_component_parameter_cpp(
    int target_kind,
    int target_index,
    int component_index,
    double base_value,
    const Scalar *target_value
) {
    if (target_value != nullptr && target_index == component_index && target_kind >= 0) {
        return *target_value;
    }
    return scalar_constant<Scalar>(base_value);
}

template <typename Scalar>
static Scalar scalar_dielc_parameter_cpp(
    int target_kind,
    int target_index,
    int component_index,
    double base_value,
    const Scalar *target_value
) {
    return scalar_component_parameter_cpp(
        target_kind,
        target_index,
        component_index,
        base_value,
        target_kind == kGenericTargetDielcLocal ? target_value : nullptr
    );
}

template <typename Scalar, typename TemperatureScalar>
static Scalar ion_diameter_scalar_cpp(
    int i,
    const TemperatureScalar &t,
    const add_args &cppargs,
    const Scalar &sigma_i,
    const Scalar &epsilon_i
) {
    if (!is_ion_species(cppargs, i)) {
        return sigma_i;
    }
    int mode = cppargs.d_ion_mode;
    if (mode == 0) {
        return sigma_i;
    }
    if (mode == 1) {
        return sigma_i * 0.88;
    }
    if (mode == 2) {
        return sigma_i * (1.0 - 0.12 * scalar_exp(-3.0 * epsilon_i / t));
    }
    throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
}

template <typename Scalar, typename TemperatureScalar>
static Scalar ion_born_radius_scalar_cpp(
    int i,
    const TemperatureScalar &t,
    const add_args &cppargs,
    const Scalar &sigma_i,
    const Scalar &epsilon_i,
    int target_kind,
    int target_index,
    const Scalar *target_value
) {
    if (!is_ion_species(cppargs, i)) {
        return sigma_i;
    }
    int mode = cppargs.d_born_mode;
    if (mode == 0) {
        return sigma_i;
    }
    if (mode == 1) {
        return sigma_i * 0.88;
    }
    if (mode == 2) {
        return sigma_i * (1.0 - 0.12 * scalar_exp(-3.0 * epsilon_i / t));
    }
    if (mode == 3) {
        if (cppargs.d_born.size() <= static_cast<size_t>(i) || cppargs.d_born[i] <= 0.0) {
            throw ValueError("d_Born_mode=fitted_param requires positive ionic params['d_born'] values.");
        }
        return scalar_component_parameter_cpp(
            target_kind,
            target_index,
            i,
            cppargs.d_born[i],
            target_kind == kGenericTargetDBornLocal ? target_value : nullptr
        );
    }
    throw ValueError("Unknown d_Born_mode. Supported values are 0, 1, 2, 3.");
}

template <typename Scalar, typename TemperatureScalar>
static MixtureStateScalar<Scalar> mixture_state_scalar_cpp(
    const TemperatureScalar &t,
    const Scalar &rho,
    const vector<Scalar> &x,
    const add_args &cppargs,
    int k_override_index = -1,
    const Scalar *k_override_value = nullptr,
    int l_override_index = -1,
    const Scalar *l_override_value = nullptr,
    int component_target_kind = -1,
    int component_target_index = -1,
    const Scalar *component_target_value = nullptr
) {
    MixtureStateScalar<Scalar> state;
    int ncomp = static_cast<int>(x.size());
    state.d.assign(ncomp, scalar_constant<Scalar>(0.0));
    state.e_ij.assign(ncomp * ncomp, 0.0);
    state.s_ij.assign(ncomp * ncomp, scalar_constant<Scalar>(0.0));
    state.den = rho * N_AV / 1.0e30;

    for (int i = 0; i < ncomp; ++i) {
        Scalar sigma_i = scalar_component_parameter_cpp(
            component_target_kind,
            component_target_index,
            i,
            cppargs.s[i],
            component_target_kind == kGenericTargetSLocal ? component_target_value : nullptr
        );
        Scalar epsilon_i = scalar_component_parameter_cpp(
            component_target_kind,
            component_target_index,
            i,
            cppargs.e[i],
            component_target_kind == kGenericTargetELocal ? component_target_value : nullptr
        );
        state.d[i] = sigma_i * (1.0 - 0.12 * scalar_exp(-3.0 * epsilon_i / t));
        if (!cppargs.z.empty() && std::abs(cppargs.z[i]) > 1e-12) {
            state.d[i] = ion_diameter_scalar_cpp(i, t, cppargs, sigma_i, epsilon_i);
        }
        state.m_avg += x[i] * scalar_constant<Scalar>(cppargs.m[i]);
    }

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            if (l_override_value != nullptr && l_override_index >= 0) {
                Scalar sigma_i = scalar_component_parameter_cpp(
                    component_target_kind,
                    component_target_index,
                    i,
                    cppargs.s[i],
                    component_target_kind == kGenericTargetSLocal ? component_target_value : nullptr
                );
                Scalar sigma_j = scalar_component_parameter_cpp(
                    component_target_kind,
                    component_target_index,
                    j,
                    cppargs.s[j],
                    component_target_kind == kGenericTargetSLocal ? component_target_value : nullptr
                );
                state.s_ij[idx] = 0.5 * (sigma_i + sigma_j);
                if (idx == l_override_index) {
                    state.s_ij[idx] *= (1.0 - *l_override_value);
                } else if (!cppargs.l_ij.empty()) {
                    state.s_ij[idx] *= (1.0 - cppargs.l_ij[static_cast<size_t>(idx)]);
                }
            } else {
                Scalar sigma_i = scalar_component_parameter_cpp(
                    component_target_kind,
                    component_target_index,
                    i,
                    cppargs.s[i],
                    component_target_kind == kGenericTargetSLocal ? component_target_value : nullptr
                );
                Scalar sigma_j = scalar_component_parameter_cpp(
                    component_target_kind,
                    component_target_index,
                    j,
                    cppargs.s[j],
                    component_target_kind == kGenericTargetSLocal ? component_target_value : nullptr
                );
                state.s_ij[idx] = 0.5 * (sigma_i + sigma_j);
                if (!cppargs.l_ij.empty()) {
                    state.s_ij[idx] *= (1.0 - cppargs.l_ij[static_cast<size_t>(idx)]);
                }
            }
            Scalar epsilon_i = scalar_component_parameter_cpp(
                component_target_kind,
                component_target_index,
                i,
                cppargs.e[i],
                component_target_kind == kGenericTargetELocal ? component_target_value : nullptr
            );
            Scalar epsilon_j = scalar_component_parameter_cpp(
                component_target_kind,
                component_target_index,
                j,
                cppargs.e[j],
                component_target_kind == kGenericTargetELocal ? component_target_value : nullptr
            );
            Scalar epsilon = scalar_sqrt(epsilon_i * epsilon_j);
            if (!cppargs.z.empty() && cppargs.z[i] * cppargs.z[j] > 0.0) {
                epsilon = scalar_constant<Scalar>(0.0);
            } else if (k_override_value != nullptr && k_override_index >= 0) {
                if (idx == k_override_index) {
                    epsilon *= (1.0 - *k_override_value);
                } else if (!cppargs.k_ij.empty()) {
                    epsilon *= (1.0 - cppargs.k_ij[static_cast<size_t>(idx)]);
                }
            } else if (!cppargs.k_ij.empty()) {
                epsilon *= (1.0 - cppargs.k_ij[static_cast<size_t>(idx)]);
            }
            state.e_ij[idx] = epsilon;
            state.m2es3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * state.e_ij[idx] / t * scalar_pow(state.s_ij[idx], 3);
            state.m2e2s3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * scalar_pow(state.e_ij[idx] / t, 2) * scalar_pow(state.s_ij[idx], 3);
        }
    }
    return state;
}

template <typename Scalar>
static Scalar hs_contact_value_scalar_cpp(const Scalar &pair_diameter, const Scalar &zeta2, const Scalar &zeta3) {
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
            Scalar pair_diameter = thermo.d[i] * thermo.d[j] / (thermo.d[i] + thermo.d[j]);
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
        for (int sites : cppargs.assoc_num) {
            if (sites > 0) {
                throw ValueError("unsupported: CppAD association recording requires implicit site-fraction sensitivities.");
            }
        }
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

template <typename Scalar>
static Scalar hs_contact_density_derivative_scalar_cpp(const Scalar &pair_diameter, const Scalar &zeta2, const Scalar &zeta3) {
    return zeta3 / scalar_pow(1.0 - zeta3, 2.0)
        + pair_diameter * (
            3.0 * zeta2 / scalar_pow(1.0 - zeta3, 2.0)
            + 6.0 * zeta2 * zeta3 / scalar_pow(1.0 - zeta3, 3.0)
        )
        + scalar_pow(pair_diameter, 2.0) * (
            4.0 * zeta2 * zeta2 / scalar_pow(1.0 - zeta3, 3.0)
            + 6.0 * zeta2 * zeta2 * zeta3 / scalar_pow(1.0 - zeta3, 4.0)
        );
}

template <typename Scalar>
static Scalar hs_contact_composition_derivative_scalar_cpp(
    const Scalar &pair_diameter,
    const Scalar &zeta2,
    const Scalar &zeta3,
    const Scalar &dzeta2_dx,
    const Scalar &dzeta3_dx
) {
    return dzeta3_dx / scalar_pow(1.0 - zeta3, 2.0)
        + pair_diameter * (
            3.0 * dzeta2_dx / scalar_pow(1.0 - zeta3, 2.0)
            + 6.0 * zeta2 * dzeta3_dx / scalar_pow(1.0 - zeta3, 3.0)
        )
        + scalar_pow(pair_diameter, 2.0) * (
            4.0 * zeta2 * dzeta2_dx / scalar_pow(1.0 - zeta3, 3.0)
            + 6.0 * zeta2 * zeta2 * dzeta3_dx / scalar_pow(1.0 - zeta3, 4.0)
        );
}

template <typename Scalar>
static Scalar association_volume_scalar_cpp(int comp_i, int comp_j, int ncomp, const vector<Scalar> &s_ij, const add_args &cppargs) {
    const int idxi = comp_i * ncomp + comp_i;
    const int idxj = comp_j * ncomp + comp_j;
    Scalar volume = scalar_constant<Scalar>(std::sqrt(cppargs.vol_a[static_cast<size_t>(comp_i)] * cppargs.vol_a[static_cast<size_t>(comp_j)]))
        * scalar_pow(
            scalar_sqrt(s_ij[static_cast<size_t>(idxi)] * s_ij[static_cast<size_t>(idxj)])
                / (0.5 * (s_ij[static_cast<size_t>(idxi)] + s_ij[static_cast<size_t>(idxj)])),
            3.0
        );
    if (!cppargs.k_hb.empty()) {
        volume *= scalar_constant<Scalar>(1.0 - cppargs.k_hb[static_cast<size_t>(comp_i * ncomp + comp_j)]);
    }
    return volume;
}

template <typename Scalar>
static vector<Scalar> solve_linear_system_scalar_cpp(vector<Scalar> matrix, vector<Scalar> rhs, int n) {
    for (int col = 0; col < n; ++col) {
        int pivot = col;
        double pivot_abs = std::abs(scalar_value(matrix[static_cast<size_t>(col * n + col)]));
        for (int row = col + 1; row < n; ++row) {
            const double candidate_abs = std::abs(scalar_value(matrix[static_cast<size_t>(row * n + col)]));
            if (candidate_abs > pivot_abs) {
                pivot = row;
                pivot_abs = candidate_abs;
            }
        }
        if (pivot_abs <= 1.0e-30) {
            throw ValueError("Association implicit sensitivity matrix is singular.");
        }
        if (pivot != col) {
            for (int j = col; j < n; ++j) {
                std::swap(matrix[static_cast<size_t>(col * n + j)], matrix[static_cast<size_t>(pivot * n + j)]);
            }
            std::swap(rhs[static_cast<size_t>(col)], rhs[static_cast<size_t>(pivot)]);
        }
        const Scalar pivot_value = matrix[static_cast<size_t>(col * n + col)];
        for (int row = col + 1; row < n; ++row) {
            const Scalar factor = matrix[static_cast<size_t>(row * n + col)] / pivot_value;
            matrix[static_cast<size_t>(row * n + col)] = scalar_constant<Scalar>(0.0);
            for (int j = col + 1; j < n; ++j) {
                matrix[static_cast<size_t>(row * n + j)] -= factor * matrix[static_cast<size_t>(col * n + j)];
            }
            rhs[static_cast<size_t>(row)] -= factor * rhs[static_cast<size_t>(col)];
        }
    }

    vector<Scalar> out(static_cast<size_t>(n), scalar_constant<Scalar>(0.0));
    for (int row = n - 1; row >= 0; --row) {
        Scalar accum = rhs[static_cast<size_t>(row)];
        for (int col = row + 1; col < n; ++col) {
            accum -= matrix[static_cast<size_t>(row * n + col)] * out[static_cast<size_t>(col)];
        }
        out[static_cast<size_t>(row)] = accum / matrix[static_cast<size_t>(row * n + row)];
    }
    return out;
}

template <typename Scalar>
static vector<Scalar> association_site_fraction_density_terms_scalar_cpp(
    const vector<Scalar> &delta_ij,
    const Scalar &den,
    const vector<Scalar> &XA,
    const vector<Scalar> &ddelta_weighted,
    const vector<Scalar> &x_assoc
) {
    const int num_sites = static_cast<int>(XA.size());
    vector<Scalar> matrix(static_cast<size_t>(num_sites * num_sites), scalar_constant<Scalar>(0.0));
    vector<Scalar> rhs(static_cast<size_t>(num_sites), scalar_constant<Scalar>(0.0));

    int ij = 0;
    for (int i = 0; i < num_sites; ++i) {
        Scalar summ = scalar_constant<Scalar>(0.0);
        for (int j = 0; j < num_sites; ++j) {
            rhs[static_cast<size_t>(i)] -= x_assoc[static_cast<size_t>(j)] * XA[static_cast<size_t>(j)] * ddelta_weighted[static_cast<size_t>(ij)];
            matrix[static_cast<size_t>(i * num_sites + j)] = x_assoc[static_cast<size_t>(j)] * delta_ij[static_cast<size_t>(ij)];
            summ += x_assoc[static_cast<size_t>(j)] * XA[static_cast<size_t>(j)] * delta_ij[static_cast<size_t>(ij)];
            ++ij;
        }
        rhs[static_cast<size_t>(i)] -= summ;
        matrix[static_cast<size_t>(i * num_sites + i)] = scalar_pow(1.0 + den * summ, 2.0) / den;
    }

    return solve_linear_system_scalar_cpp(matrix, rhs, num_sites);
}

template <typename Scalar>
static vector<Scalar> association_site_fraction_composition_terms_scalar_cpp(
    const vector<Scalar> &delta_ij,
    const Scalar &den,
    const vector<Scalar> &XA,
    const vector<Scalar> &ddelta_dx,
    const vector<int> &site_component_index,
    const vector<Scalar> &x_assoc,
    int ncomp
) {
    const int num_sites = static_cast<int>(XA.size());
    vector<Scalar> dXA_dx(static_cast<size_t>(ncomp * num_sites), scalar_constant<Scalar>(0.0));

    for (int k = 0; k < ncomp; ++k) {
        vector<Scalar> matrix(static_cast<size_t>(num_sites * num_sites), scalar_constant<Scalar>(0.0));
        vector<Scalar> rhs(static_cast<size_t>(num_sites), scalar_constant<Scalar>(0.0));

        int ij = 0;
        for (int i = 0; i < num_sites; ++i) {
            Scalar direct_sum = scalar_constant<Scalar>(0.0);
            Scalar delta_sum = scalar_constant<Scalar>(0.0);
            for (int j = 0; j < num_sites; ++j) {
                if (site_component_index[static_cast<size_t>(j)] == k) {
                    direct_sum += XA[static_cast<size_t>(j)] * delta_ij[static_cast<size_t>(ij)];
                }
                delta_sum += x_assoc[static_cast<size_t>(j)] * XA[static_cast<size_t>(j)]
                    * ddelta_dx[static_cast<size_t>(k * num_sites * num_sites + ij)];
                matrix[static_cast<size_t>(i * num_sites + j)] = x_assoc[static_cast<size_t>(j)] * delta_ij[static_cast<size_t>(ij)];
                ++ij;
            }
            rhs[static_cast<size_t>(i)] = -(direct_sum + delta_sum);
            matrix[static_cast<size_t>(i * num_sites + i)] += 1.0 / (den * XA[static_cast<size_t>(i)] * XA[static_cast<size_t>(i)]);
        }

        vector<Scalar> solution = solve_linear_system_scalar_cpp(matrix, rhs, num_sites);
        for (int i = 0; i < num_sites; ++i) {
            dXA_dx[static_cast<size_t>(k * num_sites + i)] = solution[static_cast<size_t>(i)];
        }
    }

    return dXA_dx;
}

template <typename Scalar>
static AssociationImplicitTermsScalar<Scalar> association_implicit_terms_scalar_cpp(
    const MixtureStateScalar<Scalar> &thermo,
    const HardChainStateScalar<Scalar> &hc_state,
    double t,
    const vector<Scalar> &x,
    const add_args &cppargs,
    const vector<Scalar> &XA
) {
    const int ncomp = static_cast<int>(x.size());
    AssociationImplicitTermsScalar<Scalar> out;
    out.dadx.assign(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));
    out.mu.assign(static_cast<size_t>(ncomp), scalar_constant<Scalar>(0.0));

    for (int comp = 0; comp < static_cast<int>(cppargs.assoc_num.size()); ++comp) {
        for (int site = 0; site < cppargs.assoc_num[static_cast<size_t>(comp)]; ++site) {
            out.site_component_index.push_back(comp);
            out.x_assoc.push_back(x[static_cast<size_t>(comp)]);
        }
    }

    const int num_sites = static_cast<int>(out.site_component_index.size());
    out.delta_ij.assign(static_cast<size_t>(num_sites * num_sites), scalar_constant<Scalar>(0.0));
    out.residuals.assign(static_cast<size_t>(num_sites), scalar_constant<Scalar>(0.0));
    if (num_sites == 0) {
        return out;
    }
    if (static_cast<int>(XA.size()) != num_sites) {
        throw ValueError("Association implicit sensitivity received a site-fraction vector with the wrong size.");
    }

    int idxa = 0;
    for (int i = 0; i < num_sites; ++i) {
        const int comp_i = out.site_component_index[static_cast<size_t>(i)];
        for (int j = 0; j < num_sites; ++j) {
            const int comp_j = out.site_component_index[static_cast<size_t>(j)];
            if (cppargs.assoc_matrix[static_cast<size_t>(idxa)] != 0) {
                const double eABij = 0.5 * (cppargs.e_assoc[static_cast<size_t>(comp_i)] + cppargs.e_assoc[static_cast<size_t>(comp_j)]);
                const Scalar volABij = association_volume_scalar_cpp(comp_i, comp_j, ncomp, thermo.s_ij, cppargs);
                out.delta_ij[static_cast<size_t>(idxa)] =
                    hc_state.ghs[static_cast<size_t>(comp_i * ncomp + comp_j)]
                    * scalar_constant<Scalar>(std::exp(eABij / t) - 1.0)
                    * scalar_pow(thermo.s_ij[static_cast<size_t>(comp_i * ncomp + comp_j)], 3.0)
                    * volABij;
            }
            ++idxa;
        }
    }

    for (int i = 0; i < num_sites; ++i) {
        Scalar summ = scalar_constant<Scalar>(0.0);
        for (int j = 0; j < num_sites; ++j) {
            summ += out.x_assoc[static_cast<size_t>(j)] * XA[static_cast<size_t>(j)]
                * out.delta_ij[static_cast<size_t>(i * num_sites + j)];
        }
        out.residuals[static_cast<size_t>(i)] = XA[static_cast<size_t>(i)] * (1.0 + thermo.den * summ) - 1.0;
        const int component_index = out.site_component_index[static_cast<size_t>(i)];
        out.ares += x[static_cast<size_t>(component_index)] * (
            scalar_log(XA[static_cast<size_t>(i)]) - 0.5 * XA[static_cast<size_t>(i)] + 0.5
        );
    }

    vector<Scalar> ddelta_weighted(static_cast<size_t>(num_sites * num_sites), scalar_constant<Scalar>(0.0));
    idxa = 0;
    for (int i = 0; i < num_sites; ++i) {
        const int comp_i = out.site_component_index[static_cast<size_t>(i)];
        for (int j = 0; j < num_sites; ++j) {
            const int comp_j = out.site_component_index[static_cast<size_t>(j)];
            if (cppargs.assoc_matrix[static_cast<size_t>(idxa)] != 0) {
                const Scalar pair_diameter = thermo.d[static_cast<size_t>(comp_i)] * thermo.d[static_cast<size_t>(comp_j)]
                    / (thermo.d[static_cast<size_t>(comp_i)] + thermo.d[static_cast<size_t>(comp_j)]);
                const double eABij = 0.5 * (cppargs.e_assoc[static_cast<size_t>(comp_i)] + cppargs.e_assoc[static_cast<size_t>(comp_j)]);
                const Scalar volABij = association_volume_scalar_cpp(comp_i, comp_j, ncomp, thermo.s_ij, cppargs);
                ddelta_weighted[static_cast<size_t>(idxa)] =
                    hs_contact_density_derivative_scalar_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3])
                    * scalar_constant<Scalar>(std::exp(eABij / t) - 1.0)
                    * scalar_pow(thermo.s_ij[static_cast<size_t>(comp_i * ncomp + comp_j)], 3.0)
                    * volABij;
            }
            ++idxa;
        }
    }

    vector<Scalar> dXA_density = association_site_fraction_density_terms_scalar_cpp(
        out.delta_ij, thermo.den, XA, ddelta_weighted, out.x_assoc
    );
    for (int i = 0; i < num_sites; ++i) {
        const int component_index = out.site_component_index[static_cast<size_t>(i)];
        out.zraw += x[static_cast<size_t>(component_index)] * (1.0 / XA[static_cast<size_t>(i)] - 0.5)
            * dXA_density[static_cast<size_t>(i)];
    }

    vector<Scalar> ddelta_dx(static_cast<size_t>(num_sites * num_sites * ncomp), scalar_constant<Scalar>(0.0));
    int idx_ddelta = 0;
    for (int k = 0; k < ncomp; ++k) {
        for (int i = 0; i < num_sites; ++i) {
            const int comp_i = out.site_component_index[static_cast<size_t>(i)];
            for (int j = 0; j < num_sites; ++j) {
                const int comp_j = out.site_component_index[static_cast<size_t>(j)];
                if (cppargs.assoc_matrix[static_cast<size_t>(i * num_sites + j)] != 0) {
                    const Scalar pair_diameter = thermo.d[static_cast<size_t>(comp_i)] * thermo.d[static_cast<size_t>(comp_j)]
                        / (thermo.d[static_cast<size_t>(comp_i)] + thermo.d[static_cast<size_t>(comp_j)]);
                    const Scalar dzeta2_dx = scalar_constant<Scalar>(PI / 6.0 * cppargs.m[static_cast<size_t>(k)])
                        * thermo.den
                        * scalar_pow(thermo.d[static_cast<size_t>(k)], 2.0);
                    const Scalar dzeta3_dx = scalar_constant<Scalar>(PI / 6.0 * cppargs.m[static_cast<size_t>(k)])
                        * thermo.den
                        * scalar_pow(thermo.d[static_cast<size_t>(k)], 3.0);
                    const double eABij = 0.5 * (cppargs.e_assoc[static_cast<size_t>(comp_i)] + cppargs.e_assoc[static_cast<size_t>(comp_j)]);
                    const Scalar volABij = association_volume_scalar_cpp(comp_i, comp_j, ncomp, thermo.s_ij, cppargs);
                    ddelta_dx[static_cast<size_t>(idx_ddelta)] =
                        hs_contact_composition_derivative_scalar_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3], dzeta2_dx, dzeta3_dx)
                        * scalar_constant<Scalar>(std::exp(eABij / t) - 1.0)
                        * scalar_pow(thermo.s_ij[static_cast<size_t>(comp_i * ncomp + comp_j)], 3.0)
                        * volABij;
                }
                ++idx_ddelta;
            }
        }
    }

    vector<Scalar> dXA_dx = association_site_fraction_composition_terms_scalar_cpp(
        out.delta_ij,
        thermo.den,
        XA,
        ddelta_dx,
        out.site_component_index,
        out.x_assoc,
        ncomp
    );
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < num_sites; ++j) {
            out.dadx[static_cast<size_t>(i)] += x[static_cast<size_t>(out.site_component_index[static_cast<size_t>(j)])]
                * dXA_dx[static_cast<size_t>(i * num_sites + j)] * (1.0 / XA[static_cast<size_t>(j)] - 0.5);
        }
    }
    for (int i = 0; i < num_sites; ++i) {
        const int component_index = out.site_component_index[static_cast<size_t>(i)];
        out.dadx[static_cast<size_t>(component_index)] += scalar_log(XA[static_cast<size_t>(i)])
            - 0.5 * XA[static_cast<size_t>(i)] + 0.5;
    }
    for (int i = 0; i < ncomp; ++i) {
        out.sum_x_dadx += x[static_cast<size_t>(i)] * out.dadx[static_cast<size_t>(i)];
    }
    for (int i = 0; i < ncomp; ++i) {
        out.mu[static_cast<size_t>(i)] = out.ares + out.zraw + out.dadx[static_cast<size_t>(i)] - out.sum_x_dadx;
    }
    return out;
}

// EqID: ares_dh
template <typename Scalar, typename TemperatureScalar>
static Scalar ares_ion_scalar_cpp(
    const TemperatureScalar &t,
    const MixtureStateScalar<Scalar> &thermo,
    const vector<Scalar> &x,
    const add_args &cppargs,
    int component_target_kind = -1,
    int component_target_index = -1,
    const Scalar *component_target_value = nullptr
) {
    if (cppargs.z.empty()) {
        return scalar_constant<Scalar>(0.0);
    }
    bool has_charge = false;
    for (double charge : cppargs.z) {
        if (std::abs(charge) > 1.0e-12) {
            has_charge = true;
            break;
        }
    }
    if (!has_charge) {
        return scalar_constant<Scalar>(0.0);
    }
    Scalar q2_sum = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        q2_sum += x[i] * cppargs.z[i] * cppargs.z[i];
    }
    Scalar eps = scalar_constant<Scalar>(0.0);
    if (cppargs.dielc_rule == 0) {
        if (component_target_kind == kGenericTargetDielcLocal) {
            throw ValueError("unsupported: CppAD ionic dielectric-parameter recording requires dielc_rule=1.");
        }
        eps = scalar_constant<Scalar>(*std::max_element(cppargs.dielc.begin(), cppargs.dielc.end()));
    } else if (cppargs.dielc_rule == 1) {
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            eps += x[i] * scalar_dielc_parameter_cpp(
                component_target_kind,
                component_target_index,
                i,
                cppargs.dielc[i],
                component_target_value
            );
        }
    } else {
        if (component_target_kind == kGenericTargetDielcLocal) {
            throw ValueError("unsupported: CppAD ionic dielectric-parameter recording requires dielc_rule=1.");
        }
#ifdef EPCSAFT_HAS_CPPAD
        eps = dielectric_constant_rule_cppad_cpp(cppargs.dielc_rule, x, cppargs);
#else
        throw ValueError("CppAD ionic recording requires CppAD support for composition-dependent dielectric rules.");
#endif
    }
    Scalar kappa = scalar_sqrt(thermo.den * E_CHRG * E_CHRG / kb / t / (eps * perm_vac) * q2_sum);
    Scalar chi_sum = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        Scalar ka = kappa * thermo.d[i];
        Scalar chi = 3.0 / scalar_pow(ka, 3) * (1.5 + scalar_log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * scalar_pow(1.0 + ka, 2));
        chi_sum += x[i] * cppargs.z[i] * cppargs.z[i] * chi;
    }
    Scalar K0 = scalar_constant<Scalar>(E_CHRG * E_CHRG / (12.0 * PI * kb * perm_vac)) / t;
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
template <typename Scalar, typename TemperatureScalar>
static Scalar ares_born_scalar_cpp(
    const TemperatureScalar &t,
    const vector<Scalar> &x,
    const add_args &cppargs,
    int component_target_kind = -1,
    int component_target_index = -1,
    const Scalar *component_target_value = nullptr
) {
    if (cppargs.z.empty() || cppargs.born_model == 0) {
        return scalar_constant<Scalar>(0.0);
    }
    bool has_charge = false;
    for (double charge : cppargs.z) {
        if (std::abs(charge) > 1.0e-12) {
            has_charge = true;
            break;
        }
    }
    if (!has_charge) {
        return scalar_constant<Scalar>(0.0);
    }
    if (cppargs.born_model != 1) {
        if (cppargs.born_model != 2) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }
    Scalar eps = scalar_constant<Scalar>(0.0);
    if (cppargs.dielc_rule == 0) {
        if (component_target_kind == kGenericTargetDielcLocal) {
            throw ValueError("unsupported: CppAD Born dielectric-parameter recording requires dielc_rule=1.");
        }
        eps = scalar_constant<Scalar>(*std::max_element(cppargs.dielc.begin(), cppargs.dielc.end()));
    } else if (cppargs.dielc_rule == 1) {
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            eps += x[i] * scalar_dielc_parameter_cpp(
                component_target_kind,
                component_target_index,
                i,
                cppargs.dielc[i],
                component_target_value
            );
        }
    } else {
        if (component_target_kind == kGenericTargetDielcLocal) {
            throw ValueError("unsupported: CppAD Born dielectric-parameter recording requires dielc_rule=1.");
        }
#ifdef EPCSAFT_HAS_CPPAD
        eps = dielectric_constant_rule_cppad_cpp(cppargs.dielc_rule, x, cppargs);
#else
        throw ValueError("CppAD Born recording requires CppAD support for composition-dependent dielectric rules.");
#endif
    }
    if (cppargs.born_eps_mode == 1) {
        throw ValueError("unsupported: CppAD Born reference-solvent dielectric routing requires a tape route.");
    }
    if (cppargs.born_model == 1) {
        Scalar charge_radius_sum = scalar_constant<Scalar>(0.0);
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            if (!is_ion_species(cppargs, i)) {
                continue;
            }
            Scalar sigma_i = scalar_component_parameter_cpp(
                component_target_kind,
                component_target_index,
                i,
                cppargs.s[i],
                component_target_kind == kGenericTargetSLocal ? component_target_value : nullptr
            );
            Scalar epsilon_i = scalar_component_parameter_cpp(
                component_target_kind,
                component_target_index,
                i,
                cppargs.e[i],
                component_target_kind == kGenericTargetELocal ? component_target_value : nullptr
            );
            Scalar d_born_i = ion_born_radius_scalar_cpp(
                i,
                t,
                cppargs,
                sigma_i,
                epsilon_i,
                component_target_kind,
                component_target_index,
                component_target_value
            );
            charge_radius_sum += x[i] * cppargs.z[i] * cppargs.z[i] / d_born_i;
        }
        Scalar Kborn = scalar_constant<Scalar>(E_CHRG * E_CHRG / (4.0 * PI * kb * perm_vac)) / t;
        return -Kborn * (1.0 - 1.0 / eps) * charge_radius_sum;
    }

    const bool use_ssm = (cppargs.born_solvation_shell_model != 0);
    const bool use_ds = (cppargs.born_dielectric_saturation != 0);
    const Scalar eps_ion = scalar_constant<Scalar>(8.0);
    Scalar f_mix = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        Scalar f_i = scalar_constant<Scalar>(1.0);
        if (!is_ion_species(cppargs, i) && cppargs.f_solv.size() > static_cast<size_t>(i)) {
            f_i = scalar_component_parameter_cpp(
                component_target_kind,
                component_target_index,
                i,
                cppargs.f_solv[i],
                component_target_kind == kGenericTargetFSolvLocal ? component_target_value : nullptr
            );
        }
        f_mix += x[i] * f_i;
    }

    Scalar sum_bracket = scalar_constant<Scalar>(0.0);
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        if (!is_ion_species(cppargs, i)) {
            continue;
        }
        Scalar sigma_i = scalar_component_parameter_cpp(
            component_target_kind,
            component_target_index,
            i,
            cppargs.s[i],
            component_target_kind == kGenericTargetSLocal ? component_target_value : nullptr
        );
        Scalar epsilon_i = scalar_component_parameter_cpp(
            component_target_kind,
            component_target_index,
            i,
            cppargs.e[i],
            component_target_kind == kGenericTargetELocal ? component_target_value : nullptr
        );
        Scalar d_born_i = ion_born_radius_scalar_cpp(
            i,
            t,
            cppargs,
            sigma_i,
            epsilon_i,
            component_target_kind,
            component_target_index,
            component_target_value
        );
        const double absz = std::abs(cppargs.z[i]);
        const double z2 = cppargs.z[i] * cppargs.z[i];
        Scalar delta_di = use_ssm ? ((f_mix - 1.0) * d_born_i / absz) : scalar_constant<Scalar>(0.0);
        Scalar D_i = d_born_i + delta_di;
        if (scalar_value(D_i) <= 0.0) {
            throw ValueError("Born model generated a non-positive d_born + Delta d.");
        }
        Scalar invD = 1.0 / D_i;
        Scalar gap = 1.0 / d_born_i - invD;
        Scalar base_term = (1.0 - 1.0 / eps) * invD;
        Scalar ds_term = use_ds ? ((1.0 - 1.0 / eps_ion) * gap) : scalar_constant<Scalar>(0.0);
        sum_bracket += x[i] * z2 * (base_term + ds_term);
    }
    Scalar Kborn = scalar_constant<Scalar>(E_CHRG * E_CHRG / (4.0 * PI * kb * perm_vac)) / t;
    return -Kborn * sum_bracket;
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

template <typename Scalar, typename TemperatureScalar>
static AresContributionsScalar<Scalar> ares_contributions_scalar_cpp(
    const TemperatureScalar &t,
    const Scalar &rho,
    const vector<Scalar> &x,
    const add_args &cppargs,
    int k_override_index = -1,
    const Scalar *k_override_value = nullptr,
    int l_override_index = -1,
    const Scalar *l_override_value = nullptr,
    int component_target_kind = -1,
    int component_target_index = -1,
    const Scalar *component_target_value = nullptr
) {
    AresContributionsScalar<Scalar> out;
    MixtureStateScalar<Scalar> thermo = mixture_state_scalar_cpp(
        t,
        rho,
        x,
        cppargs,
        k_override_index,
        k_override_value,
        l_override_index,
        l_override_value,
        component_target_kind,
        component_target_index,
        component_target_value
    );
    HardChainStateScalar<Scalar> hc_state = hard_chain_state_scalar_cpp(thermo, x, cppargs);
    DispersionPolynomialStateScalar<Scalar> dispersion = dispersion_polynomials_scalar_cpp(thermo.m_avg, hc_state.eta);
    out.hc = ares_hc_scalar_cpp(thermo, hc_state, x, cppargs);
    out.disp = ares_disp_scalar_cpp(thermo, dispersion);
    out.assoc = ares_assoc_scalar_cpp(x, cppargs);
    out.ion = ares_ion_scalar_cpp(
        t,
        thermo,
        x,
        cppargs,
        component_target_kind,
        component_target_index,
        component_target_value
    );
    out.born = ares_born_scalar_cpp(t, x, cppargs, component_target_kind, component_target_index, component_target_value);
    return out;
}

}  // namespace ares_detail

namespace {

std::string contribution_backend_name(int mode) {
    if (mode == 0) return "analytic";
    if (mode == 1) return "unsupported";
    if (mode == 2) return "cppad";
    if (mode == 3 || mode == 5) {
        return "analytic";
    }
    if (mode == 4) return "cppad";
    return "unknown";
}

bool has_association_sites(const add_args &cppargs) {
    for (int sites : cppargs.assoc_num) {
        if (sites > 0) {
            return true;
        }
    }
    return false;
}

std::string association_backend_name(const add_args &cppargs) {
    if (!has_association_sites(cppargs)) {
        return contribution_backend_name(cppargs.assoc_dadx_diff_mode);
    }
    if (cppargs.assoc_dadx_diff_mode == 0 || cppargs.assoc_dadx_diff_mode == 3 || cppargs.assoc_dadx_diff_mode == 5) {
        return "analytic_implicit";
    }
    return "unsupported";
}

std::map<std::string, std::string> composition_derivative_backend_map(const add_args &cppargs) {
    std::map<std::string, std::string> backends;
    backends["hc"] = contribution_backend_name(cppargs.hc_dadx_diff_mode);
    backends["disp"] = contribution_backend_name(cppargs.disp_dadx_diff_mode);
    backends["assoc"] = association_backend_name(cppargs);
    backends["ion"] = contribution_backend_name(cppargs.mu_DH_diff_mode);
    backends["born"] = contribution_backend_name(cppargs.born_diff_mode);
    return backends;
}

#ifdef EPCSAFT_HAS_CPPAD
ares_detail::AssociationDensityResponse association_density_response_cppad_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs
) {
    using CppADScalar = CppAD::AD<double>;
    const int ncomp = static_cast<int>(x.size());
    MixtureState thermo = mixture_state_cpp(t, rho, x, cppargs, false);
    HardChainState hc_state = hard_chain_state_cpp(thermo, x, cppargs);
    AssociationIntermediateState assoc_state = association_intermediate_state_cpp(thermo, hc_state, t, x, cppargs, false, false);

    ares_detail::AssociationDensityResponse out;
    out.mu.assign(static_cast<size_t>(ncomp), 0.0);
    out.dmu_drho.assign(static_cast<size_t>(ncomp), 0.0);
    if (!assoc_state.active) {
        return out;
    }

    const int num_sites = static_cast<int>(assoc_state.XA.size());
    const int var_count = 1 + num_sites;
    std::vector<CppADScalar> avars(static_cast<size_t>(var_count));
    avars[0] = rho;
    for (int i = 0; i < num_sites; ++i) {
        avars[static_cast<size_t>(1 + i)] = assoc_state.XA[static_cast<size_t>(i)];
    }
    CppAD::Independent(avars);

    std::vector<CppADScalar> ax(static_cast<size_t>(ncomp));
    for (int i = 0; i < ncomp; ++i) {
        ax[static_cast<size_t>(i)] = CppADScalar(x[static_cast<size_t>(i)]);
    }

    ares_detail::MixtureStateScalar<CppADScalar> scalar_thermo = ares_detail::mixture_state_scalar_cpp(
        t,
        avars[0],
        ax,
        cppargs
    );
    ares_detail::HardChainStateScalar<CppADScalar> scalar_hc = ares_detail::hard_chain_state_scalar_cpp(
        scalar_thermo,
        ax,
        cppargs
    );
    std::vector<CppADScalar> site_vars(static_cast<size_t>(num_sites));
    for (int i = 0; i < num_sites; ++i) {
        site_vars[static_cast<size_t>(i)] = avars[static_cast<size_t>(1 + i)];
    }
    auto assoc_terms = ares_detail::association_implicit_terms_scalar_cpp(
        scalar_thermo,
        scalar_hc,
        t,
        ax,
        cppargs,
        site_vars
    );

    std::vector<CppADScalar> ay;
    ay.reserve(static_cast<size_t>(1 + ncomp + num_sites));
    ay.push_back(assoc_terms.zraw);
    for (int i = 0; i < ncomp; ++i) {
        ay.push_back(assoc_terms.mu[static_cast<size_t>(i)]);
    }
    for (int i = 0; i < num_sites; ++i) {
        ay.push_back(assoc_terms.residuals[static_cast<size_t>(i)]);
    }

    CppAD::ADFun<double> function(avars, ay);
    std::vector<double> point(static_cast<size_t>(var_count), 0.0);
    point[0] = rho;
    for (int i = 0; i < num_sites; ++i) {
        point[static_cast<size_t>(1 + i)] = assoc_state.XA[static_cast<size_t>(i)];
    }
    auto values = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);
    const auto jac = [&](int row, int col) {
        return jacobian[static_cast<size_t>(row * var_count + col)];
    };

    vector<double> residual_matrix(static_cast<size_t>(num_sites * num_sites), 0.0);
    vector<double> residual_rhs(static_cast<size_t>(num_sites), 0.0);
    const int residual_row0 = 1 + ncomp;
    for (int i = 0; i < num_sites; ++i) {
        residual_rhs[static_cast<size_t>(i)] = -jac(residual_row0 + i, 0);
        for (int j = 0; j < num_sites; ++j) {
            residual_matrix[static_cast<size_t>(i * num_sites + j)] = jac(residual_row0 + i, 1 + j);
        }
    }
    vector<double> dXA_drho = ares_detail::solve_linear_system_scalar_cpp(residual_matrix, residual_rhs, num_sites);

    out.zraw = values[0];
    out.dzraw_drho = jac(0, 0);
    for (int j = 0; j < num_sites; ++j) {
        out.dzraw_drho += jac(0, 1 + j) * dXA_drho[static_cast<size_t>(j)];
    }
    for (int i = 0; i < ncomp; ++i) {
        const int output_row = 1 + i;
        out.mu[static_cast<size_t>(i)] = values[static_cast<size_t>(output_row)];
        out.dmu_drho[static_cast<size_t>(i)] = jac(output_row, 0);
        for (int j = 0; j < num_sites; ++j) {
            out.dmu_drho[static_cast<size_t>(i)] += jac(output_row, 1 + j) * dXA_drho[static_cast<size_t>(j)];
        }
    }
    return out;
}

ares_detail::AssociationPhaseStateResponse association_phase_state_response_cppad_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs
) {
    using CppADScalar = CppAD::AD<double>;
    const int ncomp = static_cast<int>(x.size());
    const int base_var_count = 1 + ncomp;
    MixtureState thermo = mixture_state_cpp(t, rho, x, cppargs, false);
    HardChainState hc_state = hard_chain_state_cpp(thermo, x, cppargs);
    AssociationIntermediateState assoc_state = association_intermediate_state_cpp(thermo, hc_state, t, x, cppargs, false, false);

    ares_detail::AssociationPhaseStateResponse out;
    out.ncomp = ncomp;
    out.base_var_count = base_var_count;
    out.mu.assign(static_cast<size_t>(ncomp), 0.0);
    out.dzraw_dvar.assign(static_cast<size_t>(base_var_count), 0.0);
    out.dmu_dvar_row_major.assign(static_cast<size_t>(ncomp * base_var_count), 0.0);
    out.d2zraw_dvar2_row_major.assign(static_cast<size_t>(base_var_count * base_var_count), 0.0);
    out.d2mu_dvar2_tensor_row_major.assign(
        static_cast<size_t>(ncomp * base_var_count * base_var_count),
        0.0
    );
    if (!assoc_state.active) {
        return out;
    }
    out.active = true;

    const int num_sites = static_cast<int>(assoc_state.XA.size());
    const int var_count = base_var_count + num_sites;
    std::vector<CppADScalar> avars(static_cast<size_t>(var_count));
    avars[0] = rho;
    for (int i = 0; i < ncomp; ++i) {
        avars[static_cast<size_t>(1 + i)] = x[static_cast<size_t>(i)];
    }
    for (int i = 0; i < num_sites; ++i) {
        avars[static_cast<size_t>(base_var_count + i)] = assoc_state.XA[static_cast<size_t>(i)];
    }
    CppAD::Independent(avars);

    std::vector<CppADScalar> ax(static_cast<size_t>(ncomp));
    for (int i = 0; i < ncomp; ++i) {
        ax[static_cast<size_t>(i)] = avars[static_cast<size_t>(1 + i)];
    }
    ares_detail::MixtureStateScalar<CppADScalar> scalar_thermo = ares_detail::mixture_state_scalar_cpp(
        t,
        avars[0],
        ax,
        cppargs
    );
    ares_detail::HardChainStateScalar<CppADScalar> scalar_hc = ares_detail::hard_chain_state_scalar_cpp(
        scalar_thermo,
        ax,
        cppargs
    );
    std::vector<CppADScalar> site_vars(static_cast<size_t>(num_sites));
    for (int i = 0; i < num_sites; ++i) {
        site_vars[static_cast<size_t>(i)] = avars[static_cast<size_t>(base_var_count + i)];
    }
    auto assoc_terms = ares_detail::association_implicit_terms_scalar_cpp(
        scalar_thermo,
        scalar_hc,
        t,
        ax,
        cppargs,
        site_vars
    );

    std::vector<CppADScalar> ay;
    ay.reserve(static_cast<size_t>(1 + ncomp + num_sites));
    ay.push_back(assoc_terms.zraw);
    for (int i = 0; i < ncomp; ++i) {
        ay.push_back(assoc_terms.mu[static_cast<size_t>(i)]);
    }
    for (int i = 0; i < num_sites; ++i) {
        ay.push_back(assoc_terms.residuals[static_cast<size_t>(i)]);
    }

    CppAD::ADFun<double> function(avars, ay);
    std::vector<double> point(static_cast<size_t>(var_count), 0.0);
    point[0] = rho;
    for (int i = 0; i < ncomp; ++i) {
        point[static_cast<size_t>(1 + i)] = x[static_cast<size_t>(i)];
    }
    for (int i = 0; i < num_sites; ++i) {
        point[static_cast<size_t>(base_var_count + i)] = assoc_state.XA[static_cast<size_t>(i)];
    }
    auto values = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);
    std::vector<std::vector<double>> output_hessians;
    output_hessians.reserve(static_cast<std::size_t>(1 + ncomp + num_sites));
    for (int output = 0; output < 1 + ncomp + num_sites; ++output) {
        output_hessians.push_back(ares_detail::vector_output_component_hessian_cppad(
            function,
            point,
            static_cast<std::size_t>(output)
        ));
    }
    const auto jac = [&](int row, int col) {
        return jacobian[static_cast<size_t>(row * var_count + col)];
    };
    const auto hess = [&](int output, int row, int col) {
        return output_hessians[static_cast<std::size_t>(output)][
            static_cast<std::size_t>(row * var_count + col)
        ];
    };

    const int residual_row0 = 1 + ncomp;
    vector<double> residual_matrix(static_cast<size_t>(num_sites * num_sites), 0.0);
    for (int i = 0; i < num_sites; ++i) {
        for (int j = 0; j < num_sites; ++j) {
            residual_matrix[static_cast<size_t>(i * num_sites + j)] =
                jac(residual_row0 + i, base_var_count + j);
        }
    }

    out.zraw = values[0];
    for (int i = 0; i < ncomp; ++i) {
        out.mu[static_cast<size_t>(i)] = values[static_cast<size_t>(1 + i)];
    }
    std::vector<double> dxa_dvar(static_cast<size_t>(base_var_count * num_sites), 0.0);
    for (int v = 0; v < base_var_count; ++v) {
        vector<double> rhs(static_cast<size_t>(num_sites), 0.0);
        for (int i = 0; i < num_sites; ++i) {
            rhs[static_cast<size_t>(i)] = -jac(residual_row0 + i, v);
        }
        vector<double> dxa_dv = ares_detail::solve_linear_system_scalar_cpp(residual_matrix, rhs, num_sites);
        for (int site = 0; site < num_sites; ++site) {
            dxa_dvar[static_cast<size_t>(v * num_sites + site)] = dxa_dv[static_cast<size_t>(site)];
        }
        out.dzraw_dvar[static_cast<size_t>(v)] = jac(0, v);
        for (int site = 0; site < num_sites; ++site) {
            out.dzraw_dvar[static_cast<size_t>(v)] += jac(0, base_var_count + site) * dxa_dv[static_cast<size_t>(site)];
        }
        for (int i = 0; i < ncomp; ++i) {
            double value = jac(1 + i, v);
            for (int site = 0; site < num_sites; ++site) {
                value += jac(1 + i, base_var_count + site) * dxa_dv[static_cast<size_t>(site)];
            }
            out.dmu_dvar_row_major[static_cast<size_t>(i * base_var_count + v)] = value;
        }
    }
    for (int u = 0; u < base_var_count; ++u) {
        for (int v = 0; v < base_var_count; ++v) {
            vector<double> rhs(static_cast<size_t>(num_sites), 0.0);
            for (int residual = 0; residual < num_sites; ++residual) {
                double value = -hess(residual_row0 + residual, u, v);
                for (int site = 0; site < num_sites; ++site) {
                    const double dy_v = dxa_dvar[static_cast<size_t>(v * num_sites + site)];
                    const double dy_u = dxa_dvar[static_cast<size_t>(u * num_sites + site)];
                    value -= hess(residual_row0 + residual, u, base_var_count + site) * dy_v;
                    value -= hess(residual_row0 + residual, v, base_var_count + site) * dy_u;
                }
                for (int site_i = 0; site_i < num_sites; ++site_i) {
                    const double dy_u = dxa_dvar[static_cast<size_t>(u * num_sites + site_i)];
                    if (dy_u == 0.0) {
                        continue;
                    }
                    for (int site_j = 0; site_j < num_sites; ++site_j) {
                        value -= hess(
                            residual_row0 + residual,
                            base_var_count + site_i,
                            base_var_count + site_j
                        ) * dy_u * dxa_dvar[static_cast<size_t>(v * num_sites + site_j)];
                    }
                }
                rhs[static_cast<size_t>(residual)] = value;
            }
            const vector<double> d2xa = ares_detail::solve_linear_system_scalar_cpp(residual_matrix, rhs, num_sites);
            out.d2zraw_dvar2_row_major[static_cast<size_t>(u * base_var_count + v)] = hess(0, u, v);
            for (int site = 0; site < num_sites; ++site) {
                const double dy_v = dxa_dvar[static_cast<size_t>(v * num_sites + site)];
                const double dy_u = dxa_dvar[static_cast<size_t>(u * num_sites + site)];
                out.d2zraw_dvar2_row_major[static_cast<size_t>(u * base_var_count + v)] +=
                    hess(0, u, base_var_count + site) * dy_v
                    + hess(0, v, base_var_count + site) * dy_u
                    + jac(0, base_var_count + site) * d2xa[static_cast<size_t>(site)];
            }
            for (int site_i = 0; site_i < num_sites; ++site_i) {
                const double dy_u = dxa_dvar[static_cast<size_t>(u * num_sites + site_i)];
                if (dy_u == 0.0) {
                    continue;
                }
                for (int site_j = 0; site_j < num_sites; ++site_j) {
                    out.d2zraw_dvar2_row_major[static_cast<size_t>(u * base_var_count + v)] +=
                        hess(0, base_var_count + site_i, base_var_count + site_j)
                        * dy_u
                        * dxa_dvar[static_cast<size_t>(v * num_sites + site_j)];
                }
            }
            for (int component = 0; component < ncomp; ++component) {
                double value = hess(1 + component, u, v);
                for (int site = 0; site < num_sites; ++site) {
                    const double dy_v = dxa_dvar[static_cast<size_t>(v * num_sites + site)];
                    const double dy_u = dxa_dvar[static_cast<size_t>(u * num_sites + site)];
                    value += hess(1 + component, u, base_var_count + site) * dy_v;
                    value += hess(1 + component, v, base_var_count + site) * dy_u;
                    value += jac(1 + component, base_var_count + site) * d2xa[static_cast<size_t>(site)];
                }
                for (int site_i = 0; site_i < num_sites; ++site_i) {
                    const double dy_u = dxa_dvar[static_cast<size_t>(u * num_sites + site_i)];
                    if (dy_u == 0.0) {
                        continue;
                    }
                    for (int site_j = 0; site_j < num_sites; ++site_j) {
                        value += hess(
                            1 + component,
                            base_var_count + site_i,
                            base_var_count + site_j
                        ) * dy_u * dxa_dvar[static_cast<size_t>(v * num_sites + site_j)];
                    }
                }
                out.d2mu_dvar2_tensor_row_major[
                    static_cast<size_t>(component * base_var_count * base_var_count + u * base_var_count + v)
                ] = value;
            }
        }
    }
    return out;
}
#endif

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
    BornIntermediateState born_state = born_intermediate_state_cpp(t, x, cppargs, false);

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

epcsaft::native::cppad_support::CppADDerivativeResult cppad_eos_contribution_derivatives_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs
) {
    using CppADScalar = CppAD::AD<double>;
    if (!cppargs.assoc_num.empty()) {
        for (int sites : cppargs.assoc_num) {
            if (sites > 0) {
                throw ValueError("unsupported: CppAD association recording requires implicit site-fraction sensitivities.");
            }
        }
    }
    if (!cppargs.z.empty() && cppargs.born_model > 1) {
        throw ValueError("unsupported: CppAD Born recording supports direct Born model=1 formulas only.");
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

    epcsaft::native::cppad_support::CppADDerivativeResult result;
    result.supported = true;
    result.backend = "cppad";
    result.message = "CppAD EOS contribution composition derivatives available";
    result.value = std::move(value);
    result.jacobian_row_major = std::move(jacobian);
    result.outputs = {"hc", "disp", "assoc", "ion", "born", "total"};
    result.variables.reserve(static_cast<size_t>(ncomp));
    for (int i = 0; i < ncomp; ++i) {
        result.variables.push_back("x_" + std::to_string(i));
    }
    result.rows = 6;
    result.cols = ncomp;
    return result;
}

void eos_phase_objective_derivatives_cpp(
    double t,
    double target_pressure,
    const vector<double> &amounts,
    double volume,
    const add_args &cppargs,
    double *objective,
    vector<double> *gradient,
    vector<double> *hessian_row_major
) {
#ifdef EPCSAFT_HAS_CPPAD
    using CppADScalar = CppAD::AD<double>;
    if (!objective || !gradient || !hessian_row_major) {
        throw ValueError("EOS phase objective derivative output buffers must be valid.");
    }
    if (!cppargs.assoc_num.empty()) {
        for (int sites : cppargs.assoc_num) {
            if (sites > 0) {
                throw ValueError("unsupported: EOS phase objective Hessian requires association implicit variables.");
            }
        }
    }
    if (!cppargs.z.empty() && cppargs.born_model > 1) {
        throw ValueError("unsupported: EOS phase objective Hessian supports direct Born model=1 formulas only.");
    }
    const int ncomp = static_cast<int>(amounts.size());
    const int nvars = ncomp + 1;
    std::vector<CppADScalar> variables(static_cast<std::size_t>(nvars));
    for (int index = 0; index < ncomp; ++index) {
        variables[static_cast<std::size_t>(index)] = amounts[static_cast<std::size_t>(index)];
    }
    variables[static_cast<std::size_t>(ncomp)] = volume;
    CppAD::Independent(variables);

    CppADScalar total_amount = scalar_constant<CppADScalar>(0.0);
    for (int index = 0; index < ncomp; ++index) {
        total_amount += variables[static_cast<std::size_t>(index)];
    }
    const CppADScalar phase_volume = variables[static_cast<std::size_t>(ncomp)];
    const CppADScalar rho = total_amount / phase_volume;
    std::vector<CppADScalar> composition(static_cast<std::size_t>(ncomp));
    for (int index = 0; index < ncomp; ++index) {
        composition[static_cast<std::size_t>(index)] = variables[static_cast<std::size_t>(index)] / total_amount;
    }

    CppADScalar ideal = scalar_constant<CppADScalar>(0.0);
    for (int index = 0; index < ncomp; ++index) {
        const CppADScalar amount = variables[static_cast<std::size_t>(index)];
        ideal += amount * (CppAD::log(amount / phase_volume) - scalar_constant<CppADScalar>(1.0));
    }
    const auto contributions = ares_detail::ares_contributions_scalar_cpp(t, rho, composition, cppargs);
    const CppADScalar residual = total_amount
        * (contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born);
    const CppADScalar pressure_work = target_pressure * phase_volume / (kb * N_AV * t);
    std::vector<CppADScalar> outputs(1);
    outputs[0] = ideal + residual + pressure_work;

    CppAD::ADFun<double> function(variables, outputs);
    std::vector<double> point(amounts.begin(), amounts.end());
    point.push_back(volume);
    auto value = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);
    auto hessian = function.Hessian(point, 0);

    *objective = value.at(0);
    *gradient = std::move(jacobian);
    *hessian_row_major = std::move(hessian);
#else
    (void)t;
    (void)target_pressure;
    (void)amounts;
    (void)volume;
    (void)cppargs;
    (void)objective;
    (void)gradient;
    (void)hessian_row_major;
    throw ValueError("EOS phase objective Hessian requires CppAD support.");
#endif
}

void eos_phase_temperature_variable_derivatives_cpp(
    double t,
    double target_pressure,
    const vector<double> &amounts,
    double volume,
    const add_args &cppargs,
    double *objective,
    vector<double> *gradient,
    vector<double> *hessian_row_major
) {
#ifdef EPCSAFT_HAS_CPPAD
    using CppADScalar = CppAD::AD<double>;
    if (!objective || !gradient || !hessian_row_major) {
        throw ValueError("EOS phase variable-temperature derivative output buffers must be valid.");
    }
    if (!cppargs.assoc_num.empty()) {
        for (int sites : cppargs.assoc_num) {
            if (sites > 0) {
                throw ValueError("unsupported: EOS phase temperature route Hessian requires association implicit variables.");
            }
        }
    }
    if (!cppargs.z.empty() && cppargs.born_model > 1) {
        throw ValueError("unsupported: EOS phase temperature route Hessian supports direct Born model=1 formulas only.");
    }
    const int ncomp = static_cast<int>(amounts.size());
    const int nvars = ncomp + 2;
    std::vector<CppADScalar> variables(static_cast<std::size_t>(nvars));
    for (int index = 0; index < ncomp; ++index) {
        variables[static_cast<std::size_t>(index)] = amounts[static_cast<std::size_t>(index)];
    }
    variables[static_cast<std::size_t>(ncomp)] = volume;
    variables[static_cast<std::size_t>(ncomp + 1)] = t;
    CppAD::Independent(variables);

    CppADScalar total_amount = scalar_constant<CppADScalar>(0.0);
    for (int index = 0; index < ncomp; ++index) {
        total_amount += variables[static_cast<std::size_t>(index)];
    }
    const CppADScalar phase_volume = variables[static_cast<std::size_t>(ncomp)];
    const CppADScalar temperature = variables[static_cast<std::size_t>(ncomp + 1)];
    const CppADScalar rho = total_amount / phase_volume;
    std::vector<CppADScalar> composition(static_cast<std::size_t>(ncomp));
    for (int index = 0; index < ncomp; ++index) {
        composition[static_cast<std::size_t>(index)] = variables[static_cast<std::size_t>(index)] / total_amount;
    }

    CppADScalar ideal = scalar_constant<CppADScalar>(0.0);
    for (int index = 0; index < ncomp; ++index) {
        const CppADScalar amount = variables[static_cast<std::size_t>(index)];
        ideal += amount * (CppAD::log(amount / phase_volume) - scalar_constant<CppADScalar>(1.0));
    }
    const auto contributions = ares_detail::ares_contributions_scalar_cpp(temperature, rho, composition, cppargs);
    const CppADScalar residual = total_amount
        * (contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born);
    const CppADScalar pressure_work = target_pressure * phase_volume / (kb * N_AV * temperature);
    std::vector<CppADScalar> outputs(1);
    outputs[0] = ideal + residual + pressure_work;

    CppAD::ADFun<double> function(variables, outputs);
    std::vector<double> point(amounts.begin(), amounts.end());
    point.push_back(volume);
    point.push_back(t);
    auto value = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);
    auto hessian = function.Hessian(point, 0);

    *objective = value.at(0);
    *gradient = std::move(jacobian);
    *hessian_row_major = std::move(hessian);
#else
    (void)t;
    (void)target_pressure;
    (void)amounts;
    (void)volume;
    (void)cppargs;
    (void)objective;
    (void)gradient;
    (void)hessian_row_major;
    throw ValueError("EOS phase variable-temperature Hessian requires CppAD support.");
#endif
}

PhaseStateCompositionSensitivityResult phase_state_ln_fugacity_composition_sensitivity_cpp(
    double t,
    double p,
    vector<double> x,
    int phase,
    const add_args &cppargs
) {
    PhaseStateCompositionSensitivityResult out;
    out.temperature = t;
    out.pressure = p;
    out.composition = x;
    out.rows = static_cast<int>(x.size());
    out.cols = static_cast<int>(x.size());

    DensitySolveResult density = density_solve_report_cpp(t, p, x, phase, cppargs);
    if (!density.valid || !(density.rho > 0.0) || !std::isfinite(density.rho)) {
        out.message = density.diagnostics.rejection_reason.empty()
            ? "Density root was not valid for phase-state sensitivity evaluation."
            : density.diagnostics.rejection_reason;
        return out;
    }
    out.density = density.rho;
    out.ln_fugacity = fugacity_coefficient_result_cpp(t, density.rho, x, cppargs).lnfugcoef.total;

    using CppADScalar = CppAD::AD<double>;
    const int ncomp = static_cast<int>(x.size());
    const int var_count = ncomp + 1;
    std::vector<CppADScalar> avars(static_cast<size_t>(var_count));
    avars[0] = density.rho;
    for (int i = 0; i < ncomp; ++i) {
        avars[static_cast<size_t>(1 + i)] = x[static_cast<size_t>(i)];
    }
    CppAD::Independent(avars);

    std::vector<CppADScalar> ax(static_cast<size_t>(ncomp));
    for (int i = 0; i < ncomp; ++i) {
        ax[static_cast<size_t>(i)] = avars[static_cast<size_t>(1 + i)];
    }
    const bool active_association = has_association_sites(cppargs);
    add_args recording_args = cppargs;
    if (active_association) {
        recording_args.assoc_num.clear();
        recording_args.assoc_matrix.clear();
        recording_args.e_assoc.clear();
        recording_args.vol_a.clear();
        recording_args.k_hb.clear();
    }
    auto contributions = ares_detail::ares_contributions_scalar_cpp(t, avars[0], ax, recording_args);
    std::vector<CppADScalar> ay(1);
    ay[0] = contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born;

    CppAD::ADFun<double> function(avars, ay);
    std::vector<double> point(static_cast<size_t>(var_count), 0.0);
    point[0] = density.rho;
    for (int i = 0; i < ncomp; ++i) {
        point[static_cast<size_t>(1 + i)] = x[static_cast<size_t>(i)];
    }
    auto values = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);
    auto hessian = function.Hessian(point, 0);
    std::vector<double> third_tensor =
        ares_detail::scalar_function_third_derivative_tensor_cppad(function, point);
    const auto h = [&](int row, int col) {
        return hessian[static_cast<size_t>(row * var_count + col)];
    };
    const auto h3 = [&](int leading, int row, int col) {
        return third_tensor[
            static_cast<size_t>(leading * var_count * var_count + row * var_count + col)
        ];
    };
    ares_detail::AssociationPhaseStateResponse association_response;
    if (active_association) {
        association_response = association_phase_state_response_cppad_cpp(t, density.rho, x, cppargs);
        if (!association_response.active) {
            out.message = "Association phase-state sensitivity expected active association but received no active site fractions.";
            out.density_backend = "selected_density_root";
            return out;
        }
    }

    (void)values;
    const double da_drho = jacobian[0];
    const double association_zraw = active_association ? association_response.zraw : 0.0;
    const double z_raw = density.rho * da_drho + association_zraw;
    const double z = 1.0 + z_raw;
    if (!(z > 0.0) || !std::isfinite(z)) {
        out.message = "Phase-state sensitivity produced a non-positive compressibility factor.";
        out.density_backend = "selected_density_root";
        return out;
    }

    std::vector<double> dadx(static_cast<size_t>(ncomp), 0.0);
    double sum_x_dadx = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        dadx[static_cast<size_t>(i)] = jacobian[static_cast<size_t>(1 + i)];
        sum_x_dadx += x[static_cast<size_t>(i)] * dadx[static_cast<size_t>(i)];
    }

    const auto base_dzraw_dvar = [&](int var_index) {
        double value = density.rho * h(0, var_index);
        if (var_index == 0) {
            value += da_drho;
        }
        return value;
    };
    const auto base_d2zraw_dvar2 = [&](int first, int second) {
        double value = density.rho * h3(0, first, second);
        if (first == 0) {
            value += h(0, second);
        }
        if (second == 0) {
            value += h(0, first);
        }
        return value;
    };
    const auto dzraw_dvar = [&](int var_index) {
        double value = base_dzraw_dvar(var_index);
        if (active_association) {
            value += association_response.dzraw_dvar[static_cast<size_t>(var_index)];
        }
        return value;
    };
    const auto d2zraw_dvar2 = [&](int first, int second) {
        double value = base_d2zraw_dvar2(first, second);
        if (active_association) {
            value += association_response.d2zraw_dvar2_row_major[
                static_cast<size_t>(first * var_count + second)
            ];
        }
        return value;
    };
    const auto pressure_dvar = [&](int var_index) {
        const double dzraw = dzraw_dvar(var_index);
        double value = density.rho * dzraw;
        if (var_index == 0) {
            value += z;
        }
        return kb * t * N_AV * value;
    };
    const auto pressure_d2var = [&](int first, int second) {
        double value = density.rho * d2zraw_dvar2(first, second);
        if (first == 0) {
            value += dzraw_dvar(second);
        }
        if (second == 0) {
            value += dzraw_dvar(first);
        }
        return kb * t * N_AV * value;
    };
    const auto sum_x_dadx_dvar = [&](int var_index) {
        double value = 0.0;
        if (var_index > 0) {
            value += dadx[static_cast<size_t>(var_index - 1)];
        }
        for (int k = 0; k < ncomp; ++k) {
            value += x[static_cast<size_t>(k)] * h(1 + k, var_index);
        }
        return value;
    };
    const auto d2sum_x_dadx_dvar2 = [&](int first, int second) {
        double value = 0.0;
        if (first > 0) {
            value += h(first, second);
        }
        if (second > 0) {
            value += h(second, first);
        }
        for (int k = 0; k < ncomp; ++k) {
            value += x[static_cast<size_t>(k)] * h3(1 + k, first, second);
        }
        return value;
    };
    const auto mu_direct_dvar = [&](int component, int var_index) {
        const double da_dvar = jacobian[static_cast<size_t>(var_index)];
        double dmu =
            da_dvar
            + base_dzraw_dvar(var_index)
            + h(1 + component, var_index)
            - sum_x_dadx_dvar(var_index);
        if (active_association) {
            dmu += association_response.dmu_dvar_row_major[
                static_cast<size_t>(component * var_count + var_index)
            ];
        }
        return dmu;
    };
    const auto mu_direct_d2var = [&](int component, int first, int second) {
        double value =
            h(first, second)
            + base_d2zraw_dvar2(first, second)
            + h3(1 + component, first, second)
            - d2sum_x_dadx_dvar2(first, second);
        if (active_association) {
            value += association_response.d2mu_dvar2_tensor_row_major[
                static_cast<size_t>(component * var_count * var_count + first * var_count + second)
            ];
        }
        return value;
    };
    const auto fixed_density_dlnphi_dvar = [&](int component, int var_index) {
        const double dzraw = dzraw_dvar(var_index);
        return mu_direct_dvar(component, var_index) - dzraw / z;
    };
    const auto fixed_density_d2lnphi_dvar2 = [&](int component, int first, int second) {
        const double dz_first = dzraw_dvar(first);
        const double dz_second = dzraw_dvar(second);
        const double d2z = d2zraw_dvar2(first, second);
        return mu_direct_d2var(component, first, second) - d2z / z + dz_first * dz_second / (z * z);
    };

    const double dp_drho = pressure_dvar(0);
    if (!std::isfinite(dp_drho) || std::abs(dp_drho) <= 1.0e-30) {
        out.message = "Phase-state sensitivity produced a singular pressure-density derivative.";
        out.density_backend = "selected_density_root";
        return out;
    }

    out.pressure_density_derivative = dp_drho;
    out.pressure_density_second_derivative = pressure_d2var(0, 0);
    out.density_composition_derivative.assign(static_cast<size_t>(ncomp), 0.0);
    out.density_composition_hessian_row_major.assign(static_cast<size_t>(ncomp * ncomp), 0.0);
    out.pressure_composition_fixed_density_derivative.assign(static_cast<size_t>(ncomp), 0.0);
    out.pressure_density_composition_cross_derivative.assign(static_cast<size_t>(ncomp), 0.0);
    out.pressure_composition_fixed_density_hessian_row_major.assign(
        static_cast<size_t>(ncomp * ncomp),
        0.0
    );
    out.ln_fugacity_density_derivative.assign(static_cast<size_t>(ncomp), 0.0);
    out.fixed_density_jacobian_row_major.assign(static_cast<size_t>(ncomp * ncomp), 0.0);
    out.fixed_density_hessian_tensor_row_major.assign(static_cast<size_t>(ncomp * ncomp * ncomp), 0.0);
    out.jacobian_row_major.assign(static_cast<size_t>(ncomp * ncomp), 0.0);
    out.hessian_tensor_row_major.assign(static_cast<size_t>(ncomp * ncomp * ncomp), 0.0);

    for (int j = 0; j < ncomp; ++j) {
        const int var_index = 1 + j;
        const double dp_dx_fixed = pressure_dvar(var_index);
        out.pressure_composition_fixed_density_derivative[static_cast<size_t>(j)] = dp_dx_fixed;
        out.pressure_density_composition_cross_derivative[static_cast<size_t>(j)] = pressure_d2var(0, var_index);
        out.density_composition_derivative[static_cast<size_t>(j)] = -dp_dx_fixed / dp_drho;
    }
    for (int j = 0; j < ncomp; ++j) {
        const int first = 1 + j;
        for (int k = 0; k < ncomp; ++k) {
            const int second = 1 + k;
            const double density_hessian =
                -(
                    pressure_d2var(first, second)
                    + pressure_d2var(0, first) * out.density_composition_derivative[static_cast<size_t>(k)]
                    + pressure_d2var(0, second) * out.density_composition_derivative[static_cast<size_t>(j)]
                    + pressure_d2var(0, 0)
                        * out.density_composition_derivative[static_cast<size_t>(j)]
                        * out.density_composition_derivative[static_cast<size_t>(k)]
                ) / dp_drho;
            out.density_composition_hessian_row_major[static_cast<size_t>(j * ncomp + k)] = density_hessian;
            out.pressure_composition_fixed_density_hessian_row_major[static_cast<size_t>(j * ncomp + k)] =
                pressure_d2var(first, second);
        }
    }
    for (int i = 0; i < ncomp; ++i) {
        const double dln_drho = fixed_density_dlnphi_dvar(i, 0);
        out.ln_fugacity_density_derivative[static_cast<size_t>(i)] = dln_drho;
        for (int j = 0; j < ncomp; ++j) {
            const double fixed = fixed_density_dlnphi_dvar(i, 1 + j);
            const double total = fixed + dln_drho * out.density_composition_derivative[static_cast<size_t>(j)];
            const size_t index = static_cast<size_t>(i * ncomp + j);
            out.fixed_density_jacobian_row_major[index] = fixed;
            out.jacobian_row_major[index] = total;
            for (int k = 0; k < ncomp; ++k) {
                const double fixed_hessian = fixed_density_d2lnphi_dvar2(i, 1 + j, 1 + k);
                out.fixed_density_hessian_tensor_row_major[
                    static_cast<size_t>(i * ncomp * ncomp + j * ncomp + k)
                ] = fixed_hessian;
                out.hessian_tensor_row_major[
                    static_cast<size_t>(i * ncomp * ncomp + j * ncomp + k)
                ] =
                    fixed_hessian
                    + fixed_density_d2lnphi_dvar2(i, 1 + j, 0)
                        * out.density_composition_derivative[static_cast<size_t>(k)]
                    + fixed_density_d2lnphi_dvar2(i, 1 + k, 0)
                        * out.density_composition_derivative[static_cast<size_t>(j)]
                    + fixed_density_d2lnphi_dvar2(i, 0, 0)
                        * out.density_composition_derivative[static_cast<size_t>(j)]
                        * out.density_composition_derivative[static_cast<size_t>(k)]
                    + dln_drho * out.density_composition_hessian_row_major[static_cast<size_t>(j * ncomp + k)];
            }
        }
    }

    out.supported = true;
    out.backend = "cppad_implicit";
    out.density_backend = "implicit_density_root";
    out.message = (cppargs.born_model == 2)
        ? "CppAD phase-state fugacity composition sensitivities with implicit density-root routing and SSM/DS Born terms are available."
        : "CppAD phase-state fugacity composition sensitivities with implicit density-root routing are available.";
    return out;
}

EosPhasePressureDerivativeResult eos_phase_pressure_derivatives_cpp(
    double t,
    const vector<double> &amounts,
    double volume,
    const add_args &cppargs
) {
    EosPhasePressureDerivativeResult out;
    using CppADScalar = CppAD::AD<double>;
    if (!(std::isfinite(t) && t > 0.0) || !(std::isfinite(volume) && volume > 0.0)) {
        out.message = "EOS phase pressure derivatives require positive finite temperature and volume.";
        return out;
    }
    const int ncomp = static_cast<int>(amounts.size());
    if (ncomp <= 0) {
        out.message = "EOS phase pressure derivatives require at least one species amount.";
        return out;
    }
    double total_amount = 0.0;
    for (double amount : amounts) {
        if (!(std::isfinite(amount) && amount > 0.0)) {
            out.message = "EOS phase pressure derivatives require positive finite species amounts.";
            return out;
        }
        total_amount += amount;
    }
    if (!(std::isfinite(total_amount) && total_amount > 0.0)) {
        out.message = "EOS phase pressure derivatives require a positive finite total amount.";
        return out;
    }

    const double rho = total_amount / volume;
    vector<double> x(static_cast<size_t>(ncomp), 0.0);
    for (int i = 0; i < ncomp; ++i) {
        x[static_cast<size_t>(i)] = amounts[static_cast<size_t>(i)] / total_amount;
    }

    const int var_count = ncomp + 1;
    std::vector<CppADScalar> avars(static_cast<size_t>(var_count));
    avars[0] = rho;
    for (int i = 0; i < ncomp; ++i) {
        avars[static_cast<size_t>(1 + i)] = x[static_cast<size_t>(i)];
    }
    CppAD::Independent(avars);

    std::vector<CppADScalar> ax(static_cast<size_t>(ncomp));
    for (int i = 0; i < ncomp; ++i) {
        ax[static_cast<size_t>(i)] = avars[static_cast<size_t>(1 + i)];
    }
    const bool active_association = has_association_sites(cppargs);
    add_args recording_args = cppargs;
    if (active_association) {
        recording_args.assoc_num.clear();
        recording_args.assoc_matrix.clear();
        recording_args.e_assoc.clear();
        recording_args.vol_a.clear();
        recording_args.k_hb.clear();
    }
    auto contributions = ares_detail::ares_contributions_scalar_cpp(t, avars[0], ax, recording_args);
    std::vector<CppADScalar> ay(1);
    ay[0] = contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born;

    CppAD::ADFun<double> function(avars, ay);
    std::vector<double> point(static_cast<size_t>(var_count), 0.0);
    point[0] = rho;
    for (int i = 0; i < ncomp; ++i) {
        point[static_cast<size_t>(1 + i)] = x[static_cast<size_t>(i)];
    }
    auto jacobian = function.Jacobian(point);
    auto hessian = function.Hessian(point, 0);
    const auto h = [&](int row, int col) {
        return hessian[static_cast<size_t>(row * var_count + col)];
    };

    ares_detail::AssociationPhaseStateResponse association_response;
    if (active_association) {
        association_response = association_phase_state_response_cppad_cpp(t, rho, x, cppargs);
        if (!association_response.active) {
            out.message =
                "Association pressure derivative expected active association but received no active site fractions.";
            return out;
        }
    }

    const double da_drho = jacobian[0];
    const double association_zraw = active_association ? association_response.zraw : 0.0;
    const double z_raw = rho * da_drho + association_zraw;
    const double z = 1.0 + z_raw;
    if (!(std::isfinite(z) && z > 0.0)) {
        out.message = "EOS phase pressure derivatives produced a non-positive compressibility factor.";
        return out;
    }

    const auto base_dzraw_dvar = [&](int var_index) {
        double value = rho * h(0, var_index);
        if (var_index == 0) {
            value += da_drho;
        }
        return value;
    };
    const auto dzraw_dvar = [&](int var_index) {
        double value = base_dzraw_dvar(var_index);
        if (active_association) {
            value += association_response.dzraw_dvar[static_cast<size_t>(var_index)];
        }
        return value;
    };
    const auto pressure_dvar = [&](int var_index) {
        const double dzraw = dzraw_dvar(var_index);
        double value = rho * dzraw;
        if (var_index == 0) {
            value += z;
        }
        return kb * t * N_AV * value;
    };

    const double dp_drho = pressure_dvar(0);
    if (!(std::isfinite(dp_drho) && std::abs(dp_drho) > 1.0e-30)) {
        out.message = "EOS phase pressure derivatives produced a singular pressure-density derivative.";
        return out;
    }

    vector<double> dpdx_fixed(static_cast<size_t>(ncomp), 0.0);
    double sum_x_dpdx = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        dpdx_fixed[static_cast<size_t>(i)] = pressure_dvar(1 + i);
        if (!std::isfinite(dpdx_fixed[static_cast<size_t>(i)])) {
            out.message = "EOS phase pressure composition derivative was not finite.";
            return out;
        }
        sum_x_dpdx += x[static_cast<size_t>(i)] * dpdx_fixed[static_cast<size_t>(i)];
    }

    out.pressure_jacobian_row_major.assign(static_cast<size_t>(ncomp + 1), 0.0);
    for (int species = 0; species < ncomp; ++species) {
        out.pressure_jacobian_row_major[static_cast<size_t>(species)] =
            dp_drho / volume
            + (dpdx_fixed[static_cast<size_t>(species)] - sum_x_dpdx) / total_amount;
    }
    out.pressure_jacobian_row_major[static_cast<size_t>(ncomp)] = -dp_drho * rho / volume;
    for (double value : out.pressure_jacobian_row_major) {
        if (!std::isfinite(value)) {
            out.message = "EOS phase pressure amount-volume derivative was not finite.";
            return out;
        }
    }
    out.supported = true;
    out.backend = active_association ? "cppad_implicit" : "cppad";
    out.message = active_association
        ? "CppAD pressure derivatives with implicit association sensitivities are available."
        : "CppAD pressure derivatives are available.";
    out.pressure_density_derivative = dp_drho;
    return out;
}

NeutralBinaryKijPhaseDerivatives neutral_binary_kij_phase_derivatives_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs,
    int k_index
) {
    return neutral_binary_pair_parameter_phase_derivatives_cpp(t, rho, x, cppargs, k_index, "k_ij");
}

NeutralBinaryKijPhaseDerivatives neutral_binary_pair_parameter_phase_derivatives_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs,
    int parameter_index,
    const std::string &parameter_name
) {
#ifdef EPCSAFT_HAS_CPPAD
    using CppADScalar = CppAD::AD<double>;
    const int ncomp = static_cast<int>(x.size());
    if (ncomp != 2 || cppargs.m.size() != 2 || cppargs.s.size() != 2 || cppargs.e.size() != 2) {
        throw ValueError("unsupported: native binary pair-parameter CppAD derivatives require exactly two neutral components.");
    }
    if (!cppargs.z.empty()) {
        for (double charge : cppargs.z) {
            if (std::abs(charge) > 1.0e-12) {
                throw ValueError("unsupported: native binary pair-parameter CppAD derivatives do not support ionic components.");
            }
        }
    }
    const bool is_kij = parameter_name == "k_ij";
    const bool is_lij = parameter_name == "l_ij";
    if (!is_kij && !is_lij) {
        throw ValueError("Native binary pair-parameter derivative supports only k_ij and l_ij.");
    }
    const bool active_association = has_association_sites(cppargs);
    if (active_association && !is_kij) {
        throw ValueError("unsupported: associating binary pair-parameter derivatives support k_ij only until l_ij association-size sensitivities are implemented.");
    }
    const vector<double> &parameter_matrix = is_kij ? cppargs.k_ij : cppargs.l_ij;
    if (parameter_matrix.size() != static_cast<size_t>(ncomp * ncomp)) {
        throw ValueError("unsupported: native binary pair-parameter CppAD derivatives require a dense parameter matrix.");
    }
    if (parameter_index < 0 || static_cast<size_t>(parameter_index) >= parameter_matrix.size()) {
        throw ValueError("Native binary pair-parameter derivative index is out of range.");
    }
    if (!(t > 0.0) || !(rho > 0.0) || x.size() != 2 || !(x[0] > 0.0) || !(x[1] > 0.0)) {
        throw ValueError("Native binary k_ij derivative evaluation requires positive T, rho, and composition values.");
    }

    constexpr int kRhoIndex = 0;
    constexpr int kKijIndex = 1;
    constexpr int kX0Index = 2;
    constexpr int kX1Index = 3;
    constexpr int kVarCount = 4;
    std::vector<CppADScalar> avars(kVarCount);
    avars[kRhoIndex] = rho;
    avars[kKijIndex] = parameter_matrix[static_cast<size_t>(parameter_index)];
    avars[kX0Index] = x[0];
    avars[kX1Index] = x[1];
    CppAD::Independent(avars);

    add_args recording_args = cppargs;
    if (active_association && is_kij) {
        recording_args.assoc_num.clear();
        recording_args.assoc_matrix.clear();
        recording_args.e_assoc.clear();
        recording_args.vol_a.clear();
        recording_args.k_hb.clear();
    }

    std::vector<CppADScalar> ax = {avars[kX0Index], avars[kX1Index]};
    auto contributions = ares_detail::ares_contributions_scalar_cpp(
        t,
        avars[kRhoIndex],
        ax,
        recording_args,
        is_kij ? parameter_index : -1,
        is_kij ? &avars[kKijIndex] : nullptr,
        is_lij ? parameter_index : -1,
        is_lij ? &avars[kKijIndex] : nullptr
    );
    std::vector<CppADScalar> ay(1);
    ay[0] = contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born;

    CppAD::ADFun<double> function(avars, ay);
    std::vector<double> point = {rho, parameter_matrix[static_cast<size_t>(parameter_index)], x[0], x[1]};
    auto values = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);
    auto hessian = function.Hessian(point, 0);

    const double ares = values[0];
    const double da_drho = jacobian[kRhoIndex];
    const double da_dk = jacobian[kKijIndex];
    const double da_dx0 = jacobian[kX0Index];
    const double da_dx1 = jacobian[kX1Index];
    const auto h = [&](int row, int col) {
        return hessian[static_cast<size_t>(row * kVarCount + col)];
    };
    const double d2a_drho2 = h(kRhoIndex, kRhoIndex);
    const double d2a_drho_dk = h(kRhoIndex, kKijIndex);
    const double d2a_dx0_drho = h(kX0Index, kRhoIndex);
    const double d2a_dx1_drho = h(kX1Index, kRhoIndex);
    const double d2a_dx0_dk = h(kX0Index, kKijIndex);
    const double d2a_dx1_dk = h(kX1Index, kKijIndex);

    const double z_raw = rho * da_drho;
    const double z = 1.0 + z_raw;
    if (!(z > 0.0)) {
        throw ValueError("Native binary k_ij derivative evaluation produced non-positive Z.");
    }
    const double dz_drho = da_drho + rho * d2a_drho2;
    const double dz_dk = rho * d2a_drho_dk;
    const double pressure_factor = kb * t * N_AV;
    double z_for_values = z;
    vector<double> mu_for_values;
    vector<double> lnfug_for_values;
    double dpdrho_for_values = pressure_factor * (z + rho * dz_drho);
    double assoc_dzraw_drho_for_values = 0.0;
    vector<double> assoc_dmu_drho_for_values;
    if (active_association && is_kij) {
        ResidualChemicalPotentialResult native_mu =
            residual_chemical_potential_result_cpp(t, rho, vector<double>(x.begin(), x.end()), cppargs);
        FugacityContributionResult native_fugacity =
            fugacity_coefficient_result_cpp(t, rho, vector<double>(x.begin(), x.end()), cppargs);
        ares_detail::AssociationDensityResponse assoc_density =
            association_density_response_cppad_cpp(t, rho, vector<double>(x.begin(), x.end()), cppargs);
        z_for_values = native_mu.composition.z.total;
        mu_for_values = native_mu.mu.total;
        lnfug_for_values = native_fugacity.lnfugcoef.total;
        assoc_dzraw_drho_for_values = assoc_density.dzraw_drho;
        assoc_dmu_drho_for_values = std::move(assoc_density.dmu_drho);
        dpdrho_for_values = pressure_factor * (z_for_values + rho * (dz_drho + assoc_dzraw_drho_for_values));
    }
    NeutralBinaryKijPhaseDerivatives out;
    out.rho = rho;
    out.z = z_for_values;
    out.pressure = rho * pressure_factor * z_for_values;
    out.dpdrho = dpdrho_for_values;
    out.dpdk = rho * pressure_factor * dz_dk;
    if (!(std::isfinite(out.dpdrho)) || std::abs(out.dpdrho) <= 0.0) {
        throw ValueError("Native binary k_ij derivative evaluation produced invalid dP/drho.");
    }
    out.drhodk = -out.dpdk / out.dpdrho;
    out.mu_res.assign(2, 0.0);
    out.dmu_res_dk_fixed_rho.assign(2, 0.0);
    out.lnphi.assign(2, 0.0);
    out.dlnphi_drho.assign(2, 0.0);
    out.dlnphi_dk_fixed_rho.assign(2, 0.0);
    out.dlnphi_dk_total.assign(2, 0.0);

    const std::array<double, 2> dadx = {da_dx0, da_dx1};
    const std::array<double, 2> dadx_drho = {d2a_dx0_drho, d2a_dx1_drho};
    const std::array<double, 2> dadx_dk = {d2a_dx0_dk, d2a_dx1_dk};
    const double sum_x_dadx = x[0] * dadx[0] + x[1] * dadx[1];
    const double sum_x_dadx_drho = x[0] * dadx_drho[0] + x[1] * dadx_drho[1];
    const double sum_x_dadx_dk = x[0] * dadx_dk[0] + x[1] * dadx_dk[1];
    for (int i = 0; i < 2; ++i) {
        const double mu = ares + z_raw + dadx[static_cast<size_t>(i)] - sum_x_dadx;
        const double dmu_drho = da_drho + dz_drho + dadx_drho[static_cast<size_t>(i)] - sum_x_dadx_drho;
        const double dmu_dk = da_dk + dz_dk + dadx_dk[static_cast<size_t>(i)] - sum_x_dadx_dk;
        out.mu_res[static_cast<size_t>(i)] =
            mu_for_values.empty() ? mu : mu_for_values[static_cast<size_t>(i)];
        out.dmu_res_dk_fixed_rho[static_cast<size_t>(i)] = dmu_dk;
        out.lnphi[static_cast<size_t>(i)] =
            lnfug_for_values.empty() ? mu - std::log(z) : lnfug_for_values[static_cast<size_t>(i)];
        out.dlnphi_drho[static_cast<size_t>(i)] = assoc_dmu_drho_for_values.empty()
            ? dmu_drho - dz_drho / z_for_values
            : dmu_drho + assoc_dmu_drho_for_values[static_cast<size_t>(i)]
                - (dz_drho + assoc_dzraw_drho_for_values) / z_for_values;
        out.dlnphi_dk_fixed_rho[static_cast<size_t>(i)] = dmu_dk - dz_dk / z_for_values;
        out.dlnphi_dk_total[static_cast<size_t>(i)] =
            out.dlnphi_dk_fixed_rho[static_cast<size_t>(i)] + out.dlnphi_drho[static_cast<size_t>(i)] * out.drhodk;
    }
    out.backend = active_association ? "cppad_implicit" : "cppad";
    return out;
#else
    (void)t;
    (void)rho;
    (void)x;
    (void)cppargs;
    (void)parameter_index;
    (void)parameter_name;
    throw ValueError("unsupported: CppAD support is not enabled in this native build.");
#endif
}

NeutralBinaryKijPhaseDerivatives generic_component_parameter_phase_derivatives_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs,
    int target_kind,
    int target_index
) {
#ifdef EPCSAFT_HAS_CPPAD
    using CppADScalar = CppAD::AD<double>;
    const int ncomp = static_cast<int>(x.size());
    if (ncomp <= 0 || cppargs.s.size() != x.size() || cppargs.e.size() != x.size()) {
        throw ValueError("unsupported: generic component-parameter CppAD derivatives require aligned component parameters.");
    }
    if (target_kind != ares_detail::kGenericTargetSLocal
        && target_kind != ares_detail::kGenericTargetELocal
        && target_kind != ares_detail::kGenericTargetDBornLocal
        && target_kind != ares_detail::kGenericTargetFSolvLocal
        && target_kind != ares_detail::kGenericTargetDielcLocal) {
        throw ValueError("unsupported: generic component-parameter CppAD derivatives support s, e, d_born, f_solv, and dielc only.");
    }
    if (target_index < 0 || target_index >= ncomp) {
        throw ValueError("Native generic component-parameter derivative target index is out of range.");
    }
    if (target_kind == ares_detail::kGenericTargetDBornLocal) {
        if (cppargs.d_born.size() != x.size()) {
            throw ValueError("unsupported: d_born CppAD derivatives require aligned d_born values.");
        }
        if (cppargs.born_model != 1 && cppargs.born_model != 2) {
            throw ValueError("unsupported: d_born CppAD derivatives require a direct or SSM/DS Born model.");
        }
    }
    if (target_kind == ares_detail::kGenericTargetFSolvLocal) {
        if (cppargs.f_solv.size() != x.size()) {
            throw ValueError("unsupported: f_solv CppAD derivatives require aligned f_solv values.");
        }
        if (cppargs.born_model != 2) {
            throw ValueError("unsupported: f_solv CppAD derivatives require the SSM/DS Born model.");
        }
    }
    if (target_kind == ares_detail::kGenericTargetDielcLocal) {
        if (cppargs.dielc.size() != x.size()) {
            throw ValueError("unsupported: dielc CppAD derivatives require aligned relative-permittivity values.");
        }
        if (cppargs.dielc_rule != 1) {
            throw ValueError("unsupported: dielc CppAD derivatives require linear mole-fraction relative-permittivity mixing.");
        }
    }
    if (!(t > 0.0) || !(rho > 0.0)) {
        throw ValueError("Native generic component-parameter derivative evaluation requires positive T and rho.");
    }
    for (double xi : x) {
        if (!(xi >= 0.0)) {
            throw ValueError("Native generic component-parameter derivative evaluation requires nonnegative composition.");
        }
    }

    const int rho_index = 0;
    const int theta_index = 1;
    const int x_start = 2;
    const int var_count = x_start + ncomp;
    std::vector<CppADScalar> avars(static_cast<size_t>(var_count));
    avars[static_cast<size_t>(rho_index)] = rho;
    double theta0 = 0.0;
    if (target_kind == ares_detail::kGenericTargetSLocal) {
        theta0 = cppargs.s[static_cast<size_t>(target_index)];
    } else if (target_kind == ares_detail::kGenericTargetELocal) {
        theta0 = cppargs.e[static_cast<size_t>(target_index)];
    } else if (target_kind == ares_detail::kGenericTargetDBornLocal) {
        theta0 = cppargs.d_born[static_cast<size_t>(target_index)];
    } else if (target_kind == ares_detail::kGenericTargetFSolvLocal) {
        theta0 = cppargs.f_solv[static_cast<size_t>(target_index)];
    } else {
        theta0 = cppargs.dielc[static_cast<size_t>(target_index)];
    }
    avars[static_cast<size_t>(theta_index)] = theta0;
    for (int i = 0; i < ncomp; ++i) {
        avars[static_cast<size_t>(x_start + i)] = x[static_cast<size_t>(i)];
    }
    CppAD::Independent(avars);

    std::vector<CppADScalar> ax(static_cast<size_t>(ncomp));
    for (int i = 0; i < ncomp; ++i) {
        ax[static_cast<size_t>(i)] = avars[static_cast<size_t>(x_start + i)];
    }
    const CppADScalar *no_pair_override = nullptr;
    auto contributions = ares_detail::ares_contributions_scalar_cpp(
        t,
        avars[static_cast<size_t>(rho_index)],
        ax,
        cppargs,
        -1,
        no_pair_override,
        -1,
        no_pair_override,
        target_kind,
        target_index,
        &avars[static_cast<size_t>(theta_index)]
    );
    std::vector<CppADScalar> ay(1);
    ay[0] = contributions.hc + contributions.disp + contributions.assoc + contributions.ion + contributions.born;

    CppAD::ADFun<double> function(avars, ay);
    std::vector<double> point(static_cast<size_t>(var_count), 0.0);
    point[static_cast<size_t>(rho_index)] = rho;
    point[static_cast<size_t>(theta_index)] = theta0;
    for (int i = 0; i < ncomp; ++i) {
        point[static_cast<size_t>(x_start + i)] = x[static_cast<size_t>(i)];
    }
    auto values = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);
    auto hessian = function.Hessian(point, 0);
    const auto h = [&](int row, int col) {
        return hessian[static_cast<size_t>(row * var_count + col)];
    };

    const double ares = values[0];
    const double da_drho = jacobian[static_cast<size_t>(rho_index)];
    const double da_dtheta = jacobian[static_cast<size_t>(theta_index)];
    const double d2a_drho2 = h(rho_index, rho_index);
    const double d2a_drho_dtheta = h(rho_index, theta_index);
    const double z_raw = rho * da_drho;
    const double z = 1.0 + z_raw;
    if (!(z > 0.0)) {
        throw ValueError("Native generic component-parameter derivative evaluation produced non-positive Z.");
    }
    const double dz_drho = da_drho + rho * d2a_drho2;
    const double dz_dtheta = rho * d2a_drho_dtheta;
    const double pressure_factor = kb * t * N_AV;
    NeutralBinaryKijPhaseDerivatives out;
    out.rho = rho;
    out.z = z;
    out.pressure = rho * pressure_factor * z;
    out.dpdrho = pressure_factor * (z + rho * dz_drho);
    out.dpdk = rho * pressure_factor * dz_dtheta;
    if (!(std::isfinite(out.dpdrho)) || std::abs(out.dpdrho) <= 0.0) {
        throw ValueError("Native generic component-parameter derivative evaluation produced invalid dP/drho.");
    }
    out.drhodk = -out.dpdk / out.dpdrho;
    out.mu_res.assign(static_cast<size_t>(ncomp), 0.0);
    out.dmu_res_dk_fixed_rho.assign(static_cast<size_t>(ncomp), 0.0);
    out.lnphi.assign(static_cast<size_t>(ncomp), 0.0);
    out.dlnphi_drho.assign(static_cast<size_t>(ncomp), 0.0);
    out.dlnphi_dk_fixed_rho.assign(static_cast<size_t>(ncomp), 0.0);
    out.dlnphi_dk_total.assign(static_cast<size_t>(ncomp), 0.0);

    vector<double> dadx(static_cast<size_t>(ncomp), 0.0);
    vector<double> dadx_drho(static_cast<size_t>(ncomp), 0.0);
    vector<double> dadx_dtheta(static_cast<size_t>(ncomp), 0.0);
    double sum_x_dadx = 0.0;
    double sum_x_dadx_drho = 0.0;
    double sum_x_dadx_dtheta = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        const int xi = x_start + i;
        dadx[static_cast<size_t>(i)] = jacobian[static_cast<size_t>(xi)];
        dadx_drho[static_cast<size_t>(i)] = h(xi, rho_index);
        dadx_dtheta[static_cast<size_t>(i)] = h(xi, theta_index);
        sum_x_dadx += x[static_cast<size_t>(i)] * dadx[static_cast<size_t>(i)];
        sum_x_dadx_drho += x[static_cast<size_t>(i)] * dadx_drho[static_cast<size_t>(i)];
        sum_x_dadx_dtheta += x[static_cast<size_t>(i)] * dadx_dtheta[static_cast<size_t>(i)];
    }
    for (int i = 0; i < ncomp; ++i) {
        const double mu = ares + z_raw + dadx[static_cast<size_t>(i)] - sum_x_dadx;
        const double dmu_drho = da_drho + dz_drho + dadx_drho[static_cast<size_t>(i)] - sum_x_dadx_drho;
        const double dmu_dtheta = da_dtheta + dz_dtheta + dadx_dtheta[static_cast<size_t>(i)] - sum_x_dadx_dtheta;
        out.mu_res[static_cast<size_t>(i)] = mu;
        out.dmu_res_dk_fixed_rho[static_cast<size_t>(i)] = dmu_dtheta;
        out.lnphi[static_cast<size_t>(i)] = mu - std::log(z);
        out.dlnphi_drho[static_cast<size_t>(i)] = dmu_drho - dz_drho / z;
        out.dlnphi_dk_fixed_rho[static_cast<size_t>(i)] = dmu_dtheta - dz_dtheta / z;
        out.dlnphi_dk_total[static_cast<size_t>(i)] =
            out.dlnphi_dk_fixed_rho[static_cast<size_t>(i)] + out.dlnphi_drho[static_cast<size_t>(i)] * out.drhodk;
    }
    out.backend = "cppad_implicit";
    return out;
#else
    (void)t;
    (void)rho;
    (void)x;
    (void)cppargs;
    (void)target_kind;
    (void)target_index;
    throw ValueError("unsupported: CppAD support is not enabled in this native build.");
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
    BornIntermediateState born_state = born_intermediate_state_cpp(t, x, cppargs, true);

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
    BornIntermediateState born_state = born_intermediate_state_cpp(t, x, cppargs, false);

    ContributionDadxResult hc = dadx_hc_cpp(thermo, hc_state, t, rho, x, cppargs);
    ContributionDadxResult disp = dadx_disp_cpp(thermo, hc_state, dispersion, t, rho, x, cppargs);
    ContributionDadxResult assoc = dadx_assoc_cpp(thermo, hc_state, assoc_state, t, rho, x, cppargs);
    ContributionDadxResult ion = dadx_ion_cpp(ion_state, t, rho, x, cppargs);
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
        if (item.second == "unsupported") {
            result.derivative_available = false;
            break;
        }
    }
    return result;
}
