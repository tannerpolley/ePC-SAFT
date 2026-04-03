#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <array>
#include <numeric>
#include "math.h"
#include "Eigen/Dense"

#include "pcsaft_electrolyte.h"

using std::vector;

namespace {
thread_local vector<double> g_last_mu_hc;
thread_local vector<double> g_last_mu_disp;
thread_local vector<double> g_last_mu_polar;
thread_local vector<double> g_last_mu_assoc;
thread_local vector<double> g_last_mu_ion;
thread_local vector<double> g_last_mu_born;
thread_local vector<double> g_last_mu_total;
thread_local vector<double> g_last_lnfug_hc;
thread_local vector<double> g_last_lnfug_disp;
thread_local vector<double> g_last_lnfug_polar;
thread_local vector<double> g_last_lnfug_assoc;
thread_local vector<double> g_last_lnfug_ion;
thread_local vector<double> g_last_lnfug_born;
thread_local vector<double> g_last_lnfugcoef;
thread_local vector<double> g_last_dadx_hc;
thread_local vector<double> g_last_dadx_disp;
thread_local vector<double> g_last_dadx_polar;
thread_local vector<double> g_last_dadx_assoc;
thread_local vector<double> g_last_dadx_ion;
thread_local vector<double> g_last_dadx_born;
thread_local double g_last_a_hc = 0.0;
thread_local double g_last_a_disp = 0.0;
thread_local double g_last_a_polar = 0.0;
thread_local double g_last_a_assoc = 0.0;
thread_local double g_last_a_ion = 0.0;
thread_local double g_last_a_born = 0.0;
thread_local double g_last_sum_x_dadx_hc = 0.0;
thread_local double g_last_sum_x_dadx_disp = 0.0;
thread_local double g_last_sum_x_dadx_polar = 0.0;
thread_local double g_last_sum_x_dadx_assoc = 0.0;
thread_local double g_last_sum_x_dadx_ion = 0.0;
thread_local double g_last_sum_x_dadx_born = 0.0;
thread_local double g_last_z_raw_hc = 0.0;
thread_local double g_last_z_raw_disp = 0.0;
thread_local double g_last_z_raw_polar = 0.0;
thread_local double g_last_z_raw_assoc = 0.0;
thread_local double g_last_z_raw_ion = 0.0;
thread_local double g_last_z_raw_born = 0.0;
thread_local double g_last_z_hc = 0.0;
thread_local double g_last_z_disp = 0.0;
thread_local double g_last_z_polar = 0.0;
thread_local double g_last_z_assoc = 0.0;
thread_local double g_last_z_ion = 0.0;
thread_local double g_last_z_born = 0.0;

enum class AresContributionKind {
    HC,
    DISP,
    POLAR,
    ASSOC,
    ION,
    BORN
};

struct AresContributions {
    double hc = 0.0;
    double disp = 0.0;
    double polar = 0.0;
    double assoc = 0.0;
    double ion = 0.0;
    double born = 0.0;
};
thread_local double g_last_z_total = 1.0;

constexpr std::array<double, 7> kDispersionA0 = {
    0.9105631445, 0.6361281449, 2.6861347891, -26.547362491,
    97.759208784, -159.59154087, 91.297774084
};
constexpr std::array<double, 7> kDispersionA1 = {
    -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629,
    -65.255885330, 83.318680481, -33.746922930
};
constexpr std::array<double, 7> kDispersionA2 = {
    -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131,
    -4.1302112531, 13.776631870, -8.6728470368
};
constexpr std::array<double, 7> kDispersionB0 = {
    0.7240946941, 2.2382791861, -4.0025849485, -21.003576815,
    26.855641363, 206.55133841, -355.60235612
};
constexpr std::array<double, 7> kDispersionB1 = {
    -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648,
    192.67226447, -161.82646165, -165.20769346
};
constexpr std::array<double, 7> kDispersionB2 = {
    0.0976883116, -0.2557574982, -9.1558561530, 20.642075974,
    -38.804430052, 93.626774077, -29.666905585
};
constexpr std::array<double, 5> kDipoleA0 = { 0.3043504, -0.1358588, 1.4493329, 0.3556977, -2.0653308 };
constexpr std::array<double, 5> kDipoleA1 = { 0.9534641, -1.8396383, 2.0131180, -7.3724958, 8.2374135 };
constexpr std::array<double, 5> kDipoleA2 = { -1.1610080, 4.5258607, 0.9751222, -12.281038, 5.9397575 };
constexpr std::array<double, 5> kDipoleB0 = { 0.2187939, -1.1896431, 1.1626889, 0.0, 0.0 };
constexpr std::array<double, 5> kDipoleB1 = { -0.5873164, 1.2489132, -0.5085280, 0.0, 0.0 };
constexpr std::array<double, 5> kDipoleB2 = { 3.4869576, -14.915974, 15.372022, 0.0, 0.0 };
constexpr std::array<double, 5> kDipoleC0 = { -0.0646774, 0.1975882, -0.8087562, 0.6902849, 0.0 };
constexpr std::array<double, 5> kDipoleC1 = { -0.9520876, 2.9924258, -2.3802636, -0.2701261, 0.0 };
constexpr std::array<double, 5> kDipoleC2 = { -0.6260979, 1.2924686, 1.6542783, -3.4396744, 0.0 };
constexpr double kDipoleConversion = 7242.702976750923;

template <size_t N>
double evaluate_polynomial(const std::array<double, N> &coeffs, double x) {
    double value = 0.0;
    for (size_t i = 0; i < N; ++i) {
        value += coeffs[i] * std::pow(x, static_cast<int>(i));
    }
    return value;
}

template <size_t N>
double evaluate_polynomial_derivative(const std::array<double, N> &coeffs, double x) {
    double value = 0.0;
    for (size_t i = 1; i < N; ++i) {
        value += coeffs[i] * static_cast<double>(i) * std::pow(x, static_cast<int>(i - 1));
    }
    return value;
}

template <size_t N>
double evaluate_eta_weighted_derivative(const std::array<double, N> &coeffs, double x) {
    double value = 0.0;
    for (size_t i = 0; i < N; ++i) {
        value += coeffs[i] * static_cast<double>(i + 1) * std::pow(x, static_cast<int>(i));
    }
    return value;
}

struct ThermoCommonState {
    vector<double> d;
    vector<double> dd_dt;
    vector<double> zeta;
    vector<double> dzeta_dt;
    vector<double> e_ij;
    vector<double> s_ij;
    vector<double> ghs;
    vector<double> denghs;
    double den = 0.0;
    double eta = 0.0;
    double m_avg = 0.0;
    double m2es3 = 0.0;
    double m2e2s3 = 0.0;
};

struct DispersionPolynomialState {
    std::array<double, 7> a{};
    std::array<double, 7> b{};
    double I1 = 0.0;
    double I2 = 0.0;
    double dI1_deta = 0.0;
    double dI2_deta = 0.0;
    double dEtaI1_deta = 0.0;
    double dEtaI2_deta = 0.0;
    double C1 = 0.0;
    double C2 = 0.0;
};

struct AssociationSetup {
    vector<int> site_component_index;
    vector<double> x_assoc;
    vector<double> delta_ij;
};

double compute_pair_sigma(size_t idx, int i, int j, const add_args &cppargs) {
    double sigma = 0.5 * (cppargs.s[i] + cppargs.s[j]);
    if (!cppargs.l_ij.empty()) {
        sigma *= (1.0 - cppargs.l_ij[idx]);
    }
    return sigma;
}

double compute_pair_epsilon(size_t idx, int i, int j, const add_args &cppargs) {
    if (!cppargs.z.empty() && cppargs.z[i] * cppargs.z[j] > 0.0) {
        return 0.0;
    }
    double epsilon = std::sqrt(cppargs.e[i] * cppargs.e[j]);
    if (!cppargs.k_ij.empty()) {
        epsilon *= (1.0 - cppargs.k_ij[idx]);
    }
    return epsilon;
}

double compute_pair_diameter(double d_i, double d_j) {
    return d_i * d_j / (d_i + d_j);
}

double compute_hs_contact_value(double pair_diameter, double zeta2, double zeta3) {
    return 1.0 / (1.0 - zeta3)
        + pair_diameter * 3.0 * zeta2 / std::pow(1.0 - zeta3, 2.0)
        + std::pow(pair_diameter, 2.0) * 2.0 * zeta2 * zeta2 / std::pow(1.0 - zeta3, 3.0);
}

double compute_hs_contact_density_derivative(double pair_diameter, double zeta2, double zeta3) {
    return zeta3 / std::pow(1.0 - zeta3, 2.0)
        + pair_diameter * (3.0 * zeta2 / std::pow(1.0 - zeta3, 2.0) + 6.0 * zeta2 * zeta3 / std::pow(1.0 - zeta3, 3.0))
        + std::pow(pair_diameter, 2.0) * (4.0 * zeta2 * zeta2 / std::pow(1.0 - zeta3, 3.0) + 6.0 * zeta2 * zeta2 * zeta3 / std::pow(1.0 - zeta3, 4.0));
}

double compute_hs_contact_time_derivative(
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

ThermoCommonState build_thermo_common_state(double t, double rho, const vector<double> &x, const add_args &cppargs, bool include_dt) {
    ThermoCommonState state;
    int ncomp = static_cast<int>(x.size());
    state.d.assign(ncomp, 0.0);
    if (include_dt) {
        state.dd_dt.assign(ncomp, 0.0);
        state.dzeta_dt.assign(4, 0.0);
    }
    state.zeta.assign(4, 0.0);
    state.e_ij.assign(ncomp * ncomp, 0.0);
    state.s_ij.assign(ncomp * ncomp, 0.0);
    state.ghs.assign(ncomp * ncomp, 0.0);
    state.denghs.assign(ncomp * ncomp, 0.0);
    state.den = rho * N_AV / 1.0e30;

    for (int i = 0; i < ncomp; ++i) {
        state.d[i] = cppargs.s[i] * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
    }

    if (include_dt) {
        for (int i = 0; i < ncomp; ++i) {
            state.dd_dt[i] = -0.36 * cppargs.s[i] * cppargs.e[i] * std::exp(-3.0 * cppargs.e[i] / t) / (t * t);
        }
    }

    for (int k = 0; k < 4; ++k) {
        double summ = 0.0;
        for (int j = 0; j < ncomp; ++j) {
            summ += x[j] * cppargs.m[j] * std::pow(state.d[j], k);
        }
        state.zeta[k] = PI / 6.0 * state.den * summ;
    }

    if (include_dt) {
        for (int k = 1; k < 4; ++k) {
            double summ = 0.0;
            for (int j = 0; j < ncomp; ++j) {
                summ += x[j] * cppargs.m[j] * k * state.dd_dt[j] * std::pow(state.d[j], k - 1);
            }
            state.dzeta_dt[k] = PI / 6.0 * state.den * summ;
        }
    }

    state.eta = state.zeta[3];
    state.m_avg = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        state.m_avg += x[i] * cppargs.m[i];
    }

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        if (!cppargs.z.empty() && std::abs(cppargs.z[i]) > 1e-12) {
            int mode = cppargs.d_ion_mode;
            double sigma_i = cppargs.s[i];
            if (sigma_i <= 0.0) {
                throw ValueError("DH/ion diameter requires positive ionic sigma_i.");
            }
            if (mode == 0) {
                state.d[i] = sigma_i;
            }
            else if (mode == 1) {
                state.d[i] = sigma_i * (1.0 - 0.12);
            }
            else if (mode == 2) {
                state.d[i] = sigma_i * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
            }
            else {
                throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
            }
        }
        for (int j = 0; j < ncomp; ++j) {
            idx += 1;
            state.s_ij[idx] = compute_pair_sigma(static_cast<size_t>(idx), i, j, cppargs);
            state.e_ij[idx] = compute_pair_epsilon(static_cast<size_t>(idx), i, j, cppargs);
            state.m2es3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * state.e_ij[idx] / t * std::pow(state.s_ij[idx], 3);
            state.m2e2s3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * std::pow(state.e_ij[idx] / t, 2) * std::pow(state.s_ij[idx], 3);

            double pair_diameter = compute_pair_diameter(state.d[i], state.d[j]);
            state.ghs[idx] = compute_hs_contact_value(pair_diameter, state.zeta[2], state.zeta[3]);
            state.denghs[idx] = compute_hs_contact_density_derivative(pair_diameter, state.zeta[2], state.zeta[3]);
        }
    }

    if (include_dt) {
        for (int i = 0; i < ncomp; ++i) {
            if (!cppargs.z.empty() && std::abs(cppargs.z[i]) > 1e-12 && cppargs.d_ion_mode != 2) {
                state.dd_dt[i] = 0.0;
            }
        }
    }

    return state;
}

DispersionPolynomialState build_dispersion_polynomials(double m_avg, double eta) {
    DispersionPolynomialState state;
    double c1 = (m_avg - 1.0) / m_avg;
    double c2 = (m_avg - 2.0) / m_avg;
    for (size_t i = 0; i < state.a.size(); ++i) {
        state.a[i] = kDispersionA0[i] + c1 * kDispersionA1[i] + c1 * c2 * kDispersionA2[i];
        state.b[i] = kDispersionB0[i] + c1 * kDispersionB1[i] + c1 * c2 * kDispersionB2[i];
        state.I1 += state.a[i] * std::pow(eta, static_cast<int>(i));
        state.I2 += state.b[i] * std::pow(eta, static_cast<int>(i));
        if (i > 0) {
            state.dI1_deta += state.a[i] * static_cast<double>(i) * std::pow(eta, static_cast<int>(i - 1));
            state.dI2_deta += state.b[i] * static_cast<double>(i) * std::pow(eta, static_cast<int>(i - 1));
        }
        state.dEtaI1_deta += state.a[i] * static_cast<double>(i + 1) * std::pow(eta, static_cast<int>(i));
        state.dEtaI2_deta += state.b[i] * static_cast<double>(i + 1) * std::pow(eta, static_cast<int>(i));
    }
    state.C1 = 1.0 / (1.0
        + m_avg * (8.0 * eta - 2.0 * eta * eta) / std::pow(1.0 - eta, 4.0)
        + (1.0 - m_avg) * (20.0 * eta - 27.0 * eta * eta + 12.0 * std::pow(eta, 3.0) - 2.0 * std::pow(eta, 4.0))
            / std::pow((1.0 - eta) * (2.0 - eta), 2.0));
    state.C2 = -state.C1 * state.C1 * (
        m_avg * (-4.0 * eta * eta + 20.0 * eta + 8.0) / std::pow(1.0 - eta, 5.0)
        + (1.0 - m_avg) * (2.0 * std::pow(eta, 3.0) + 12.0 * eta * eta - 48.0 * eta + 40.0)
            / std::pow((1.0 - eta) * (2.0 - eta), 3.0));
    return state;
}

vector<double> build_dipole_coefficients(const std::array<double, 5> &c0, const std::array<double, 5> &c1, const std::array<double, 5> &c2, double m) {
    vector<double> coeffs(5, 0.0);
    double a = (m - 1.0) / m;
    double b = (m - 2.0) / m;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        coeffs[i] = c0[i] + a * c1[i] + a * b * c2[i];
    }
    return coeffs;
}

double compute_dh_kappa(double den, double t, double eps, double q2_sum) {
    return std::sqrt(den * E_CHRG * E_CHRG / kb / t / (eps * perm_vac) * q2_sum);
}

double compute_dh_chi(double kappa, double diameter) {
    double ka = kappa * diameter;
    return 3.0 / std::pow(ka, 3.0) * (1.5 + std::log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * std::pow(1.0 + ka, 2.0));
}

AssociationSetup build_association_setup(
    const vector<double> &x,
    const add_args &cppargs,
    const vector<double> &s_ij,
    const vector<double> &ghs,
    double t
) {
    int ncomp = static_cast<int>(x.size());
    AssociationSetup setup;
    setup.site_component_index.reserve(ncomp);
    setup.x_assoc.reserve(ncomp);

    for (std::vector<int>::const_iterator it = cppargs.assoc_num.begin(); it != cppargs.assoc_num.end(); ++it) {
        for (int i = 0; i < *it; ++i) {
            setup.site_component_index.push_back(static_cast<int>(it - cppargs.assoc_num.begin()));
            setup.x_assoc.push_back(x[setup.site_component_index.back()]);
        }
    }

    int num_sites = static_cast<int>(setup.site_component_index.size());
    setup.delta_ij.assign(static_cast<size_t>(num_sites * num_sites), 0.0);

    int idxa = 0;
    for (int i = 0; i < num_sites; ++i) {
        int comp_i = setup.site_component_index[i];
        int idxi = comp_i * ncomp + comp_i;
        for (int j = 0; j < num_sites; ++j) {
            int comp_j = setup.site_component_index[j];
            int idxj = comp_j * ncomp + comp_j;
            if (cppargs.assoc_matrix[idxa] != 0) {
                double eABij = 0.5 * (cppargs.e_assoc[comp_i] + cppargs.e_assoc[comp_j]);
                double volABij = sqrt(cppargs.vol_a[comp_i] * cppargs.vol_a[comp_j]) * std::pow(
                    sqrt(s_ij[idxi] * s_ij[idxj]) / (0.5 * (s_ij[idxi] + s_ij[idxj])),
                    3.0
                );
                if (!cppargs.k_hb.empty()) {
                    volABij *= (1.0 - cppargs.k_hb[comp_i * ncomp + comp_j]);
                }
                setup.delta_ij[idxa] = ghs[comp_i * ncomp + comp_j] * (std::exp(eABij / t) - 1.0)
                    * std::pow(s_ij[comp_i * ncomp + comp_j], 3.0) * volABij;
            }
            ++idxa;
        }
    }

    return setup;
}

int gcd_int(int a, int b) {
    a = std::abs(a);
    b = std::abs(b);
    while (b != 0) {
        int t = a % b;
        a = b;
        b = t;
    }
    return a == 0 ? 1 : a;
}

struct ChargeGroups {
    vector<int> cations;
    vector<int> anions;
    vector<int> solvents;
};

ChargeGroups collect_charge_groups(const add_args& args, size_t ncomp) {
    ChargeGroups groups;
    groups.cations.reserve(ncomp);
    groups.anions.reserve(ncomp);
    groups.solvents.reserve(ncomp);
    for (size_t i = 0; i < ncomp; ++i) {
        if (i >= args.z.size()) {
            throw ValueError("Composition and charge vectors must be aligned.");
        }
        if (std::abs(args.z[i]) < 1e-12) {
            groups.solvents.push_back(static_cast<int>(i));
        }
        else if (args.z[i] > 0.0) {
            groups.cations.push_back(static_cast<int>(i));
        }
        else {
            groups.anions.push_back(static_cast<int>(i));
        }
    }
    return groups;
}

double stable_logz_over_zminus1(double Z) {
    double dz = Z - 1.0;
    if (std::abs(dz) < 1e-8) {
        double dz2 = dz * dz;
        double dz3 = dz2 * dz;
        return 1.0 - 0.5 * dz + dz2 / 3.0 - 0.25 * dz3;
    }
    return std::log(Z) / dz;
}

vector<double> normalize_z_contributions(const vector<double> &z_raw, double Z_total) {
    vector<double> out = z_raw;
    double target = Z_total - 1.0;
    double raw_sum = 0.0;
    for (double value : z_raw) {
        raw_sum += value;
    }

    if (std::abs(target) <= 1e-12 && std::abs(raw_sum) <= 1e-12) {
        std::fill(out.begin(), out.end(), 0.0);
        return out;
    }
    if (std::abs(raw_sum) <= 1e-14) {
        throw ValueError("Could not normalize contribution Z terms because their sum is ~0 while Z-1 is non-zero.");
    }
    double scale = target / raw_sum;
    for (double &value : out) {
        value *= scale;
    }
    return out;
}

struct DensityScanPoint {
    double nu = 0.0;
    double rho = 0.0;
    double resid = 0.0;
    bool finite = false;
};

struct DensityBracket {
    double nu_lo = 0.0;
    double nu_hi = 0.0;
    double resid_lo = 0.0;
    double resid_hi = 0.0;
};

struct DensityRootCandidate {
    double rho_sort = 0.0;
    double rho = 0.0;
    double gres = 0.0;
    double rel_resid = 0.0;
    double abs_p_error = 0.0;
    double dpdrho = 0.0;
    bool valid = false;
};

struct FlashPhaseState {
    double p = 0.0;
    double rhol = 0.0;
    double rhov = 0.0;
    bool valid = false;
};

vector<double> build_density_scan_grid() {
    const double nu_min = 1e-13;
    const double nu_log_max = 5e-3;
    const double nu_linear_start = 5e-3;
    const double nu_max = 0.7405 - 1e-4;
    const int n_log = 24;
    const int n_linear = 256;

    vector<double> grid;
    grid.reserve(1 + n_log + n_linear);
    grid.push_back(nu_min);

    if (n_log > 0) {
        const double log_min = std::log(nu_min);
        const double log_max = std::log(nu_log_max);
        for (int i = 0; i < n_log; i++) {
            double frac = static_cast<double>(i) / static_cast<double>(n_log - 1);
            double nu = std::exp(log_min + frac * (log_max - log_min));
            if (nu > grid.back()) {
                grid.push_back(nu);
            }
        }
    }

    if (n_linear > 0) {
        for (int i = 0; i < n_linear; i++) {
            double frac = static_cast<double>(i) / static_cast<double>(n_linear - 1);
            double nu = nu_linear_start + frac * (nu_max - nu_linear_start);
            if (nu > grid.back()) {
                grid.push_back(nu);
            }
        }
    }

    return grid;
}

DensityScanPoint evaluate_density_scan_point(double nu, double t, int ncomp, const vector<double> &x, double p, const add_args &cppargs) {
    DensityScanPoint point;
    point.nu = nu;
    point.rho = reduced_to_molar(nu, t, ncomp, x, cppargs);
    try {
        point.resid = resid_rho(point.rho, t, p, x, cppargs);
        point.finite = std::isfinite(point.resid);
    }
    catch (const std::exception&) {
        point.resid = 0.0;
        point.finite = false;
    }
    return point;
}

vector<DensityBracket> find_density_brackets(const vector<DensityScanPoint> &points) {
    vector<DensityBracket> brackets;
    if (points.size() < 2) {
        return brackets;
    }

    for (size_t i = 1; i < points.size(); i++) {
        const DensityScanPoint &lo = points[i - 1];
        const DensityScanPoint &hi = points[i];
        if (!lo.finite || !hi.finite) {
            continue;
        }
        if (lo.resid * hi.resid < 0.0) {
            DensityBracket bracket;
            bracket.nu_lo = lo.nu;
            bracket.nu_hi = hi.nu;
            bracket.resid_lo = lo.resid;
            bracket.resid_hi = hi.resid;
            brackets.push_back(bracket);
        }
    }

    return brackets;
}

void append_refined_density_brackets(
    const DensityBracket &coarse,
    double t,
    int ncomp,
    const vector<double> &x,
    double p,
    const add_args &cppargs,
    vector<DensityBracket> &refined_brackets
) {
    const int refine_segments = 256;
    const double discontinuity_threshold = 1e12;

    vector<DensityScanPoint> refined_points;
    refined_points.reserve(refine_segments + 1);
    for (int i = 0; i <= refine_segments; i++) {
        double frac = static_cast<double>(i) / static_cast<double>(refine_segments);
        double nu = coarse.nu_lo + frac * (coarse.nu_hi - coarse.nu_lo);
        refined_points.push_back(evaluate_density_scan_point(nu, t, ncomp, x, p, cppargs));
    }

    vector<DensityBracket> local_brackets = find_density_brackets(refined_points);
    for (const DensityBracket &candidate : local_brackets) {
        double nu_mid = 0.5 * (candidate.nu_lo + candidate.nu_hi);
        DensityScanPoint mid = evaluate_density_scan_point(nu_mid, t, ncomp, x, p, cppargs);
        if (!mid.finite) {
            continue;
        }

        double min_abs = std::min(std::abs(candidate.resid_lo), std::abs(candidate.resid_hi));
        min_abs = std::min(min_abs, std::abs(mid.resid));
        if (min_abs > discontinuity_threshold) {
            continue;
        }

        refined_brackets.push_back(candidate);
    }
}

bool validate_density_root(
    double t,
    double p,
    const vector<double> &x,
    const add_args &cppargs,
    double rho,
    DensityRootCandidate *candidate
) {
    const double rel_tol = 1e-6;
    const double abs_tol = 1e-7;
    const double ultra_low_pressure_rel_tol = 2.5e-4;
    const double ultra_low_pressure_cutoff = 1e-3;

    double p_calc = 0.0;
    try {
        p_calc = pcsaft_p_cpp(t, rho, x, cppargs);
    }
    catch (const std::exception&) {
        return false;
    }
    if (!std::isfinite(p_calc)) {
        return false;
    }

    double abs_p_error = std::abs(p_calc - p);
    double rel_resid = abs_p_error / std::max(std::abs(p), 1e-300);
    bool pressure_ok = (rel_resid <= rel_tol || abs_p_error <= abs_tol);
    if (!pressure_ok && p <= ultra_low_pressure_cutoff) {
        pressure_ok = (rel_resid <= ultra_low_pressure_rel_tol);
    }
    if (!pressure_ok) {
        return false;
    }

    double h = std::max(1e-12, std::abs(rho) * 1e-6);
    double dpdrho = 0.0;
    try {
        double p_plus = pcsaft_p_cpp(t, rho + h, x, cppargs);
        if (!std::isfinite(p_plus)) {
            return false;
        }
        if (rho > h) {
            double p_minus = pcsaft_p_cpp(t, rho - h, x, cppargs);
            if (!std::isfinite(p_minus)) {
                return false;
            }
            dpdrho = (p_plus - p_minus) / (2.0 * h);
        }
        else {
            dpdrho = (p_plus - p_calc) / h;
        }
    }
    catch (const std::exception&) {
        return false;
    }
    if (!(std::isfinite(dpdrho) && dpdrho > 0.0)) {
        return false;
    }

    double gres = 0.0;
    try {
        gres = pcsaft_gres_cpp(t, rho, x, cppargs);
    }
    catch (const std::exception&) {
        return false;
    }
    if (!std::isfinite(gres)) {
        return false;
    }

    candidate->rho = rho;
    candidate->gres = gres;
    candidate->rel_resid = rel_resid;
    candidate->abs_p_error = abs_p_error;
    candidate->dpdrho = dpdrho;
    candidate->valid = true;
    return true;
}

bool try_flash_phase_state(
    double t,
    double p,
    const vector<double> &xl,
    const vector<double> &xv,
    const add_args &cppargs,
    FlashPhaseState *state
) {
    try {
        state->p = p;
        state->rhol = pcsaft_den_cpp(t, p, xl, 0, cppargs);
        state->rhov = pcsaft_den_cpp(t, p, xv, 1, cppargs);
    }
    catch (const std::exception&) {
        state->valid = false;
        return false;
    }

    state->valid = std::isfinite(state->rhol) && std::isfinite(state->rhov) && ((state->rhol - state->rhov) >= 1e-4);
    return state->valid;
}

FlashPhaseState find_nearby_flash_phase_state(
    double t,
    double base_p,
    const vector<double> &xl,
    const vector<double> &xv,
    const add_args &cppargs,
    int preferred_direction
) {
    FlashPhaseState state;
    const std::array<double, 10> relative_steps = {1e-2, 5e-3, 2e-3, 1e-3, 5e-4, 2e-4, 1e-4, 5e-5, 2e-5, 1e-5};
    const std::array<int, 2> directions = {preferred_direction >= 0 ? 1 : -1, preferred_direction >= 0 ? -1 : 1};

    for (int direction : directions) {
        for (double rel_step : relative_steps) {
            double candidate_p = base_p * (1.0 + direction * rel_step);
            if (!(candidate_p > 0.0) || std::abs(candidate_p - base_p) <= DBL_EPSILON * std::max(1.0, std::abs(base_p))) {
                continue;
            }
            if (try_flash_phase_state(t, candidate_p, xl, xv, cppargs, &state)) {
                return state;
            }
        }
    }

    state.valid = false;
    return state;
}

bool is_effectively_pure_feed(const vector<double> &x) {
    int nonzero = 0;
    for (double xi : x) {
        if (std::abs(xi) > 1e-12) {
            nonzero += 1;
        }
    }
    return nonzero == 1;
}

double pure_component_flash_resid(double t, double p, const vector<double> &x, const add_args &cppargs) {
    FlashPhaseState state;
    if (!try_flash_phase_state(t, p, x, x, cppargs, &state)) {
        throw SolutionError("Pure-component flash residual could not evaluate both phase densities.");
    }

    vector<double> fugcoef_l = pcsaft_fugcoef_cpp(t, state.rhol, x, cppargs);
    vector<double> fugcoef_v = pcsaft_fugcoef_cpp(t, state.rhov, x, cppargs);
    int comp = 0;
    for (int i = 0; i < static_cast<int>(x.size()); ++i) {
        if (std::abs(x[i]) > 1e-12) {
            comp = i;
            break;
        }
    }

    double phi_l = fugcoef_l[comp];
    double phi_v = fugcoef_v[comp];
    if (!(phi_l > 0.0) || !(phi_v > 0.0) || !std::isfinite(phi_l) || !std::isfinite(phi_v)) {
        throw SolutionError("Pure-component flash residual encountered a non-positive fugacity coefficient.");
    }
    return std::log(phi_l) - std::log(phi_v);
}

vector<double> pure_component_flashTQ(double t, vector<double> x, const add_args &cppargs) {
    const double log_p_min = -6.0;
    const double log_p_max = 9.0;
    const int scan_points = 301;
    const int refine_levels = 6;
    const int refine_points = 121;
    const double resid_tol = 1e-8;
    vector<std::pair<double, double>> valid_points;
    valid_points.reserve(scan_points);

    for (int i = 0; i < scan_points; ++i) {
        double frac = static_cast<double>(i) / static_cast<double>(scan_points - 1);
        double log_p = log_p_min + frac * (log_p_max - log_p_min);
        double p = std::pow(10.0, log_p);

        double f = 0.0;
        try {
            f = pure_component_flash_resid(t, p, x, cppargs);
        }
        catch (const std::exception&) {
            continue;
        }
        if (!std::isfinite(f)) {
            continue;
        }
        valid_points.emplace_back(p, f);
    }

    if (valid_points.empty()) {
        throw SolutionError("No valid pure-component saturation states were found during the pressure scan.");
    }

    auto find_sign_change = [](const vector<std::pair<double, double>> &samples, double *p_lo, double *f_lo, double *p_hi, double *f_hi) {
        for (size_t i = 1; i < samples.size(); ++i) {
            const auto &left = samples[i - 1];
            const auto &right = samples[i];
            if (left.second == 0.0) {
                *p_lo = left.first;
                *f_lo = left.second;
                *p_hi = left.first;
                *f_hi = left.second;
                return true;
            }
            if (right.second == 0.0 || left.second * right.second < 0.0) {
                *p_lo = left.first;
                *f_lo = left.second;
                *p_hi = right.first;
                *f_hi = right.second;
                return true;
            }
        }
        return false;
    };

    double p_lo = 0.0, p_hi = 0.0, f_lo = 0.0, f_hi = 0.0;
    if (!find_sign_change(valid_points, &p_lo, &f_lo, &p_hi, &f_hi)) {
        throw SolutionError("A valid sign-change bracket could not be found for the pure-component flash.");
    }

    if (p_lo == p_hi || std::abs(f_lo) <= resid_tol) {
        vector<double> result;
        result.push_back(p_lo);
        result.insert(result.end(), x.begin(), x.end());
        result.insert(result.end(), x.begin(), x.end());
        return result;
    }

    for (int level = 0; level < refine_levels; ++level) {
        vector<std::pair<double, double>> refined;
        refined.reserve(refine_points);
        const double log_lo = std::log(p_lo);
        const double log_hi = std::log(p_hi);
        for (int i = 0; i < refine_points; ++i) {
            double frac = static_cast<double>(i) / static_cast<double>(refine_points - 1);
            double p = std::exp(log_lo + frac * (log_hi - log_lo));
            double f = 0.0;
            try {
                f = pure_component_flash_resid(t, p, x, cppargs);
            }
            catch (const std::exception&) {
                continue;
            }
            if (!std::isfinite(f)) {
                continue;
            }
            refined.emplace_back(p, f);
        }

        if (!refined.empty()) {
            double new_lo = 0.0, new_hi = 0.0, new_f_lo = 0.0, new_f_hi = 0.0;
            if (find_sign_change(refined, &new_lo, &new_f_lo, &new_hi, &new_f_hi)) {
                p_lo = new_lo;
                f_lo = new_f_lo;
                p_hi = new_hi;
                f_hi = new_f_hi;
                if (p_lo == p_hi || std::min(std::abs(f_lo), std::abs(f_hi)) <= resid_tol) {
                    break;
                }
            }
        }
    }

    double p_est = p_lo;
    double f_est = f_lo;
    if (p_hi != p_lo && f_hi != f_lo) {
        double log_lo = std::log(p_lo);
        double log_hi = std::log(p_hi);
        double log_est = log_lo - f_lo * (log_hi - log_lo) / (f_hi - f_lo);
        p_est = std::exp(log_est);
        try {
            f_est = pure_component_flash_resid(t, p_est, x, cppargs);
        }
        catch (const std::exception&) {
            if (std::abs(f_hi) < std::abs(f_lo)) {
                p_est = p_hi;
                f_est = f_hi;
            }
        }
    }
    if (std::abs(f_est) > std::min(std::abs(f_lo), std::abs(f_hi))) {
        if (std::abs(f_hi) < std::abs(f_lo)) {
            p_est = p_hi;
        }
        else {
            p_est = p_lo;
        }
    }

    vector<double> result;
    result.push_back(p_est);
    result.insert(result.end(), x.begin(), x.end());
    result.insert(result.end(), x.begin(), x.end());
    return result;
}
}

#if defined(HUGE_VAL) && !defined(_HUGE)
    # define _HUGE HUGE_VAL
#else
    // GCC Version of huge value macro
    #if defined(HUGE) && !defined(_HUGE)
    #  define _HUGE HUGE
    #endif
#endif

struct DebugFlagGuard {
    int &flag;
    int old;
    DebugFlagGuard(int &flag_ref, int new_value) : flag(flag_ref), old(flag_ref) { flag = new_value; }
    ~DebugFlagGuard() { flag = old; }
};

struct BornSSMDSData {
    vector<double> d_born;
    vector<double> D;
    vector<double> ddelta_prefac;
    vector<double> f_k;
    vector<double> bracket;
    double sum_bracket;
    double sum_invD;
    double sum_gap;
    double sum_dpref_over_D2;
};

struct DielcState {
    double eps;
    vector<double> deps_dx;
};

inline bool is_ion_species(const add_args &cppargs, int i) {
    return std::abs(cppargs.z[i]) > 1e-12;
}

double compute_ion_diameter(int i, double t, const add_args &cppargs) {
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
        return sigma_i*(1.0 - 0.12);
    }
    if (mode == 2) {
        return sigma_i*(1.0 - 0.12*std::exp(-3.0*cppargs.e[i]/t));
    }
    throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
}

double compute_ion_diameter_dt(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return 0.0;
    }
    int mode = cppargs.d_ion_mode;
    if (mode == 2) {
        double sigma_i = cppargs.s[i];
        double expo = std::exp(-3.0*cppargs.e[i]/t);
        return -0.36*sigma_i*cppargs.e[i]*expo/(t*t);
    }
    return 0.0;
}

double compute_ion_born_radius(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return cppargs.s[i];
    }
    int mode = cppargs.d_born_mode;
    double sigma_i = cppargs.s[i];
    if (sigma_i <= 0.0) {
        throw ValueError("Born term requires positive ionic sigma_i.");
    }
    if (mode == 0) {
        return sigma_i;
    }
    if (mode == 1) {
        return sigma_i*(1.0 - 0.12);
    }
    if (mode == 2) {
        return sigma_i*(1.0 - 0.12*std::exp(-3.0*cppargs.e[i]/t));
    }
    if (mode == 3) {
        if (cppargs.d_born.size() <= static_cast<size_t>(i) || cppargs.d_born[i] <= 0.0) {
            throw ValueError("d_Born_mode=fitted_param requires positive ionic params['d_born'] values.");
        }
        return cppargs.d_born[i];
    }
    throw ValueError("Unknown d_Born_mode. Supported values are 0, 1, 2, 3.");
}

double compute_ion_born_radius_dt(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return 0.0;
    }
    int mode = cppargs.d_born_mode;
    if (mode == 2) {
        double sigma_i = cppargs.s[i];
        double expo = std::exp(-3.0*cppargs.e[i]/t);
        return -0.36*sigma_i*cppargs.e[i]*expo/(t*t);
    }
    return 0.0;
}

double compute_eps_rule(int rule, const vector<double> &x, const add_args &cppargs);
vector<double> compute_deps_rule_fd(int rule, const vector<double> &x, const add_args &cppargs);
DielcState evaluate_dielc_state(const vector<double> &x, const add_args &cppargs);

double compute_eps_aqueous_organic_mixed(const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule=8 requires params['z'] as an array with length equal to ncomp.");
    }
    if (cppargs.mixed_rel_perm_a.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_b.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_c.size() != static_cast<size_t>(ncomp) ||
        cppargs.mixed_rel_perm_mask.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule=8 requires mixed relative-permittivity arrays with length equal to ncomp.");
    }

    double x_sol = 0.0;
    double x_water = 0.0;
    double x_org = 0.0;
    double eps_org_num = 0.0;
    double a_num = 0.0;
    double b_num = 0.0;
    double c_num = 0.0;
    bool needs_coeffs = false;

    int water_idx = cppargs.mixed_rel_perm_water_index;
    bool has_water_component = (
        water_idx >= 0 &&
        water_idx < ncomp &&
        std::abs(cppargs.z[water_idx]) <= 1e-12
    );

    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) > 1e-12) {
            continue;
        }
        x_sol += x[i];
        if (has_water_component && i == water_idx) {
            x_water += x[i];
            continue;
        }
        x_org += x[i];
        eps_org_num += x[i] * cppargs.dielc[i];
        if (x[i] > 0.0) {
            needs_coeffs = true;
            if (cppargs.mixed_rel_perm_mask[i] == 0) {
                if (x_water > 0.0 || has_water_component) {
                    throw ValueError("dielc_rule=8 is missing mixed relative-permittivity coefficients for an organic solvent.");
                }
            } else {
                a_num += x[i] * cppargs.mixed_rel_perm_a[i];
                b_num += x[i] * cppargs.mixed_rel_perm_b[i];
                c_num += x[i] * cppargs.mixed_rel_perm_c[i];
            }
        }
    }

    if (x_sol <= 0.0) {
        throw ValueError("dielc_rule=8 requires at least one solvent species (z=0).");
    }
    if (x_org <= DBL_EPSILON) {
        if (!has_water_component) {
            throw ValueError("dielc_rule=8 requires at least one organic solvent or water component.");
        }
        return cppargs.dielc[water_idx];
    }

    double eps_org = eps_org_num / x_org;
    if (!has_water_component || x_water <= DBL_EPSILON) {
        return eps_org;
    }
    if (x_water >= x_sol - DBL_EPSILON) {
        return cppargs.dielc[water_idx];
    }
    if (!needs_coeffs) {
        throw ValueError("dielc_rule=8 requires mixed relative-permittivity coefficients for the organic solvent phase.");
    }

    double xw_sf = x_water / x_sol;
    double a_eff = a_num / x_org;
    double b_eff = b_num / x_org;
    double c_eff = c_num / x_org;
    return eps_org + ((a_eff * xw_sf + b_eff) * xw_sf + c_eff) * xw_sf;
}

double compute_eps_solvent_reference(const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
        return compute_eps_rule(cppargs.dielc_rule, x, cppargs);
    }
    double x_sol = 0.0;
    double eps_sol_num = 0.0;
    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) <= 1e-12) {
            x_sol += x[i];
            eps_sol_num += x[i]*cppargs.dielc[i];
        }
    }
    if (x_sol <= 0.0) {
        return compute_eps_rule(cppargs.dielc_rule, x, cppargs);
    }
    return eps_sol_num/x_sol;
}

vector<double> compute_deps_solvent_reference(const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    vector<double> deps(ncomp, 0.0);
    if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
        return deps;
    }
    double x_sol = 0.0;
    double eps_sol_num = 0.0;
    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) <= 1e-12) {
            x_sol += x[i];
            eps_sol_num += x[i]*cppargs.dielc[i];
        }
    }
    if (x_sol <= 0.0) {
        return deps;
    }
    double inv_xsol2 = 1.0/(x_sol*x_sol);
    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) <= 1e-12) {
            deps[i] = (cppargs.dielc[i]*x_sol - eps_sol_num)*inv_xsol2;
        }
    }
    return deps;
}

BornSSMDSData build_born_ssmds_data(vector<double> x, const add_args &cppargs, double t, double eps_r, double eps_r_ion) {
    int ncomp = static_cast<int>(x.size());
    const bool use_ssm = (cppargs.born_solvation_shell_model != 0);
    const bool use_ds = (cppargs.born_dielectric_saturation != 0);

    BornSSMDSData data;
    data.d_born.assign(ncomp, 1.0);
    data.D.assign(ncomp, 1.0);
    data.ddelta_prefac.assign(ncomp, 0.0);
    data.f_k.assign(ncomp, 1.0);
    data.bracket.assign(ncomp, 0.0);
    data.sum_bracket = 0.0;
    data.sum_invD = 0.0;
    data.sum_gap = 0.0;
    data.sum_dpref_over_D2 = 0.0;

    double f_mix = 0.0;
    for (int i = 0; i < ncomp; i++) {
        bool is_ion = is_ion_species(cppargs, i);
        double fi = 1.0;
        if (!is_ion && cppargs.f_solv.size() > static_cast<size_t>(i)) {
            fi = cppargs.f_solv[i];
        }
        data.f_k[i] = fi;
        f_mix += x[i]*fi;

        if (is_ion) {
            data.d_born[i] = compute_ion_born_radius(i, t, cppargs);
        }
        else if (cppargs.d_born.size() > static_cast<size_t>(i) && cppargs.d_born[i] > 0.0) {
            data.d_born[i] = cppargs.d_born[i];
        }
        else if (cppargs.s[i] > 0.0) {
            data.d_born[i] = cppargs.s[i];
        }
        else {
            throw ValueError("Born model requires positive solvent diameter.");
        }

        if (is_ion) {
            data.ddelta_prefac[i] = data.d_born[i]/std::abs(cppargs.z[i]);
        }
    }

    for (int i = 0; i < ncomp; i++) {
        bool is_ion = std::abs(cppargs.z[i]) > 1e-12;
        if (!is_ion) {
            data.D[i] = data.d_born[i];
            continue;
        }

        double delta_di = use_ssm ? ((f_mix - 1.0)*data.ddelta_prefac[i]) : 0.0;
        data.D[i] = data.d_born[i] + delta_di;
        if (data.D[i] <= 0.0) {
            throw ValueError("Born model generated a non-positive d_born + Delta d.");
        }

        double z2 = cppargs.z[i]*cppargs.z[i];
        double invD = 1.0/data.D[i];
        double gap = (1.0/data.d_born[i] - invD);
        double base_term = (1.0 - 1.0/eps_r)*invD;
        double ds_term = use_ds ? ((1.0 - 1.0/eps_r_ion)*gap) : 0.0;

        data.bracket[i] = base_term + ds_term;
        data.sum_bracket += x[i]*z2*data.bracket[i];
        data.sum_invD += x[i]*z2*invD;
        data.sum_gap += x[i]*z2*gap;
        if (use_ssm) {
            data.sum_dpref_over_D2 += x[i]*z2*data.ddelta_prefac[i]*invD*invD;
        }
    }
    return data;
}

double compute_born_ares_only(double t, const vector<double> &x, const add_args &cppargs) {
    if (cppargs.born_model == 0) {
        return 0.0;
    }
    double eps_mix = compute_eps_rule(cppargs.dielc_rule, x, cppargs);
    double eps_born = (cppargs.born_eps_mode == 1) ? compute_eps_solvent_reference(x, cppargs) : eps_mix;
    if (cppargs.born_model == 1) {
        double born_sum = 0.0;
        for (int i = 0; i < static_cast<int>(x.size()); i++) {
            if (is_ion_species(cppargs, i)) {
                double d_born_i = compute_ion_born_radius(i, t, cppargs);
                born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
            }
        }
        return -E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac)*(1.0 - 1.0/eps_born)*born_sum;
    }
    if (cppargs.born_model == 2) {
        const double eps_r_ion = 8.0;
        const double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
        BornSSMDSData born = build_born_ssmds_data(x, cppargs, t, eps_born, eps_r_ion);
        return -Kborn*born.sum_bracket;
    }
    throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
}

vector<double> compute_born_dadx_fd(double t, const vector<double> &x, const add_args &cppargs, double a0) {
    int ncomp = static_cast<int>(x.size());
    vector<double> dadx_born(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++) {
        double h = 1e-6*std::max(1.0, std::abs(x[i]));
        vector<double> xp = x;
        xp[i] += h;
        double fp = compute_born_ares_only(t, xp, cppargs);
        if (x[i] - h >= 0.0) {
            vector<double> xm = x;
            xm[i] -= h;
            double fm = compute_born_ares_only(t, xm, cppargs);
            dadx_born[i] = (fp - fm)/(2.0*h);
        }
        else {
            dadx_born[i] = (fp - a0)/h;
        }
        if (!std::isfinite(dadx_born[i])) {
            throw ValueError("Non-finite Born finite-difference derivative.");
        }
    }
    return dadx_born;
}

double compute_dh_ares_only(double t, double rho, const vector<double> &x, const add_args &cppargs) {
    if (cppargs.z.empty()) {
        return 0.0;
    }
    int ncomp = static_cast<int>(x.size());
    vector<double> d(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++) {
        d[i] = cppargs.s[i]*(1.0 - 0.12*std::exp(-3.0*cppargs.e[i]/t));
        if (is_ion_species(cppargs, i)) {
            d[i] = compute_ion_diameter(i, t, cppargs);
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

    DielcState dielc_state = evaluate_dielc_state(x, cppargs);
    double eps = dielc_state.eps;
    double kappa = compute_dh_kappa(den, t, eps, Qsum);
    if (kappa == 0.0) {
        return 0.0;
    }

    double S = 0.0;
    for (int i = 0; i < ncomp; i++) {
        S += x[i]*cppargs.z[i]*cppargs.z[i]*compute_dh_chi(kappa, d[i]);
    }

    double K0 = E_CHRG*E_CHRG/(12.0*PI*kb*t*perm_vac);
    return -K0*kappa/eps*S;
}

vector<double> compute_dh_dadx_fd(double t, double rho, const vector<double> &x, const add_args &cppargs, double a0) {
    int ncomp = static_cast<int>(x.size());
    vector<double> dadx_dh(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++) {
        double h = 1e-6*std::max(1.0, std::abs(x[i]));
        vector<double> xp = x;
        xp[i] += h;
        double fp = compute_dh_ares_only(t, rho, xp, cppargs);
        if (x[i] - h >= 0.0) {
            vector<double> xm = x;
            xm[i] -= h;
            double fm = compute_dh_ares_only(t, rho, xm, cppargs);
            dadx_dh[i] = (fp - fm)/(2.0*h);
        }
        else {
            dadx_dh[i] = (fp - a0)/h;
        }
        if (!std::isfinite(dadx_dh[i])) {
            throw ValueError("Non-finite DH finite-difference derivative.");
        }
    }
    return dadx_dh;
}

double get_ares_contribution_value(const AresContributions &terms, AresContributionKind kind) {
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

AresContributions compute_ares_contributions_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs) {
    AresContributions out;
    int ncomp = static_cast<int>(x.size());
    ThermoCommonState thermo = build_thermo_common_state(t, rho, x, cppargs, false);
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
    double summ = 0.0;
    DispersionPolynomialState dispersion = build_dispersion_polynomials(m_avg, eta);
    double I1 = dispersion.I1;
    double I2 = dispersion.I2;
    double C1 = dispersion.C1;
    double ares_hs = 1.0 / zeta[0] * (3.0 * zeta[1] * zeta[2] / (1.0 - zeta[3])
            + std::pow(zeta[2], 3.0) / (zeta[3] * std::pow(1.0 - zeta[3], 2.0))
            + (std::pow(zeta[2], 3.0) / std::pow(zeta[3], 2.0) - zeta[0]) * std::log(1.0 - zeta[3]));

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i]*(cppargs.m[i]-1)*log(ghs[i*ncomp+i]);
    }
    out.hc = m_avg*ares_hs - summ;
    out.disp = -2*PI*den*I1*m2es3 - PI*den*m_avg*C1*I2*m2e2s3;

    if (!cppargs.dipm.empty()) {
        double A2 = 0.0;
        double A3 = 0.0;
        vector<double> dipmSQ(ncomp, 0.0);

        const static double conv = 7242.702976750923;
        for (int i = 0; i < ncomp; i++) {
            dipmSQ[i] = pow(cppargs.dipm[i], 2.)/(cppargs.m[i]*cppargs.e[i]*pow(cppargs.s[i],3.))*conv;
        }

        vector<double> adip(5, 0.0);
        vector<double> bdip(5, 0.0);
        vector<double> cdip(5, 0.0);
        double J2, J3, m_ij, m_ijk;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < ncomp; j++) {
                m_ij = sqrt(cppargs.m[i]*cppargs.m[j]);
                if (m_ij > 2) {
                    m_ij = 2;
                }
                vector<double> adip = build_dipole_coefficients(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = build_dipole_coefficients(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
                J2 = 0.0;
                for (int l = 0; l < 5; l++) {
                    J2 += (adip[l] + bdip[l]*e_ij[j*ncomp+j]/t)*pow(eta, l);
                }
                A2 += x[i]*x[j]*e_ij[i*ncomp+i]/t*e_ij[j*ncomp+j]/t*pow(s_ij[i*ncomp+i],3)*pow(s_ij[j*ncomp+j],3)/
                    pow(s_ij[i*ncomp+j],3)*cppargs.dip_num[i]*cppargs.dip_num[j]*dipmSQ[i]*dipmSQ[j]*J2;

                for (int k = 0; k < ncomp; k++) {
                    m_ijk = pow((cppargs.m[i]*cppargs.m[j]*cppargs.m[k]),1/3.);
                    if (m_ijk > 2) {
                        m_ijk = 2;
                    }
                    vector<double> cdip = build_dipole_coefficients(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
                    J3 = 0.0;
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
        int num_sites = 0;
        vector<int> iA;
        for (std::vector<int>::const_iterator it = cppargs.assoc_num.begin(); it != cppargs.assoc_num.end(); ++it) {
            num_sites += *it;
            for (int i = 0; i < *it; i++) {
                iA.push_back(static_cast<int>(it - cppargs.assoc_num.begin()));
            }
        }

        vector<double> x_assoc(num_sites);
        for (int i = 0; i < num_sites; i++) {
            x_assoc[i] = x[iA[i]];
        }

        vector<double> XA(num_sites, 0.0);
        vector<double> delta_ij(num_sites * num_sites, 0.0);
        int idxa = 0;
        int idxi = 0;
        int idxj = 0;
        for (int i = 0; i < num_sites; i++) {
            idxi = iA[i]*ncomp+iA[i];
            for (int j = 0; j < num_sites; j++) {
                idxj = iA[j]*ncomp+iA[j];
                if (cppargs.assoc_matrix[idxa] != 0) {
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
                    delta_ij[idxa] = ghs[iA[i]*ncomp+iA[j]]*(exp(eABij/t)-1)*pow(s_ij[iA[i]*ncomp+iA[j]], 3)*volABij;
                }
                idxa += 1;
            }
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        int ctr = 0;
        double dif = 1000.0;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-15)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.0;
            for (int i = 0; i < num_sites; i++) {
                dif += std::abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        out.assoc = 0.0;
        for (int i = 0; i < num_sites; i++) {
            out.assoc += x[iA[i]]*(log(XA[i])-0.5*XA[i] + 0.5);
        }
    }

    if (!cppargs.z.empty()) {
        DielcState dielc_state = evaluate_dielc_state(x, cppargs);
        double eps = dielc_state.eps;
        double eps_born = (cppargs.born_eps_mode == 1) ? compute_eps_solvent_reference(x, cppargs) : eps;
        out.ion = compute_dh_ares_only(t, rho, x, cppargs);

        if (cppargs.born_model == 1) {
            double born_sum = 0.0;
            for (int i = 0; i < ncomp; i++) {
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = compute_ion_born_radius(i, t, cppargs);
                    born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
                }
            }
            out.born = -E_CHRG*E_CHRG/(4.*PI*kb*t*perm_vac)*(1.-1./eps_born)*born_sum;
        }
        else if (cppargs.born_model == 2) {
            const double eps_r_ion = 8.0;
            const double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
            BornSSMDSData born = build_born_ssmds_data(x, cppargs, t, eps_born, eps_r_ion);
            out.born = -Kborn*born.sum_bracket;
        }
        else if (cppargs.born_model != 0) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }

    return out;
}

vector<double> compute_contribution_dadx_fd(AresContributionKind kind, double t, double rho, const vector<double> &x, const add_args &cppargs, double a0) {
    int ncomp = static_cast<int>(x.size());
    vector<double> dadx(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++) {
        double h = 1e-6*std::max(1.0, std::abs(x[i]));
        vector<double> xp = x;
        xp[i] += h;
        double fp = get_ares_contribution_value(compute_ares_contributions_cpp(t, rho, xp, cppargs), kind);
        if (x[i] - h >= 0.0) {
            vector<double> xm = x;
            xm[i] -= h;
            double fm = get_ares_contribution_value(compute_ares_contributions_cpp(t, rho, xm, cppargs), kind);
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

void validate_dielc_inputs(const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.dielc.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("params['dielc'] must be an array with length equal to ncomp.");
    }
    if (cppargs.dielc_diff_mode != 0 && cppargs.dielc_diff_mode != 1) {
        throw ValueError("Unknown dielc_diff_mode. Supported values are 0 (analytic) and 1 (finite-diff).");
    }
    if (cppargs.hc_dadx_diff_mode != 0 && cppargs.hc_dadx_diff_mode != 1) {
        throw ValueError("Unknown hc_model dadx_differential_mode. Supported values are analytical/numerical (0/1).");
    }
    if (cppargs.disp_dadx_diff_mode != 0 && cppargs.disp_dadx_diff_mode != 1) {
        throw ValueError("Unknown disp_model dadx_differential_mode. Supported values are analytical/numerical (0/1).");
    }
    if (cppargs.assoc_dadx_diff_mode != 0 && cppargs.assoc_dadx_diff_mode != 1) {
        throw ValueError("Unknown assoc_model dadx_differential_mode. Supported values are analytical/numerical (0/1).");
    }
    if (cppargs.polar_dadx_diff_mode != 0 && cppargs.polar_dadx_diff_mode != 1) {
        throw ValueError("Unknown polar_model dadx_differential_mode. Supported values are analytical/numerical (0/1).");
    }
    if (cppargs.born_diff_mode != 0 && cppargs.born_diff_mode != 1 && cppargs.born_diff_mode != 2 && cppargs.born_diff_mode != 3) {
        throw ValueError("Unknown born_diff_mode. Supported values are 0 (analytic), 1 (finite-diff), 2 (Eq.133-style), and 3 (no dielectric-concentration term).");
    }
    if (cppargs.d_ion_mode < 0 || cppargs.d_ion_mode > 2) {
        throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
    }
    if (cppargs.mu_DH_diff_mode != 0 && cppargs.mu_DH_diff_mode != 1) {
        throw ValueError("Unknown mu_DH differential_mode. Supported values are analytical/numerical (0/1).");
    }
    if (cppargs.mu_DH_comp_dep_rel_perm != 0 && cppargs.mu_DH_comp_dep_rel_perm != 1) {
        throw ValueError("mu_DH comp_dep_rel_perm must be 0 or 1.");
    }
    if (cppargs.mu_DH_include_sum_term != 0 && cppargs.mu_DH_include_sum_term != 1) {
        throw ValueError("mu_DH include_sum_term must be 0 or 1.");
    }
    if (cppargs.include_born_model != 0 && cppargs.include_born_model != 1) {
        throw ValueError("include_born_model must be 0 or 1.");
    }
    if (cppargs.d_born_mode < 0 || cppargs.d_born_mode > 3) {
        throw ValueError("Unknown d_Born_mode. Supported values are 0, 1, 2, 3.");
    }
    if (cppargs.born_bulk_mode != 0 && cppargs.born_bulk_mode != 1) {
        throw ValueError("Unknown born bulk_mode. Supported values are mix/solvent (0/1).");
    }
    if (cppargs.mu_born_diff_mode != 0 && cppargs.mu_born_diff_mode != 1) {
        throw ValueError("Unknown mu_born differential_mode. Supported values are analytical/numerical (0/1).");
    }
    if (cppargs.born_eps_mode != 0 && cppargs.born_eps_mode != 1) {
        throw ValueError("Unknown born_eps_mode. Supported values are 0 (eps_r,mix) and 1 (eps_r,solvent).");
    }
    if (cppargs.born_model < 0 || cppargs.born_model > 2) {
        throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
    }
    if (cppargs.born_model > 0 && cppargs.z.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("Born contribution requires params['z'] as an array with length equal to ncomp.");
    }
    if (cppargs.born_radius_model < 1 || cppargs.born_radius_model > 5) {
        throw ValueError("Unknown born_radius_model. Supported values are 1, 2, 3, 4, 5.");
    }
    if (cppargs.born_model == 1 && cppargs.born_radius_model == 5) {
        throw ValueError("born_model=1 supports born_radius_model values 1, 2, 3, 4.");
    }
    if (cppargs.born_model > 0 && (cppargs.born_radius_model == 4 || cppargs.born_radius_model == 5)) {
        if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
            throw ValueError("born_radius_model 4/5 requires ionic charge array params['z'] with length ncomp.");
        }
        for (int i = 0; i < ncomp; i++) {
            if (is_ion_species(cppargs, i) &&
                (cppargs.d_born.size() <= static_cast<size_t>(i) || cppargs.d_born[i] <= 0.0)) {
                throw ValueError("born_radius_model 4/5 requires positive ionic params['d_born'] values.");
            }
        }
    }
    int rule = cppargs.dielc_rule;
    if (rule < 0 || rule > 8) {
        throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8.");
    }
    if ((rule == 2 || rule == 3 || rule == 4 || rule == 5) &&
        cppargs.mw.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule requires params['MW'] as an array with length equal to ncomp.");
    }
    if ((rule == 3 || rule == 4 || rule == 5 || rule == 6) &&
        cppargs.z.size() != static_cast<size_t>(ncomp)) {
        throw ValueError("dielc_rule requires params['z'] as an array with length equal to ncomp.");
    }
    if (rule == 8) {
        if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
            throw ValueError("dielc_rule=8 requires params['z'] as an array with length equal to ncomp.");
        }
        if (cppargs.mixed_rel_perm_a.size() != static_cast<size_t>(ncomp) ||
            cppargs.mixed_rel_perm_b.size() != static_cast<size_t>(ncomp) ||
            cppargs.mixed_rel_perm_c.size() != static_cast<size_t>(ncomp) ||
            cppargs.mixed_rel_perm_mask.size() != static_cast<size_t>(ncomp)) {
            throw ValueError("dielc_rule=8 requires mixed relative-permittivity arrays with length equal to ncomp.");
        }
    }
}

double compute_eps_rule(int rule, const vector<double> &x, const add_args &cppargs) {
    const double alpha = 7.01;
    int ncomp = static_cast<int>(x.size());
    if (rule == 0) {
        return *std::max_element(cppargs.dielc.begin(), cppargs.dielc.end());
    }
    if (rule == 1) {
        double eps = 0.0;
        for (int i = 0; i < ncomp; i++) {
            eps += x[i]*cppargs.dielc[i];
        }
        return eps;
    }
    if (rule == 7) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=7 requires at least one solvent species (z=0).");
        }
        if (idx_ion.empty()) {
            throw ValueError("dielc_rule=7 requires at least one ionic species (z!=0).");
        }
        double x_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            x_sol += x[idx];
            eps_sol_num += x[idx] * cppargs.dielc[idx];
        }
        if (x_sol < 0.0 || x_sol > 1.0) {
            throw ValueError("dielc_rule=7 encountered invalid solvent mole fraction.");
        }
        double eps_sol = 0.0;
        if (x_sol > 1.0e-16) {
            eps_sol = eps_sol_num / x_sol;
        }
        else {
            for (int idx : idx_sol) eps_sol += cppargs.dielc[idx];
            eps_sol /= static_cast<double>(idx_sol.size());
        }
        double eps_salt = 0.0;
        for (int idx : idx_ion) eps_salt += cppargs.dielc[idx];
        eps_salt /= static_cast<double>(idx_ion.size());
        return eps_sol * x_sol + eps_salt * (1.0 - x_sol);
    }
    if (rule == 8) {
        return compute_eps_aqueous_organic_mixed(x, cppargs);
    }
    if (rule == 2) {
        double mw_bar = 0.0;
        double num = 0.0;
        for (int i = 0; i < ncomp; i++) {
            mw_bar += x[i]*cppargs.mw[i];
            num += x[i]*cppargs.mw[i]*cppargs.dielc[i];
        }
        if (mw_bar <= 0.0) {
            throw ValueError("Average molecular weight must be positive for dielc_rule=2.");
        }
        return num/mw_bar;
    }
    if (rule == 3) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=3 requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx]*cppargs.mw[idx];
            eps_sol_num += x[idx]*cppargs.mw[idx]*cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.");
        }
        double eps_sol_w = eps_sol_num/mw_sol;
        double x_sol = 0.0;
        double eps_ion = 0.0;
        for (int idx : idx_sol) x_sol += x[idx];
        for (int idx : idx_ion) eps_ion += x[idx]*cppargs.dielc[idx];
        return x_sol*eps_sol_w + eps_ion;
    }
    if (rule == 4 || rule == 5) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sf_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx]*cppargs.mw[idx];
            eps_sf_num += x[idx]*cppargs.mw[idx]*cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule.");
        }
        double eps_sf = eps_sf_num/mw_sol;
        double x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        return eps_sf/(1.0 + alpha*x_ion);
    }
    if (rule == 6) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=6 requires at least one solvent species (z=0).");
        }
        double eps_sf_const = 0.0;
        for (int idx : idx_sol) eps_sf_const += cppargs.dielc[idx];
        eps_sf_const /= static_cast<double>(idx_sol.size());
        double x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        return eps_sf_const/(1.0 + alpha*x_ion);
    }
    throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8.");
}

vector<double> compute_deps_rule_analytic(int rule, const vector<double> &x, const add_args &cppargs) {
    const double alpha = 7.01;
    int ncomp = static_cast<int>(x.size());
    vector<double> deps_dx(ncomp, 0.0);
    if (rule == 0) {
        return deps_dx;
    }
    if (rule == 1) {
        return cppargs.dielc;
    }
    if (rule == 7) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=7 requires at least one solvent species (z=0).");
        }
        if (idx_ion.empty()) {
            throw ValueError("dielc_rule=7 requires at least one ionic species (z!=0).");
        }
        double x_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            x_sol += x[idx];
            eps_sol_num += x[idx] * cppargs.dielc[idx];
        }
        double eps_sol = 0.0;
        if (x_sol > 1.0e-16) {
            eps_sol = eps_sol_num / x_sol;
        }
        else {
            for (int idx : idx_sol) eps_sol += cppargs.dielc[idx];
            eps_sol /= static_cast<double>(idx_sol.size());
        }
        double eps_salt = 0.0;
        for (int idx : idx_ion) eps_salt += cppargs.dielc[idx];
        eps_salt /= static_cast<double>(idx_ion.size());
        vector<double> deps_dx(ncomp, 0.0);
        for (int idx : idx_sol) deps_dx[idx] = eps_sol;
        for (int idx : idx_ion) deps_dx[idx] = eps_salt;
        return deps_dx;
    }
    if (rule == 8) {
        return compute_deps_rule_fd(rule, x, cppargs);
    }
    if (rule == 2) {
        double mw_bar = 0.0;
        double eps_mix_num = 0.0;
        for (int i = 0; i < ncomp; i++) {
            mw_bar += x[i]*cppargs.mw[i];
            eps_mix_num += x[i]*cppargs.mw[i]*cppargs.dielc[i];
        }
        if (mw_bar <= 0.0) {
            throw ValueError("Average molecular weight must be positive for dielc_rule=2.");
        }
        double eps_mix = eps_mix_num/mw_bar;
        for (int i = 0; i < ncomp; i++) {
            deps_dx[i] = (cppargs.mw[i]/mw_bar)*(cppargs.dielc[i] - eps_mix);
        }
        return deps_dx;
    }
    if (rule == 3) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=3 requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sol_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx]*cppargs.mw[idx];
            eps_sol_num += x[idx]*cppargs.mw[idx]*cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=3.");
        }
        double eps_sol_w = eps_sol_num/mw_sol;
        double x_sol = 0.0;
        for (int idx : idx_sol) x_sol += x[idx];
        for (int idx : idx_sol) {
            deps_dx[idx] = eps_sol_w + x_sol*(cppargs.mw[idx]/mw_sol)*(cppargs.dielc[idx] - eps_sol_w);
        }
        for (int idx : idx_ion) {
            deps_dx[idx] = cppargs.dielc[idx];
        }
        return deps_dx;
    }
    if (rule == 4) {
        vector<int> idx_sol;
        vector<int> idx_ion;
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) idx_sol.push_back(i);
            else idx_ion.push_back(i);
        }
        if (idx_sol.empty()) {
            throw ValueError("dielc_rule=4 requires at least one solvent species (z=0).");
        }
        double mw_sol = 0.0;
        double eps_sf_num = 0.0;
        for (int idx : idx_sol) {
            mw_sol += x[idx]*cppargs.mw[idx];
            eps_sf_num += x[idx]*cppargs.mw[idx]*cppargs.dielc[idx];
        }
        if (mw_sol <= 0.0) {
            throw ValueError("Solvent molecular-weight denominator must be positive for dielc_rule=4.");
        }
        double eps_sf = eps_sf_num/mw_sol;
        double x_ion = 0.0;
        for (int idx : idx_ion) x_ion += x[idx];
        double den = 1.0 + alpha*x_ion;
        for (int idx : idx_sol) {
            deps_dx[idx] = (1.0/den)*(cppargs.mw[idx]/mw_sol)*(cppargs.dielc[idx] - eps_sf);
        }
        for (int idx : idx_ion) {
            deps_dx[idx] = -alpha*eps_sf/(den*den);
        }
        return deps_dx;
    }
    if (rule == 5) {
        return cppargs.dielc;
    }
    if (rule == 6) {
        for (int i = 0; i < ncomp; i++) {
            if (std::abs(cppargs.z[i]) <= 1e-12) deps_dx[i] = cppargs.dielc[i];
            else deps_dx[i] = 0.0;
        }
        return deps_dx;
    }
    throw ValueError("Unknown dielc_rule. Supported rules are 0, 1, 2, 3, 4, 5, 6, 7, 8.");
}

vector<double> compute_deps_rule_fd(int rule, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    vector<double> deps_dx(ncomp, 0.0);
    double f0 = compute_eps_rule(rule, x, cppargs);
    for (int i = 0; i < ncomp; i++) {
        double h = 1e-6*std::max(1.0, std::abs(x[i]));
        vector<double> xp = x;
        xp[i] += h;
        double fp = compute_eps_rule(rule, xp, cppargs);
        if (x[i] - h >= 0.0) {
            vector<double> xm = x;
            xm[i] -= h;
            double fm = compute_eps_rule(rule, xm, cppargs);
            deps_dx[i] = (fp - fm)/(2.0*h);
        }
        else {
            deps_dx[i] = (fp - f0)/h;
        }
        if (!std::isfinite(deps_dx[i])) {
            throw ValueError("Non-finite dielectric finite-difference derivative.");
        }
    }
    return deps_dx;
}

DielcState evaluate_dielc_state(const vector<double> &x, const add_args &cppargs) {
    validate_dielc_inputs(x, cppargs);
    DielcState state;
    state.eps = compute_eps_rule(cppargs.dielc_rule, x, cppargs);
    if (cppargs.dielc_diff_mode == 0 && cppargs.dielc_rule != 8) {
        state.deps_dx = compute_deps_rule_analytic(cppargs.dielc_rule, x, cppargs);
    }
    else {
        state.deps_dx = compute_deps_rule_fd(cppargs.dielc_rule, x, cppargs);
    }
    return state;
}


vector<double> XA_find(vector<double> XA_guess, vector<double> delta_ij, double den,
    vector<double> x) {
    /**Iterate over this function in order to solve for XA*/
    int num_sites = static_cast<int>(XA_guess.size());
    vector<double> XA = XA_guess;

    int idxij = -1; // index for delta_ij
    for (int i = 0; i < num_sites; i++) {
        double summ = 0.;
        for (int j = 0; j < num_sites; j++) {
            idxij += 1;
            summ += den*x[j]*XA_guess[j]*delta_ij[idxij];
        }
        XA[i] = 1./(1.+summ);
    }

    return XA;
}


vector<double> dXAdt_find(vector<double> delta_ij, double den,
    vector<double> XA, vector<double> ddelta_dt, vector<double> x) {
    /**Solve for the derivative of XA with respect to temperature.*/
    int num_sites = static_cast<int>(XA.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_sites, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites, num_sites);

    double summ;
    int ij = 0;
    for (int i = 0; i < num_sites; i++) {
        summ = 0;
        for (int j = 0; j < num_sites; j++) {
            B(i) -= x[j]*XA[j]*ddelta_dt[ij];
            A(i,j) = x[j]*delta_ij[ij];
            summ += x[j]*XA[j]*delta_ij[ij];
            ij += 1;
        }
        A(i,i) = pow(1+den*summ, 2.)/den;
    }

    Eigen::MatrixXd solution = A.lu().solve(B); //Solves linear system of equations
    vector<double> dXA_dt(num_sites);
    for (int i = 0; i < num_sites; i++) {
        dXA_dt[i] = solution(i);
    }
    return dXA_dt;
}


vector<double> dXAdx_find(vector<int> assoc_num, vector<double> delta_ij,
    double den, vector<double> XA, vector<double> ddelta_dx, vector<double> x) {
    /**Solve for the derivative of XA with respect to composition, or actually
    rho_i (the molar density of component i, which equals x_i * rho).*/
    int num_sites = static_cast<int>(XA.size());
    int ncomp = static_cast<int>(assoc_num.size());
    Eigen::MatrixXd B(num_sites*ncomp, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites*ncomp, num_sites*ncomp);

    double sum1, sum2;
    int idx1 = 0;
    int ij = 0;
    for (int i = 0; i < ncomp; i++) {
        for (int j = 0; j < num_sites; j++) {
            sum1 = 0;
            for (int k = 0; k < num_sites; k++) {
                sum1 = sum1 + den*x[k]*(XA[k]*ddelta_dx[i*num_sites*num_sites + j*num_sites + k]);
                A(ij,i*num_sites+k) = XA[j]*XA[j]*den*x[k]*delta_ij[j*num_sites+k];
            }

            sum2 = 0;
            for (int l = 0; l < assoc_num[i]; l++) {
                sum2 = sum2 + XA[idx1+l]*delta_ij[idx1*num_sites+l*num_sites+j];
            }

            A(ij,ij) = A(ij,ij) + 1;
            B(ij) = -1*XA[j]*XA[j]*(sum1 + sum2);
            ij += 1;
        }
        idx1 += assoc_num[i];
    }

    Eigen::MatrixXd solution = A.lu().solve(B); //Solves linear system of equations
    vector<double> dXA_dx(num_sites*ncomp);
    for (int i = 0; i < num_sites*ncomp; i++) {
        dXA_dx[i] = solution(i);
    }
    return dXA_dx;
}


double pcsaft_Z_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate the compressibility factor.
    */
    int ncomp = static_cast<int>(x.size()); // number of components
    ThermoCommonState thermo = build_thermo_common_state(t, rho, x, cppargs, false);
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
    DispersionPolynomialState dispersion = build_dispersion_polynomials(m_avg, eta);
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

    double Zid = 1.0;
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
                vector<double> adip = build_dipole_coefficients(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = build_dipole_coefficients(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
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
                    vector<double> cdip = build_dipole_coefficients(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
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
        AssociationSetup assoc = build_association_setup(x, cppargs, s_ij, ghs, t);
        const vector<int> &iA = assoc.site_component_index;
        const vector<double> &x_assoc = assoc.x_assoc;
        const vector<double> &delta_ij = assoc.delta_ij;
        int num_sites = static_cast<int>(iA.size());

        vector<double> XA(num_sites, 0.0);
        for (int i = 0; i < num_sites; i++) {
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        vector<double> ddelta_dx(num_sites * num_sites * ncomp, 0);
        int idx_ddelta = 0;
        for (int k = 0; k < ncomp; k++) {
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

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-15)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < num_sites; i++) {
                dif += std::abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        vector<double> dXA_dx(num_sites*ncomp, 0);
        dXA_dx = dXAdx_find(cppargs.assoc_num, delta_ij, den, XA, ddelta_dx, x_assoc);

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
        DielcState dielc_state = evaluate_dielc_state(x, cppargs);
        double eps = dielc_state.eps;
        double eps_born = (cppargs.born_eps_mode == 1) ? compute_eps_solvent_reference(x, cppargs) : eps;
        vector<double> q(cppargs.z.begin(), cppargs.z.end());
        for (int i = 0; i < ncomp; i++) {
            q[i] = q[i]*E_CHRG;
        }

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            summ += cppargs.z[i]*cppargs.z[i]*x[i];
        }

        double kappa = compute_dh_kappa(den, t, eps, summ); // the inverse Debye screening length. Equation 4 in Held et al. 2008.

        if (kappa != 0) {
            double chi, sigma_k;
            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi = compute_dh_chi(kappa, d[i]);
                sigma_k = -2*chi+3/(1+kappa*d[i]);
                summ += q[i]*q[i]*x[i]*sigma_k;
            }
            Zion = -1*kappa/24./PI/kb/t/(eps*perm_vac)*summ;
        }
    }

    // Born term (Bulow 2021a, non-SSM+DS) has no explicit density dependence.
    double Z = Zid + Zhc + Zdisp + Zpolar + Zassoc + Zion + Zborn;
    if (cppargs.debug) {
        std::cout << std::fixed << std::setprecision(10)
                  << "[DEBUG Zres] t=" << t << " rho=" << rho
                  << " Zhc=" << Zhc
                  << " Zdisp=" << Zdisp
                  << " Zpolar=" << Zpolar
                  << " Zassoc=" << Zassoc
                  << " Zion=" << Zion
                  << " Zborn=" << Zborn
                  << " Zres=" << (Z - 1.0)
                  << " Z=" << Z << std::endl;
    }
    return Z;
}


vector<double> pcsaft_lnfug_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate the natural logarithm of the fugacity coefficients for one phase of the system.
    */
    int ncomp = static_cast<int>(x.size()); // number of components
    vector<double> d (ncomp);
    for (int i = 0; i < ncomp; i++) {
        d[i] = cppargs.s[i]*(1-0.12*exp(-3*cppargs.e[i]/t));
        if (!cppargs.z.empty() && is_ion_species(cppargs, i)) {
            d[i] = compute_ion_diameter(i, t, cppargs);
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
    double Zhs = zeta[3]/(1-zeta[3]) + 3.*zeta[1]*zeta[2]/zeta[0]/(1.-zeta[3])/(1.-zeta[3]) +
        (3.*pow(zeta[2], 3.) - zeta[3]*pow(zeta[2], 3.))/zeta[0]/pow(1.-zeta[3], 3.);

    DispersionPolynomialState dispersion = build_dispersion_polynomials(m_avg, eta);
    const auto &a = dispersion.a;
    const auto &b = dispersion.b;
    double detI1_det = dispersion.dEtaI1_deta;
    double detI2_det = dispersion.dEtaI2_deta;
    double I1 = dispersion.I1;
    double I2 = dispersion.I2;
    double C1 = dispersion.C1;
    double C2 = dispersion.C2;

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i]*(cppargs.m[i]-1)*log(ghs[i*ncomp+i]);
    }

    double ares_hc = m_avg*ares_hs - summ;
    double ares_disp = -2*PI*den*I1*m2es3 - PI*den*m_avg*C1*I2*m2e2s3;

    summ = 0.0;
    for (int i = 0; i < ncomp; i++) {
        summ += x[i]*(cppargs.m[i]-1)/ghs[i*ncomp+i]*denghs[i*ncomp+i];
    }

    double Zhc = m_avg*Zhs - summ;
    double Zdisp = -2*PI*den*detI1_det*m2es3 - PI*den*m_avg*(C1*detI2_det + C2*eta*I2)*m2e2s3;

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
        dahc_dx = compute_contribution_dadx_fd(AresContributionKind::HC, t, rho, x, cppargs, ares_hc);
    }
    if (cppargs.disp_dadx_diff_mode == 1) {
        dadisp_dx = compute_contribution_dadx_fd(AresContributionKind::DISP, t, rho, x, cppargs, ares_disp);
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
    double Zpolar = 0.0;
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
                vector<double> adip = build_dipole_coefficients(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = build_dipole_coefficients(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
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
                    vector<double> cdip = build_dipole_coefficients(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
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
            Zpolar = eta*((dA2_det*(1-A3/A2)+(dA3_det*A2-A3*dA2_det)/A2)/(1-A3/A2)/(1-A3/A2));
            for (int i = 0; i < ncomp; i++) {
                dapolar_dx[i] = (dA2_dx[i]*(1-A3/A2) + (dA3_dx[i]*A2 - A3*dA2_dx[i])/A2)/pow(1-A3/A2,2);
            }
            if (cppargs.polar_dadx_diff_mode == 1) {
                dapolar_dx = compute_contribution_dadx_fd(AresContributionKind::POLAR, t, rho, x, cppargs, ares_polar);
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
    double Zassoc = 0.0;
    vector<double> mu_assoc(ncomp, 0);
    vector<double> daassoc_dx(ncomp, 0.0);
    double ares_assoc = 0.0;
    double sum_x_daassoc_dx = 0.0;
    if (!cppargs.e_assoc.empty()) {
        AssociationSetup assoc = build_association_setup(x, cppargs, s_ij, ghs, t);
        const vector<int> &iA = assoc.site_component_index;
        const vector<double> &x_assoc = assoc.x_assoc;
        const vector<double> &delta_ij = assoc.delta_ij;
        int num_sites = static_cast<int>(iA.size());

        vector<double> XA(num_sites, 0.0);
        for (int i = 0; i < num_sites; i++) {
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        vector<double> ddelta_dx(num_sites * num_sites * ncomp, 0);
        int idx_ddelta = 0;
        for (int k = 0; k < ncomp; k++) {
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

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-15)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < num_sites; i++) {
                dif += std::abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        vector<double> dXA_dx(num_sites*ncomp, 0);
        dXA_dx = dXAdx_find(cppargs.assoc_num, delta_ij, den, XA, ddelta_dx, x_assoc);

        int ij = 0;
        double assoc_summ = 0.0;
        for (int i = 0; i < ncomp; i++) {
            for (int j = 0; j < num_sites; j++) {
                daassoc_dx[i] += x[iA[j]]*den*dXA_dx[ij]*(1/XA[j]-0.5);
                assoc_summ += x[i]*x[iA[j]]*den*dXA_dx[ij]*(1/XA[j]-0.5);
                ij += 1;
            }
        }

        for (int i = 0; i < num_sites; i++) {
            daassoc_dx[iA[i]] += log(XA[i]) - 0.5*XA[i] + 0.5;
            ares_assoc += x[iA[i]]*(log(XA[i]) - 0.5*XA[i] + 0.5);
        }
        Zassoc = assoc_summ;
        if (cppargs.assoc_dadx_diff_mode == 1) {
            daassoc_dx = compute_contribution_dadx_fd(AresContributionKind::ASSOC, t, rho, x, cppargs, ares_assoc);
        }
        for (int i = 0; i < ncomp; i++) {
            sum_x_daassoc_dx += x[i]*daassoc_dx[i];
        }
        for (int i = 0; i < ncomp; i++) {
            mu_assoc[i] = ares_assoc + Zassoc + daassoc_dx[i] - sum_x_daassoc_dx;
        }
    }

    // Ion terms --------------------------------------------------------------
    double Zion = 0.0;
    double Zborn = 0.0;
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
        DielcState dielc_state = evaluate_dielc_state(x, cppargs);
        double eps = dielc_state.eps; // mixed dielectric constant (relative)
        vector<double> deps_dx = dielc_state.deps_dx; // d(eps_r)/dx_i
        double eps_born = eps;
        vector<double> deps_dx_born = deps_dx;
        if ((cppargs.born_model >= 1) && (cppargs.born_eps_mode == 1)) {
            eps_born = compute_eps_solvent_reference(x, cppargs);
            deps_dx_born = compute_deps_solvent_reference(x, cppargs);
        }

        double kappa = compute_dh_kappa(den, t, eps, Qsum); // inverse Debye screening length
        if ((kappa != 0) && (Qsum != 0)) {
            vector<double> chi(ncomp, 0.0);
            vector<double> sigma_k(ncomp, 0.0);
            double S = 0.;
            double Tsum = 0.;

            for (int i = 0; i < ncomp; i++) {
                double ka = kappa*d[i];
                chi[i] = compute_dh_chi(kappa, d[i]);
                sigma_k[i] = -2*chi[i] + 3/(1+ka);

                S += x[i]*cppargs.z[i]*cppargs.z[i]*chi[i];
                Tsum += x[i]*cppargs.z[i]*cppargs.z[i]*sigma_k[i];
            }

            double K0 = E_CHRG*E_CHRG/(12.0*PI*kb*t*perm_vac); // without epsilon
            double a_DH = -K0*kappa/(eps)*S;
            double Z_DH = -(K0/2.0)*kappa/(eps)*Tsum;
            Zion = Z_DH;
            a_ion = a_DH;

            vector<double> dadx(ncomp, 0.0);
            vector<double> dkappa_dx(ncomp, 0.0);
            vector<double> dS_dx(ncomp, 0.0);
            const bool use_dh_deps = (cppargs.mu_DH_comp_dep_rel_perm != 0);
            const double dh_deps_multiplier = (cppargs.mu_DH_include_sum_term != 0) ? Qsum : 1.0;
            if (cppargs.mu_DH_diff_mode == 1) {
                dadx = compute_dh_dadx_fd(t, rho, x, cppargs, a_DH);
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
            if (cppargs.debug) {
                std::cout << std::fixed << std::setprecision(10)
                          << (cppargs.mu_DH_diff_mode == 1 ? "[DEBUG DH_fd]" : "[DEBUG DH_unified]")
                          << " model=" << dh_model
                          << " eps=" << eps
                          << " kappa=" << kappa
                          << " Qsum=" << Qsum
                          << " S=" << S
                          << " Tsum=" << Tsum
                          << " use_deps=" << use_dh_deps
                          << " deps_multiplier=" << dh_deps_multiplier
                          << " sum_x_dadx=" << sum_x_dadx*8.314*t/1000.0
                          << std::endl;
            }
        }

        if (cppargs.born_model == 1) {
            // Born term (Bulow 2021a, non-SSM+DS), using user-selected born_radius_model for D_born
            double born_sum = 0.;
            for (int i = 0; i < ncomp; i++) {
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = compute_ion_born_radius(i, t, cppargs);
                    born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
                }
            }
            double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
            double a_born = -Kborn*(1.0 - 1.0/eps_born)*born_sum;

            Zborn = 0.0;

            vector<double> dadx_born(ncomp, 0.0);
            vector<double> ion_part_vec(ncomp, 0.0);
            vector<double> eps_part_vec(ncomp, 0.0);
            if (cppargs.born_diff_mode == 1) {
                dadx_born = compute_born_dadx_fd(t, x, cppargs, a_born);
            }
            else {
                for (int i = 0; i < ncomp; i++) {
                    double ion_part = 0.0;
                    if (is_ion_species(cppargs, i)) {
                        double d_born_i = compute_ion_born_radius(i, t, cppargs);
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
            if (cppargs.debug) {
                if (cppargs.born_diff_mode == 1) {
                    std::cout << std::fixed << std::setprecision(10)
                              << "[DEBUG born_model1_fd] eps=" << eps_born
                              << " born_sum=" << born_sum
                              << " a_born=" << a_born*8.314*t/1000.0
                              << " sum_x_dadx=" << sum_x_dadx_born*8.314*t/1000.0
                              << std::endl;
                    for (int i = 0; i < ncomp; i++) {
                        std::cout << "  k=" << i
                                  << " dadx_born=" << dadx_born[i]*8.314*t/1000.0
                                  << " mu_born=" << mu_born[i]*8.314*t/1000.0
                                  << std::endl;
                    }
                }
                else {
                std::cout << std::fixed << std::setprecision(10)
                          << "[DEBUG born_model1_m" << cppargs.born_diff_mode << "] eps=" << eps_born
                          << " born_sum=" << born_sum
                          << " a_born=" << a_born*8.314*t/1000.0
                          << " sum_x_dadx=" << sum_x_dadx_born*8.314*t/1000.0
                          << std::endl;
                for (int i = 0; i < ncomp; i++) {
                    std::cout << "  k=" << i
                              << " ion_part=" << ion_part_vec[i]
                              << " eps_part=" << eps_part_vec[i]
                              << " dadx_born=" << dadx_born[i]*8.314*t/1000.0
                              << " mu_born=" << mu_born[i]*8.314*t/1000.0
                              << std::endl;
                }
                }
            }
        }
        else if (cppargs.born_model == 2) {
            const double eps_r_ion = 8.0;
            const double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
            BornSSMDSData born = build_born_ssmds_data(x, cppargs, t, eps_born, eps_r_ion);
            double a_born = -Kborn*born.sum_bracket;
            Zborn = 0.0;

            vector<double> dadx_born(ncomp, 0.0);
            vector<double> direct_part_vec(ncomp, 0.0);
            vector<double> deps_part_vec(ncomp, 0.0);
            vector<double> ddelta_part_vec(ncomp, 0.0);
            if (cppargs.born_diff_mode == 1) {
                dadx_born = compute_born_dadx_fd(t, x, cppargs, a_born);
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
            if (cppargs.debug) {
                double f_mix_dbg = 0.0;
                for (int i = 0; i < ncomp; i++) {
                    f_mix_dbg += x[i]*born.f_k[i];
                }
                if (cppargs.born_diff_mode == 1) {
                    std::cout << std::fixed << std::setprecision(10)
                              << "[DEBUG born_model" << cppargs.born_model << "_fd] eps=" << eps_born
                              << " eps_ion=" << eps_r_ion
                              << " f_mix=" << f_mix_dbg
                              << " sum_bracket=" << born.sum_bracket
                              << " a_born=" << a_born*8.314*t/1000.0
                              << " sum_x_dadx=" << sum_x_dadx_born*8.314*t/1000.0
                              << std::endl;
                    for (int i = 0; i < ncomp; i++) {
                        std::cout << "  k=" << i
                                  << " z=" << cppargs.z[i]
                                  << " d_born=" << born.d_born[i]
                                  << " D=" << born.D[i]
                                  << " f_k=" << born.f_k[i]
                                  << " bracket=" << born.bracket[i]
                                  << " dadx_born=" << dadx_born[i]*8.314*t/1000.0
                                  << " mu_born=" << mu_born[i]*8.314*t/1000.0
                                  << std::endl;
                    }
                }
                else {
                    std::cout << std::fixed << std::setprecision(10)
                              << "[DEBUG born_model" << cppargs.born_model << "] eps=" << eps_born
                              << " eps_ion=" << eps_r_ion
                              << " f_mix=" << f_mix_dbg
                              << " sum_bracket=" << born.sum_bracket
                              << " sum_invD=" << born.sum_invD
                              << " sum_gap=" << born.sum_gap
                              << " sum_dpref_over_D2=" << born.sum_dpref_over_D2
                              << " a_born=" << a_born*8.314*t/1000.0
                              << " sum_x_dadx=" << sum_x_dadx_born*8.314*t/1000.0
                              << std::endl;
                    for (int i = 0; i < ncomp; i++) {
                        std::cout << "  k=" << i
                                  << " z=" << cppargs.z[i]
                                  << " d_born=" << born.d_born[i]
                                  << " D=" << born.D[i]
                                  << " f_k=" << born.f_k[i]
                                  << " bracket=" << born.bracket[i]
                                  << " direct=" << direct_part_vec[i]
                                  << " deps=" << deps_part_vec[i]
                                  << " ddelta=" << ddelta_part_vec[i]
                                  << " dadx_born=" << dadx_born[i]*8.314*t/1000.0
                                  << " mu_born=" << mu_born[i]*8.314*t/1000.0
                                  << std::endl;
                    }
                }
            }
        }
        else if (cppargs.born_model != 0) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }
    double Z = pcsaft_Z_cpp(t, rho, x, cppargs);

    vector<double> mu(ncomp, 0);
    vector<double> lnfugcoef(ncomp, 0);
    vector<double> raw_z_terms = {Zhc, Zdisp, Zpolar, Zassoc, Zion, Zborn};
    vector<double> norm_z_terms = normalize_z_contributions(raw_z_terms, Z);
    double z_weight = stable_logz_over_zminus1(Z);
    vector<double> lnfug_hc(ncomp, 0.0);
    vector<double> lnfug_disp(ncomp, 0.0);
    vector<double> lnfug_polar(ncomp, 0.0);
    vector<double> lnfug_assoc(ncomp, 0.0);
    vector<double> lnfug_ion(ncomp, 0.0);
    vector<double> lnfug_born(ncomp, 0.0);
    for (int i = 0; i < ncomp; i++) {
        mu[i] = mu_hc[i] + mu_disp[i] + mu_polar[i] + mu_assoc[i] + mu_ion[i] + mu_born[i];
        lnfug_hc[i] = mu_hc[i] - norm_z_terms[0] * z_weight;
        lnfug_disp[i] = mu_disp[i] - norm_z_terms[1] * z_weight;
        lnfug_polar[i] = mu_polar[i] - norm_z_terms[2] * z_weight;
        lnfug_assoc[i] = mu_assoc[i] - norm_z_terms[3] * z_weight;
        lnfug_ion[i] = mu_ion[i] - norm_z_terms[4] * z_weight;
        lnfug_born[i] = mu_born[i] - norm_z_terms[5] * z_weight;
        lnfugcoef[i] = lnfug_hc[i] + lnfug_disp[i] + lnfug_polar[i] + lnfug_assoc[i] + lnfug_ion[i] + lnfug_born[i];
    }

    // Cache per-term residual chemical-potential contributions for structured API access.
    g_last_mu_hc = mu_hc;
    g_last_mu_disp = mu_disp;
    g_last_mu_polar = mu_polar;
    g_last_mu_assoc = mu_assoc;
    g_last_mu_ion = mu_ion;
    g_last_mu_born = mu_born;
    g_last_mu_total = mu;
    g_last_lnfug_hc = lnfug_hc;
    g_last_lnfug_disp = lnfug_disp;
    g_last_lnfug_polar = lnfug_polar;
    g_last_lnfug_assoc = lnfug_assoc;
    g_last_lnfug_ion = lnfug_ion;
    g_last_lnfug_born = lnfug_born;
    g_last_lnfugcoef = lnfugcoef;
    g_last_dadx_hc = dahc_dx;
    g_last_dadx_disp = dadisp_dx;
    g_last_dadx_polar = dapolar_dx;
    g_last_dadx_assoc = daassoc_dx;
    g_last_dadx_ion = dadx_ion;
    g_last_dadx_born = dadx_born_diag;
    g_last_a_hc = ares_hc;
    g_last_a_disp = ares_disp;
    g_last_a_polar = ares_polar;
    g_last_a_assoc = ares_assoc;
    g_last_a_ion = a_ion;
    g_last_a_born = a_born_diag;
    g_last_sum_x_dadx_hc = sum_x_dahc_dx;
    g_last_sum_x_dadx_disp = sum_x_dadisp_dx;
    g_last_sum_x_dadx_polar = sum_x_dapolar_dx;
    g_last_sum_x_dadx_assoc = sum_x_daassoc_dx;
    g_last_sum_x_dadx_ion = sum_x_dadx_ion;
    g_last_sum_x_dadx_born = sum_x_dadx_born_diag;
    g_last_z_raw_hc = Zhc;
    g_last_z_raw_disp = Zdisp;
    g_last_z_raw_polar = Zpolar;
    g_last_z_raw_assoc = Zassoc;
    g_last_z_raw_ion = Zion;
    g_last_z_raw_born = Zborn;
    g_last_z_hc = norm_z_terms[0];
    g_last_z_disp = norm_z_terms[1];
    g_last_z_polar = norm_z_terms[2];
    g_last_z_assoc = norm_z_terms[3];
    g_last_z_ion = norm_z_terms[4];
    g_last_z_born = norm_z_terms[5];
    g_last_z_total = Z;
    if (cppargs.debug) {
        std::cout << std::fixed << std::setprecision(10)
                  << "[DEBUG mu_res] t=" << t << " rho=" << rho << std::endl;
        for (int i = 0; i < ncomp; i++) {
            std::cout << "  i=" << i
                      << " mu_hc=" << mu_hc[i]*8.314*t/1000
                      << " mu_disp=" << mu_disp[i]*8.314*t/1000
                      << " mu_polar=" << mu_polar[i]*8.314*t/1000
                      << " mu_assoc=" << mu_assoc[i]*8.314*t/1000
                      << " mu_ion=" << mu_ion[i]*8.314*t/1000
                      << " mu_born=" << mu_born[i]*8.314*t/1000
                      << " mu_res=" << mu[i]*8.314*t/1000
                      << " lnphi=" << lnfugcoef[i]*8.314*t/1000 << std::endl;
        }
    }

    return lnfugcoef;
}


vector<double> pcsaft_lnfug_terms_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate per-term residual chemical-potential contributions and ln fugacity coefficients
    for one phase of the system.

    Output layout (flattened blocks, each of length ncomp):
      [mu_hc, mu_disp, mu_polar, mu_assoc, mu_ion, mu_born, mu_total, lnfugcoef,
       lnfug_hc, lnfug_disp, lnfug_polar, lnfug_assoc, lnfug_ion, lnfug_born,
       dadx_hc, dadx_disp, dadx_polar, dadx_assoc, dadx_ion, dadx_born,
       a_hc, a_disp, a_polar, a_assoc, a_ion, a_born,
       sum_x_dadx_hc, sum_x_dadx_disp, sum_x_dadx_polar, sum_x_dadx_assoc, sum_x_dadx_ion, sum_x_dadx_born,
       Zraw_hc, Zraw_disp, Zraw_polar, Zraw_assoc, Zraw_ion, Zraw_born,
       Z_hc, Z_disp, Z_polar, Z_assoc, Z_ion, Z_born, Z_total]
    */
    vector<double> lnfug = pcsaft_lnfug_cpp(t, rho, x, cppargs);
    int ncomp = static_cast<int>(x.size());

    if ((static_cast<int>(g_last_mu_hc.size()) != ncomp) ||
        (static_cast<int>(g_last_mu_disp.size()) != ncomp) ||
        (static_cast<int>(g_last_mu_polar.size()) != ncomp) ||
        (static_cast<int>(g_last_mu_assoc.size()) != ncomp) ||
        (static_cast<int>(g_last_mu_ion.size()) != ncomp) ||
        (static_cast<int>(g_last_mu_born.size()) != ncomp) ||
        (static_cast<int>(g_last_mu_total.size()) != ncomp) ||
        (static_cast<int>(g_last_lnfug_hc.size()) != ncomp) ||
        (static_cast<int>(g_last_lnfug_disp.size()) != ncomp) ||
        (static_cast<int>(g_last_lnfug_polar.size()) != ncomp) ||
        (static_cast<int>(g_last_lnfug_assoc.size()) != ncomp) ||
        (static_cast<int>(g_last_lnfug_ion.size()) != ncomp) ||
        (static_cast<int>(g_last_lnfug_born.size()) != ncomp) ||
        (static_cast<int>(g_last_dadx_hc.size()) != ncomp) ||
        (static_cast<int>(g_last_dadx_disp.size()) != ncomp) ||
        (static_cast<int>(g_last_dadx_polar.size()) != ncomp) ||
        (static_cast<int>(g_last_dadx_assoc.size()) != ncomp) ||
        (static_cast<int>(g_last_dadx_ion.size()) != ncomp) ||
        (static_cast<int>(g_last_dadx_born.size()) != ncomp) ||
        (static_cast<int>(g_last_lnfugcoef.size()) != ncomp) ||
        (static_cast<int>(lnfug.size()) != ncomp)) {
        throw ValueError("Internal lnfug term cache size mismatch.");
    }

    vector<double> out(20 * ncomp + 25, 0.0);
    auto copy_block = [&](int block, const vector<double> &src) {
        for (int i = 0; i < ncomp; i++) {
            out[block*ncomp + i] = src[i];
        }
    };

    copy_block(0, g_last_mu_hc);
    copy_block(1, g_last_mu_disp);
    copy_block(2, g_last_mu_polar);
    copy_block(3, g_last_mu_assoc);
    copy_block(4, g_last_mu_ion);
    copy_block(5, g_last_mu_born);
    copy_block(6, g_last_mu_total);
    copy_block(7, g_last_lnfugcoef);
    copy_block(8, g_last_lnfug_hc);
    copy_block(9, g_last_lnfug_disp);
    copy_block(10, g_last_lnfug_polar);
    copy_block(11, g_last_lnfug_assoc);
    copy_block(12, g_last_lnfug_ion);
    copy_block(13, g_last_lnfug_born);
    copy_block(14, g_last_dadx_hc);
    copy_block(15, g_last_dadx_disp);
    copy_block(16, g_last_dadx_polar);
    copy_block(17, g_last_dadx_assoc);
    copy_block(18, g_last_dadx_ion);
    copy_block(19, g_last_dadx_born);
    int scalar_offset = 20 * ncomp;
    out[scalar_offset + 0] = g_last_a_hc;
    out[scalar_offset + 1] = g_last_a_disp;
    out[scalar_offset + 2] = g_last_a_polar;
    out[scalar_offset + 3] = g_last_a_assoc;
    out[scalar_offset + 4] = g_last_a_ion;
    out[scalar_offset + 5] = g_last_a_born;
    out[scalar_offset + 6] = g_last_sum_x_dadx_hc;
    out[scalar_offset + 7] = g_last_sum_x_dadx_disp;
    out[scalar_offset + 8] = g_last_sum_x_dadx_polar;
    out[scalar_offset + 9] = g_last_sum_x_dadx_assoc;
    out[scalar_offset + 10] = g_last_sum_x_dadx_ion;
    out[scalar_offset + 11] = g_last_sum_x_dadx_born;
    out[scalar_offset + 12] = g_last_z_raw_hc;
    out[scalar_offset + 13] = g_last_z_raw_disp;
    out[scalar_offset + 14] = g_last_z_raw_polar;
    out[scalar_offset + 15] = g_last_z_raw_assoc;
    out[scalar_offset + 16] = g_last_z_raw_ion;
    out[scalar_offset + 17] = g_last_z_raw_born;
    out[scalar_offset + 18] = g_last_z_hc;
    out[scalar_offset + 19] = g_last_z_disp;
    out[scalar_offset + 20] = g_last_z_polar;
    out[scalar_offset + 21] = g_last_z_assoc;
    out[scalar_offset + 22] = g_last_z_ion;
    out[scalar_offset + 23] = g_last_z_born;
    out[scalar_offset + 24] = g_last_z_total;
    return out;
}


vector<double> pcsaft_fugcoef_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate the fugacity coefficients for one phase of the system.
    */
    int ncomp = static_cast<int>(x.size()); // number of components
    vector<double> lnfug = pcsaft_lnfug_cpp(t, rho, x, cppargs);
    vector<double> fugcoef(ncomp, 0);
    for (int i = 0; i < ncomp; i++) {
        fugcoef[i] = exp(lnfug[i]); // the fugacity coefficients
    }

    return fugcoef;
}


double pcsaft_p_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate pressure
    */
    double den = rho*N_AV/1.0e30;

    double Z = pcsaft_Z_cpp(t, rho, x, cppargs);
    double P = Z*kb*t*den*1.0e30; // Pa
    return P;
}


double pcsaft_ares_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate the residual Helmholtz energy
    */
    int ncomp = static_cast<int>(x.size()); // number of components
    double summ = 0.0;
    ThermoCommonState thermo = build_thermo_common_state(t, rho, x, cppargs, false);
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
    DispersionPolynomialState dispersion = build_dispersion_polynomials(m_avg, eta);
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
                vector<double> adip = build_dipole_coefficients(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = build_dipole_coefficients(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
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
                    vector<double> cdip = build_dipole_coefficients(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
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
        AssociationSetup assoc = build_association_setup(x, cppargs, s_ij, ghs, t);
        const vector<int> &iA = assoc.site_component_index;
        const vector<double> &x_assoc = assoc.x_assoc;
        const vector<double> &delta_ij = assoc.delta_ij;
        int num_sites = static_cast<int>(iA.size());

        vector<double> XA(num_sites, 0.0);
        for (int i = 0; i < num_sites; i++) {
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-15)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < num_sites; i++) {
                dif += std::abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        ares_assoc = 0.;
        for (int i = 0; i < num_sites; i++) {
            ares_assoc += x[iA[i]]*(log(XA[i])-0.5*XA[i] + 0.5);
        }
    }

    // Ion term ---------------------------------------------------------------
    double ares_ion = 0.;
    double ares_born = 0.;
    if (!cppargs.z.empty()) {
        DielcState dielc_state = evaluate_dielc_state(x, cppargs);
        double eps = dielc_state.eps;
        double eps_born = (cppargs.born_eps_mode == 1) ? compute_eps_solvent_reference(x, cppargs) : eps;
        ares_ion = compute_dh_ares_only(t, rho, x, cppargs);

        if (cppargs.born_model == 1) {
            // Born term (Bulow 2021a, non-SSM+DS): use d_born,i in denominator
            double born_sum = 0.;
            for (int i = 0; i < ncomp; i++) {
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = compute_ion_born_radius(i, t, cppargs);
                    born_sum += x[i]*cppargs.z[i]*cppargs.z[i]/d_born_i;
                }
            }
            ares_born = -E_CHRG*E_CHRG/(4.*PI*kb*t*perm_vac)*(1.-1./eps_born)*born_sum;
        }
        else if (cppargs.born_model == 2) {
            const double eps_r_ion = 8.0;
            const double Kborn = E_CHRG*E_CHRG/(4.0*PI*kb*t*perm_vac);
            BornSSMDSData born = build_born_ssmds_data(x, cppargs, t, eps_born, eps_r_ion);
            ares_born = -Kborn*born.sum_bracket;
        }
        else if (cppargs.born_model != 0) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }

    double ares = ares_hc + ares_disp + ares_polar + ares_assoc + ares_ion + ares_born;
    if (cppargs.debug) {
        std::cout << std::fixed << std::setprecision(10)
                  << "[DEBUG ares] t=" << t << " rho=" << rho
                  << " ares_hc=" << ares_hc
                  << " ares_disp=" << ares_disp
                  << " ares_polar=" << ares_polar
                  << " ares_assoc=" << ares_assoc
                  << " ares_ion=" << ares_ion
                  << " ares_born=" << ares_born
                  << " ares_total=" << ares << std::endl;
    }
    return ares;
}


double pcsaft_dadt_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate the temperature derivative of the residual Helmholtz energy at
    constant density.
    */
    int ncomp = static_cast<int>(x.size()); // number of components
    double summ = 0.0;
    ThermoCommonState thermo = build_thermo_common_state(t, rho, x, cppargs, true);
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
    DispersionPolynomialState dispersion = build_dispersion_polynomials(m_avg, eta);
    vector<double> dghs_dt(ncomp * ncomp, 0.0);
    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            double pair_diameter = compute_pair_diameter(d[i], d[j]);
            double pair_diameter_dt = pair_diameter * (
                dd_dt[i] / d[i] + dd_dt[j] / d[j] - (dd_dt[i] + dd_dt[j]) / (d[i] + d[j])
            );
            dghs_dt[idx] = compute_hs_contact_time_derivative(
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
                vector<double> adip = build_dipole_coefficients(kDipoleA0, kDipoleA1, kDipoleA2, m_ij);
                vector<double> bdip = build_dipole_coefficients(kDipoleB0, kDipoleB1, kDipoleB2, m_ij);
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
                    vector<double> cdip = build_dipole_coefficients(kDipoleC0, kDipoleC1, kDipoleC2, m_ijk);
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
        AssociationSetup assoc = build_association_setup(x, cppargs, s_ij, ghs, t);
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
            XA[i] = (-1 + sqrt(1+8*den*delta_ij[i*num_sites+i]))/(4*den*delta_ij[i*num_sites+i]);
            if (!std::isfinite(XA[i])) {
                XA[i] = 0.02;
            }
        }

        int ctr = 0;
        double dif = 1000.;
        vector<double> XA_old = XA;
        while ((ctr < 100) && (dif > 1e-15)) {
            ctr += 1;
            XA = XA_find(XA_old, delta_ij, den, x_assoc);
            dif = 0.;
            for (int i = 0; i < num_sites; i++) {
                dif += std::abs(XA[i] - XA_old[i]);
            }
            for (int i = 0; i < num_sites; i++) {
                XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
            }
        }

        vector<double> dXA_dt(num_sites, 0);
        dXA_dt = dXAdt_find(delta_ij, den, XA, ddelta_dt, x_assoc);

        for (int i = 0; i < num_sites; i++) {
            dadt_assoc += x[iA[i]]*(1/XA[i]-0.5)*dXA_dt[i];
        }
    }

    // Ion term ---------------------------------------------------------------
    double dadt_ion = 0.;
    double dadt_born = 0.;
    if (!cppargs.z.empty()) {
        DielcState dielc_state = evaluate_dielc_state(x, cppargs);
        double eps = dielc_state.eps;
        double eps_born = (cppargs.born_eps_mode == 1) ? compute_eps_solvent_reference(x, cppargs) : eps;
        vector<double> q(cppargs.z.begin(), cppargs.z.end());
        for (int i = 0; i < ncomp; i++) {
            q[i] = q[i]*E_CHRG;
        }

        summ = 0.;
        for (int i = 0; i < ncomp; i++) {
            summ += cppargs.z[i]*cppargs.z[i]*x[i];
        }
        double kappa = compute_dh_kappa(den, t, eps, summ); // the inverse Debye screening length. Equation 4 in Held et al. 2008.

        double dkappa_dt;
        if (kappa != 0) {
            vector<double> chi(ncomp);
            vector<double> dchikap_dk(ncomp);
            summ = 0.;
            for (int i = 0; i < ncomp; i++) {
                chi[i] = compute_dh_chi(kappa, d[i]);
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
                    double d_born_i = compute_ion_born_radius(i, t, cppargs);
                    double d_born_dt = compute_ion_born_radius_dt(i, t, cppargs);
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
            BornSSMDSData born = build_born_ssmds_data(x, cppargs, t, eps_born, eps_r_ion);
            dadt_born = E_CHRG*E_CHRG/(4.*PI*kb*perm_vac*t*t)*born.sum_bracket;
        }
        else if (cppargs.born_model != 0) {
            throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
        }
    }

    double dadt = dadt_hc + dadt_disp + dadt_assoc + dadt_polar + dadt_ion + dadt_born;
    return dadt;
}


double pcsaft_hres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate the residual enthalpy for one phase of the system.
    */
    double Z = pcsaft_Z_cpp(t, rho, x, cppargs);
    double dares_dt = pcsaft_dadt_cpp(t, rho, x, cppargs);

    double hres = (-t*dares_dt + (Z-1))*kb*N_AV*t; // Equation A.46 from Gross and Sadowski 2001
    return hres;
}


double pcsaft_sres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate the residual entropy (constant volume) for one phase of the system.
    */
    double gres = pcsaft_gres_cpp(t, rho, x, cppargs);
    double hres = pcsaft_hres_cpp(t, rho, x, cppargs);

    double sres = (hres - gres)/t;
    return sres;
}

double pcsaft_gres_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    /**
    Calculate the residual Gibbs energy for one phase of the system.
    */
    double ares = pcsaft_ares_cpp(t, rho, x, cppargs);
    double Z = pcsaft_Z_cpp(t, rho, x, cppargs);

    double gres = (ares + (Z - 1) - log(Z))*kb*N_AV*t; // Equation A.50 from Gross and Sadowski 2001
    return gres;
}


vector<double> flashTQ_cpp(double t, double Q, vector<double> x, const add_args &cppargs) {
    bool solution_found = false;
    double p_guess;
    vector<double> output;
    try {
        p_guess = estimate_flash_p(t, Q, x, cppargs);
        output = outerTQ(p_guess, t, Q, x, cppargs);
        solution_found = true;
    }
    catch (const SolutionError&) {}
    catch (const ValueError&) {}

    // if solution hasn't been found, try cycling through a range of pressures
    if (!solution_found) {
        double p_lbound = -6; // here we're using log10 of the pressure
        double p_ubound = 9;
        double p_step = 0.1;
        p_guess = p_lbound;
        while (p_guess < p_ubound && !solution_found) {
            try {
                output = outerTQ(pow(10, p_guess), t, Q, x, cppargs);
                solution_found = true;
            } catch (const SolutionError&) {
                p_guess += p_step;
            } catch (const ValueError&) {
                p_guess += p_step;
            }
        }
    }

    if (!solution_found) {
        if (is_effectively_pure_feed(x)) {
            return pure_component_flashTQ(t, x, cppargs);
        }
        throw SolutionError("solution could not be found for TQ flash");
    }
    return output;
}


vector<double> flashTQ_cpp(double t, double Q, vector<double> x, const add_args &cppargs, double p_guess) {
    vector<double> output;
    try {
        output = outerTQ(p_guess, t, Q, x, cppargs);
    } catch (const SolutionError& ex) {
        if (is_effectively_pure_feed(x)) {
            return pure_component_flashTQ(t, x, cppargs);
        }
        throw ex;
    }

    return output;
}


vector<double> flashPQ_cpp(double p, double Q, vector<double> x, const add_args &cppargs){
    bool solution_found = false;
    double t_guess;
    vector<double> output;
    try {
        t_guess = estimate_flash_t(p, Q, x, cppargs);
        output = outerPQ(t_guess, p, Q, x, cppargs);
        solution_found = true;
    }
    catch (const SolutionError&) {}
    catch (const ValueError&) {}

    // if solution hasn't been found, try calling the flash function directly with a range of initial temperatures
    if (!solution_found) {
        double t_lbound = 1;
        double t_ubound = 800;
        double t_step = 10;
        if (!cppargs.z.empty()) {
            t_lbound = 264;
            t_ubound = 350;
        }
        t_guess = t_ubound;
        while (t_guess > t_lbound && !solution_found) {
            try {
                output = outerPQ(t_guess, p, Q, x, cppargs);
                solution_found = true;
            } catch (const SolutionError&) {
                t_guess -= t_step;
            } catch (const ValueError&) {
                t_guess -= t_step;
            }
        }
    }

    if (!solution_found) {
        throw SolutionError("solution could not be found for PQ flash");
    }
    return output;
}


vector<double> flashPQ_cpp(double p, double Q, vector<double> x, const add_args &cppargs, double t_guess){
    vector<double> output;
    try {
        output = outerPQ(t_guess, p, Q, x, cppargs);
    } catch (const SolutionError&) {
        output = flashPQ_cpp(p, Q, x, cppargs); // call function without an initial guess
    }

    return output;
}


vector<double> outerPQ(double t_guess, double p, double Q, vector<double> x, add_args cppargs) {
    // Based on the algorithm proposed in H. A. J. Watson, M. Vikse, T. Gundersen, and P. I. Barton, “Reliable Flash Calculations: Part 1. Nonsmooth Inside-Out Algorithms,” Ind. Eng. Chem. Res., vol. 56, no. 4, pp. 960–973, Feb. 2017, doi: 10.1021/acs.iecr.6b03956.
    int ncomp = static_cast<int>(x.size());
    double TOL = 1e-8;
    double MAXITER = 200;

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (!cppargs.z.empty() && cppargs.z[i] != 0) {
            x_ions += x[i];
        }
    }

    // initialize variables
    vector<double> fugcoef_l(ncomp), fugcoef_v(ncomp), k(ncomp), u(ncomp), kprime(ncomp);
    double rhol, rhov;
    double Tref = t_guess - 1;
    double Tprime = t_guess + 1;
    double t = t_guess;

    // calculate sigma for water, if it is present
    std::vector<double>::iterator water_iter = std::find(cppargs.e.begin(), cppargs.e.end(), 353.9449);
    int water_idx = -1;
    if (water_iter != cppargs.e.end()) {
        water_idx = static_cast<int>(std::distance(cppargs.e.begin(), water_iter));
        cppargs.s[water_idx] = calc_sigma(t, &calc_water_sigma);
    }

    // calculate initial guess for compositions based on fugacity coefficients and Raoult's Law.
    rhol = pcsaft_den_cpp(t, p, x, 0, cppargs);
    rhov = pcsaft_den_cpp(t, p, x, 1, cppargs);
    if ((rhol - rhov) < 1e-4) {
        throw SolutionError("liquid and vapor densities are the same.");
    }
    fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, x, cppargs);
    fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, x, cppargs);

    for (int i = 0; i < ncomp; i++) {
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            k[i] = fugcoef_l[i] / fugcoef_v[i];
        } else {
            k[i] = 0; // set k to 0 for ionic components
        }
    }

    vector<double> xl(ncomp);
    vector<double> xv(ncomp);
    double xv_sum = 0;
    double xl_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        xl[i] = x[i] / (1 + Q * (k[i] - 1));
        xl_sum += xl[i];
        xv[i] = k[i] * x[i] / (1 + Q * (k[i] - 1));
        xv_sum += xv[i];
    }

    if (xv_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            xv[i] = xv[i] / xv_sum;
        }
    }

    if (xl_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            xl[i] = xl[i] / xl_sum;
        }
    }

    rhol = pcsaft_den_cpp(t, p, xl, 0, cppargs);
    fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, xl, cppargs);
    rhov = pcsaft_den_cpp(t, p, xv, 1, cppargs);
    fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, xv, cppargs);
    for (int i = 0; i < ncomp; i++) {
        k[i] = fugcoef_l[i] / fugcoef_v[i];
        u[i] = std::log(k[i] / kb);
    }

    if (water_idx >= 0) {
        cppargs.s[water_idx] = calc_sigma(Tprime, &calc_water_sigma);
    }
    rhol = pcsaft_den_cpp(Tprime, p, xl, 0, cppargs);
    fugcoef_l = pcsaft_fugcoef_cpp(Tprime, rhol, xl, cppargs);
    rhov = pcsaft_den_cpp(Tprime, p, xv, 1, cppargs);
    fugcoef_v = pcsaft_fugcoef_cpp(Tprime, rhov, xv, cppargs);
    for (int i = 0; i < ncomp; i++) {
        kprime[i] = fugcoef_l[i] / fugcoef_v[i];
    }

    vector<double> t_weight(ncomp);
    double t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Tprime - t);
        t_weight[i] = xv[i] * dlnk_dt / (1 + Q * (k[i] - 1));
        t_sum += t_weight[i];
    }

    double kb = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            kb += wi * std::log(k[i]);
        }
    }
    kb = std::exp(kb);

    t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Tprime - t);
        t_weight[i] = xv[i] * dlnk_dt / (1 + Q * (kprime[i] - 1));
        t_sum += t_weight[i];
    }

    double kbprime = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            kbprime += wi * std::log(kprime[i]);
        }
    }
    kbprime = std::exp(kbprime);
    double kb0 = kbprime;

    for (int i = 0; i < ncomp; i++) {
        u[i] = std::log(k[i] / kb);
    }

    double B = std::log(kbprime / kb) / (1/Tprime - 1/t);
    double A = std::log(kb) - B * (1/t - 1/Tref);
    if (B > 0) {
        throw SolutionError("B > 0 in outerPQ");
    }

    // solve
    vector<double> pp(ncomp);
    double maxdif = 1e10 * TOL;
    int itr = 0;
    double Rmin = 0, Rmax = 1;
    while (maxdif > TOL && itr < MAXITER) {
        // save previous values for calculating the difference at the end of the iteration
        vector<double> u_old = u;
        double A_old = A;

        double R0 = kb * Q / (kb * Q + kb0 * (1 - Q));
        double R = BoundedSecantInner(kb0, Q, u, x, cppargs, R0, Rmin, Rmax, DBL_EPSILON, 1e-8, 200);

        double pp_sum = 0;
        double eupp_sum = 0;
        for (int i = 0; i < ncomp; i++) {
            pp[i] = x[i] / (1 - R + kb0 * R * std::exp(u[i]));
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                pp_sum += pp[i];
                eupp_sum += std::exp(u[i]) * pp[i];
            }
        }
        kb = pp_sum / eupp_sum;

        t = 1 / (1 / Tref + (std::log(kb) - A) / B);
        for (int i = 0; i < ncomp; i++) {
            if (x_ions == 0) {
                xl[i] = pp[i] / pp_sum;
                xv[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else if (cppargs.z.empty() || cppargs.z[i] == 0) {
                xl[i] = pp[i] / pp_sum * (1 - x_ions/(1-Q));
                xv[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else {
                xl[i] = x[i] / (1 - Q);
                xv[i] = 0;
            }
        }

        if (water_idx >= 0) {
            cppargs.s[water_idx] = calc_sigma(t, &calc_water_sigma);
        }
        rhol = pcsaft_den_cpp(t, p, xl, 0, cppargs);
        fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, xl, cppargs);
        rhov = pcsaft_den_cpp(t, p, xv, 1, cppargs);
        fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, xv, cppargs);
        for (int i = 0; i < ncomp; i++) {
            k[i] = fugcoef_l[i] / fugcoef_v[i];
            u[i] = std::log(k[i] / kb);
        }

        if (itr == 0) {
            B = std::log(kbprime / kb) / (1/Tprime - 1/t);
            if (B > 0) {
                throw SolutionError("B > 0 in outerPQ");
            }
        }
        A = std::log(kb) - B * (1/t - 1/Tref);

        maxdif = std::abs(A - A_old);
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                double dif = std::abs(u[i] - u_old[i]);
                if (dif > maxdif) {
                    maxdif = dif;
                }
            }
        }

        itr += 1;
    }

    if (!std::isfinite(t) || maxdif > 1e-3 || t < 0) {
        throw SolutionError("outerPQ did not converge to a solution");
    }

    vector<double> result;
    result.push_back(t);
    result.insert(result.end(), xl.begin(), xl.end());
    result.insert(result.end(), xv.begin(), xv.end());
    return result;
}

vector<double> outerTQ(double p_guess, double t, double Q, vector<double> x, add_args cppargs) {
    // Based on the algorithm proposed in H. A. J. Watson, M. Vikse, T. Gundersen, and P. I. Barton, “Reliable Flash Calculations: Part 1. Nonsmooth Inside-Out Algorithms,” Ind. Eng. Chem. Res., vol. 56, no. 4, pp. 960–973, Feb. 2017, doi: 10.1021/acs.iecr.6b03956.
    int ncomp = static_cast<int>(x.size());
    double TOL = 1e-8;
    double MAXITER = 200;

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (!cppargs.z.empty() && cppargs.z[i] != 0) {
            x_ions += x[i];
        }
    }

    // initialize variables
    vector<double> fugcoef_l(ncomp), fugcoef_v(ncomp), k(ncomp), u(ncomp), kprime(ncomp);
    double rhol, rhov;
    double Pref = p_guess - 0.01 * p_guess;
    double Pprime = p_guess + 0.01 * p_guess;
    double p = p_guess;

    // calculate initial guess for compositions based on fugacity coefficients and Raoult's Law.
    rhol = pcsaft_den_cpp(t, p, x, 0, cppargs);
    rhov = pcsaft_den_cpp(t, p, x, 1, cppargs);
    if ((rhol - rhov) < 1e-4) {
        throw SolutionError("liquid and vapor densities are the same.");
    }
    fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, x, cppargs);
    fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, x, cppargs);

    for (int i = 0; i < ncomp; i++) {
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            k[i] = fugcoef_l[i] / fugcoef_v[i];
        } else {
            k[i] = 0; // set k to 0 for ionic components
        }
    }

    vector<double> xl(ncomp);
    vector<double> xv(ncomp);
    double xv_sum = 0;
    double xl_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        xl[i] = x[i] / (1 + Q * (k[i] - 1));
        xl_sum += xl[i];
        xv[i] = k[i] * x[i] / (1 + Q * (k[i] - 1));
        xv_sum += xv[i];
    }

    if (xv_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            xv[i] = xv[i] / xv_sum;
        }
    }

    if (xl_sum != 1) {
        for (int i = 0; i < ncomp; i++) {
            xl[i] = xl[i] / xl_sum;
        }
    }

    FlashPhaseState base_state;
    if (!try_flash_phase_state(t, p, xl, xv, cppargs, &base_state)) {
        throw SolutionError("No valid density roots found for outerTQ base state.");
    }
    rhol = base_state.rhol;
    fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, xl, cppargs);
    rhov = base_state.rhov;
    fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, xv, cppargs);
    for (int i = 0; i < ncomp; i++) {
        k[i] = fugcoef_l[i] / fugcoef_v[i];
        u[i] = std::log(k[i] / kb);
    }

    FlashPhaseState prime_state = find_nearby_flash_phase_state(t, p, xl, xv, cppargs, p_guess > 1e6 ? -1 : 1);
    if (!prime_state.valid) {
        throw SolutionError("No valid nearby pressure state found for outerTQ.");
    }
    Pprime = prime_state.p;
    double delta_p = Pprime - p;
    Pref = p - delta_p;
    if (!(Pref > 0.0)) {
        Pref = std::max(0.5 * p, p * (1.0 - 1e-4));
    }

    rhol = prime_state.rhol;
    fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, xl, cppargs);
    rhov = prime_state.rhov;
    fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, xv, cppargs);
    for (int i = 0; i < ncomp; i++) {
        kprime[i] = fugcoef_l[i] / fugcoef_v[i];
    }

    vector<double> t_weight(ncomp);
    double t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Pprime - p);
        t_weight[i] = xv[i] * dlnk_dt / (1 + Q * (k[i] - 1));
        t_sum += t_weight[i];
    }

    double kb = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            kb += wi * std::log(k[i]);
        }
    }
    kb = std::exp(kb);

    t_sum = 0;
    for (int i = 0; i < ncomp; i++) {
        double dlnk_dt = (kprime[i] - k[i]) / (Pprime - p);
        t_weight[i] = xv[i] * dlnk_dt / (1 + Q * (kprime[i] - 1));
        t_sum += t_weight[i];
    }

    double kbprime = 0;
    for (int i = 0; i < ncomp; i++) {
        double wi = t_weight[i] / t_sum;
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            kbprime += wi * std::log(kprime[i]);
        }
    }
    kbprime = std::exp(kbprime);
    double kb0 = kbprime;

    for (int i = 0; i < ncomp; i++) {
        u[i] = std::log(k[i] / kb);
    }

    double B = std::log(kbprime / kb) / (1/Pprime - 1/p);
    double A = std::log(kb) - B * (1/p - 1/Pref);
    if (B < 0) {
        throw SolutionError("B < 0 in outerTQ");
    }

    // solve
    vector<double> pp(ncomp);
    double maxdif = 1e10 * TOL;
    int itr = 0;
    double Rmin = 0, Rmax = 1;
    while (maxdif > TOL && itr < MAXITER) {
        // save previous values for calculating the difference at the end of the iteration
        vector<double> u_old = u;
        double A_old = A;

        double R0 = kb * Q / (kb * Q + kb0 * (1 - Q));
        double R = BoundedSecantInner(kb0, Q, u, x, cppargs, R0, Rmin, Rmax, DBL_EPSILON, 1e-8, 200);

        double pp_sum = 0;
        double eupp_sum = 0;
        for (int i = 0; i < ncomp; i++) {
            pp[i] = x[i] / (1 - R + kb0 * R * std::exp(u[i]));
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                pp_sum += pp[i];
                eupp_sum += std::exp(u[i]) * pp[i];
            }
        }
        kb = pp_sum / eupp_sum;

        p = 1 / (1 / Pref + (std::log(kb) - A) / B);
        for (int i = 0; i < ncomp; i++) {
            if (x_ions == 0) {
                xl[i] = pp[i] / pp_sum;
                xv[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else if (cppargs.z.empty() || cppargs.z[i] == 0) {
                xl[i] = pp[i] / pp_sum * (1 - x_ions/(1-Q));
                xv[i] = std::exp(u[i]) * pp[i] / eupp_sum;
            }
            else {
                xl[i] = x[i] / (1 - Q);
                xv[i] = 0;
            }
        }

        rhol = pcsaft_den_cpp(t, p, xl, 0, cppargs);
        fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, xl, cppargs);
        rhov = pcsaft_den_cpp(t, p, xv, 1, cppargs);
        fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, xv, cppargs);
        for (int i = 0; i < ncomp; i++) {
            k[i] = fugcoef_l[i] / fugcoef_v[i];
            u[i] = std::log(k[i] / kb);
        }

        if (itr == 0) {
            B = std::log(kbprime / kb) / (1/Pprime - 1/p);
            if (B < 0) {
                throw SolutionError("B < 0 in outerTQ");
            }
        }
        A = std::log(kb) - B * (1/p - 1/Pref);

        maxdif = std::abs(A - A_old);
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                double dif = std::abs(u[i] - u_old[i]);
                if (dif > maxdif) {
                    maxdif = dif;
                } else if (!std::isfinite(dif)) {
                    maxdif = dif;
                }
            }
        }

        itr += 1;
    }

    if (!std::isfinite(p) || !std::isfinite(maxdif) || maxdif > 0.1 || p < 0) {
        throw SolutionError("outerTQ did not converge to a solution");
    }

    vector<double> result;
    result.push_back(p);
    result.insert(result.end(), xl.begin(), xl.end());
    result.insert(result.end(), xv.begin(), xv.end());
    return result;
}

double resid_inner(double R, double kb0, double Q, vector<double> u, vector<double> x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    double error = 0;

    vector<double> pp(ncomp);
    double L = 0;
    for (int i = 0; i < ncomp; i++) {
        if (cppargs.z.empty() || cppargs.z[i] == 0) {
            pp[i] = x[i] / (1 - R + kb0 * R * exp(u[i]));
            L += pp[i];
        } else {
            L += x[i];
        }
    }
    L = (1 - R) * L;

    error = pow((L + Q - 1), 2.);
    return error;
}


double pcsaft_den_cpp(double t, double p, vector<double> x, int phase, const add_args &cppargs) {
    /**
    Solve for the molar density when temperature and pressure are given.

    Parameters
    ----------
    t : double
        Temperature (K)
    p : double
        Pressure (Pa)
    x : vector<double>, shape (n,)
        Mole fractions of each component. It has a length of n, where n is
        the number of components in the system.
    phase : int
        The phase for which the calculation is performed. Options: 0 (liquid),
        1 (vapor).
    cppargs : add_args
        A struct containing additional arguments that can be passed for
        use in PC-SAFT:

        m : vector<double>, shape (n,)
            Segment number for each component.
        s : vector<double>, shape (n,)
            Segment diameter for each component. For ions this is the diameter of
            the hydrated ion. Units of Angstrom.
        e : vector<double>, shape (n,)
            Dispersion energy of each component. For ions this is the dispersion
            energy of the hydrated ion. Units of K.
        k_ij : vector<double>, shape (n*n,)
            Binary interaction parameters between components in the mixture.
            (dimensions: ncomp x ncomp)
        e_assoc : vector<double>, shape (n,)
            Association energy of the associating components. For non associating
            compounds this is set to 0. Units of K.
        vol_a : vector<double>, shape (n,)
            Effective association volume of the associating components. For non
            associating compounds this is set to 0.
        dipm : vector<double>, shape (n,)
            Dipole moment of the polar components. For components where the dipole
            term is not used this is set to 0. Units of Debye.
        dip_num : vector<double>, shape (n,)
            The effective number of dipole functional groups on each component
            molecule. Some implementations use this as an adjustable parameter
            that is fit to data.
        z : vector<double>, shape (n,)
            Charge number of the ions
        dielc : double
            Dielectric constant of the medium to be used for electrolyte
            calculations.

    Returns
    -------
    rho : double
        Molar density (mol m^-3)
    */
    // suppress debug output during density root-finding to avoid repeated prints.
    int debug_flag = cppargs.debug;
    DebugFlagGuard debug_guard(debug_flag, 0);

    int ncomp = static_cast<int>(x.size());
    vector<double> scan_grid = build_density_scan_grid();
    vector<DensityScanPoint> scan_points;
    scan_points.reserve(scan_grid.size());
    for (double nu : scan_grid) {
        scan_points.push_back(evaluate_density_scan_point(nu, t, ncomp, x, p, cppargs));
    }

    vector<DensityBracket> coarse_brackets = find_density_brackets(scan_points);
    vector<DensityBracket> refined_brackets;
    for (const DensityBracket &coarse : coarse_brackets) {
        append_refined_density_brackets(coarse, t, ncomp, x, p, cppargs, refined_brackets);
    }

    if (refined_brackets.empty()) {
        throw SolutionError("No continuous density root brackets were found for the requested state.");
    }

    vector<DensityRootCandidate> candidates;
    candidates.reserve(refined_brackets.size());
    for (const DensityBracket &bracket : refined_brackets) {
        DensityRootCandidate candidate;
        candidate.rho_sort = reduced_to_molar(0.5 * (bracket.nu_lo + bracket.nu_hi), t, ncomp, x, cppargs);

        try {
            double rho_lo = reduced_to_molar(bracket.nu_lo, t, ncomp, x, cppargs);
            double rho_hi = reduced_to_molar(bracket.nu_hi, t, ncomp, x, cppargs);
            double rho_root = BrentRho(t, p, x, phase, cppargs, rho_lo, rho_hi, DBL_EPSILON, 1e-14, 200);
            validate_density_root(t, p, x, cppargs, rho_root, &candidate);
        }
        catch (const std::exception&) {
            candidate.valid = false;
        }

        candidates.push_back(candidate);
    }

    if (candidates.empty()) {
        throw SolutionError("Density solver did not produce any candidate roots.");
    }

    std::sort(candidates.begin(), candidates.end(), [](const DensityRootCandidate &a, const DensityRootCandidate &b) {
        return a.rho_sort < b.rho_sort;
    });

    const double rho_tol = 1e-8;
    if (phase == 1) {
        const double rho_extreme = candidates.front().rho_sort;
        DensityRootCandidate *best = nullptr;
        for (DensityRootCandidate &candidate : candidates) {
            if (std::abs(candidate.rho_sort - rho_extreme) > rho_tol * std::max(1.0, std::abs(rho_extreme))) {
                break;
            }
            if (candidate.valid && (best == nullptr || candidate.gres < best->gres)) {
                best = &candidate;
            }
        }
        if (best != nullptr) {
            return best->rho;
        }
        throw SolutionError("No valid density root found for vapor phase.");
    }

    const double rho_extreme = candidates.back().rho_sort;
    DensityRootCandidate *best = nullptr;
    for (auto it = candidates.rbegin(); it != candidates.rend(); ++it) {
        if (std::abs(it->rho_sort - rho_extreme) > rho_tol * std::max(1.0, std::abs(rho_extreme))) {
            break;
        }
        if (it->valid && (best == nullptr || it->gres < best->gres)) {
            best = &(*it);
        }
    }
    if (best != nullptr) {
        return best->rho;
    }
    throw SolutionError("No valid density root found for liquid phase.");
}

double estimate_flash_t(double p, double Q, vector<double> x, add_args cppargs) {
    /**
    Get a quick estimate of the temperature at which VLE occurs
    */
    double t_guess = _HUGE;
    int ncomp = static_cast<int>(x.size());

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (!cppargs.z.empty() && cppargs.z[i] != 0) {
            x_ions += x[i];
        }
    }

    bool guess_found = false;
    double t_step = 30;
    double t_start = 400;
    double t_lbound = 1;
    if (!cppargs.z.empty()) {
        t_step = 15;
        t_start = 350;
        t_lbound = 264;
    }
    while (!guess_found && t_start > t_lbound) {
        // initialize variables
        double Tprime = t_start - 50;
        double t = t_start;

        // calculate sigma for water, if it is present
        std::vector<double>::iterator water_iter = std::find(cppargs.e.begin(), cppargs.e.end(), 353.9449);
        int water_idx = -1;
        if (water_iter != cppargs.e.end()) {
            water_idx = static_cast<int>(std::distance(cppargs.e.begin(), water_iter));
            cppargs.s[water_idx] = calc_sigma(t, &calc_water_sigma);
        }

        try {
            double p1 = estimate_flash_p(t, Q, x, cppargs);
            double p2 = estimate_flash_p(Tprime, Q, x, cppargs);

            double slope = (std::log10(p1) - std::log10(p2)) / (1/t - 1/Tprime);
            double intercept = std::log10(p1) - slope * (1/t);
            t_guess = slope / (std::log10(p) - intercept);
            guess_found = true;
        } catch (const SolutionError&) {
            t_start -= t_step;
        }
    }

    if (!guess_found) {
        throw SolutionError("an estimate for the VLE temperature could not be found");
    }

    return t_guess;
}

double estimate_flash_p(double t, double Q, vector<double> x, const add_args &cppargs) {
    /**
    Get a quick estimate of the pressure at which VLE occurs
    */
    double p_guess = _HUGE;
    int ncomp = static_cast<int>(x.size());

    double x_ions = 0.; // overall mole fraction of ions in the system
    for (int i = 0; i < ncomp; i++) {
        if (!cppargs.z.empty() && cppargs.z[i] != 0) {
            x_ions += x[i];
        }
    }

    bool guess_found = false;
    double p_start = 10000;
    while (!guess_found && p_start < 1e7) {
        // initialize variables
        vector<double> fugcoef_l(ncomp), fugcoef_v(ncomp), k(ncomp), u(ncomp), kprime(ncomp);
        double rhol, rhov;
        double Pprime = 0.99 * p_start;
        double p = p_start;

        // calculate initial guess for compositions based on fugacity coefficients and Raoult's Law.
        rhol = pcsaft_den_cpp(t, p, x, 0, cppargs);
        rhov = pcsaft_den_cpp(t, p, x, 1, cppargs);
        if ((rhol - rhov) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, x, cppargs);
        fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, x, cppargs);

        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                k[i] = fugcoef_l[i] / fugcoef_v[i];
            } else {
                k[i] = 0; // set k to 0 for ionic components
            }
        }

        vector<double> xl(ncomp);
        vector<double> xv(ncomp);
        double xv_sum = 0;
        double xl_sum = 0;
        for (int i = 0; i < ncomp; i++) {
            xl[i] = x[i] / (1 + Q * (k[i] - 1));
            xl_sum += xl[i];
            xv[i] = k[i] * x[i] / (1 + Q * (k[i] - 1));
            xv_sum += xv[i];
        }

        if (xv_sum != 1) {
            for (int i = 0; i < ncomp; i++) {
                xv[i] = xv[i] / xv_sum;
            }
        }

        if (xl_sum != 1) {
            for (int i = 0; i < ncomp; i++) {
                xl[i] = xl[i] / xl_sum;
            }
        }

        rhol = pcsaft_den_cpp(t, p, xl, 0, cppargs);
        rhov = pcsaft_den_cpp(t, p, xv, 1, cppargs);
        if ((rhol - rhov) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, xl, cppargs);
        fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, xv, cppargs);
        double numer = 0;
        double denom = 0;
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                numer += xl[i] * fugcoef_l[i];
                denom += xv[i] * fugcoef_v[i];
            }
        }
        double ratio = numer / denom;

        rhol = pcsaft_den_cpp(t, Pprime, xl, 0, cppargs);
        rhov = pcsaft_den_cpp(t, Pprime, xv, 1, cppargs);
        if ((rhol - rhov) < 1e-4) {
            p_start = p_start + 2e5;
            continue;
        }
        fugcoef_l = pcsaft_fugcoef_cpp(t, rhol, xl, cppargs);
        fugcoef_v = pcsaft_fugcoef_cpp(t, rhov, xv, cppargs);
        numer = 0;
        denom = 0;
        for (int i = 0; i < ncomp; i++) {
            if (cppargs.z.empty() || cppargs.z[i] == 0) {
                numer += xl[i] * fugcoef_l[i];
                denom += xv[i] * fugcoef_v[i];
            }
        }
        double ratio_prime = numer / denom;

        double slope = (std::log10(ratio) - std::log10(ratio_prime)) / (std::log10(p) - std::log10(Pprime));
        double intercept = std::log10(ratio) - slope * std::log10(p);
        p_guess = pow(10, -intercept / slope);

        guess_found = true;
    }

    if (!guess_found) {
        throw SolutionError("an estimate for the VLE pressure could not be found");
    }

    return p_guess;
}


double reduced_to_molar(double nu, double t, int ncomp, vector<double> x, const add_args &cppargs) {
    vector<double> d(ncomp);
    double summ = 0.;
    for (int i = 0; i < ncomp; i++) {
        d[i] = cppargs.s[i]*(1-0.12*std::exp(-3*cppargs.e[i] / t));
        if (!cppargs.z.empty() && is_ion_species(cppargs, i)) {
            d[i] = compute_ion_diameter(i, t, cppargs);
        }
        summ += x[i]*cppargs.m[i]*pow(d[i],3.);
    }

    return 6/PI*nu/summ*1.0e30/N_AV;
}

double pcsaft_dielc_eps_cpp(vector<double> x, const add_args &cppargs) {
    DielcState state = evaluate_dielc_state(x, cppargs);
    return state.eps;
}

vector<double> pcsaft_dielc_diff_cpp(vector<double> x, const add_args &cppargs) {
    DielcState state = evaluate_dielc_state(x, cppargs);
    return state.deps_dx;
}

double dielc_water(double t) {
    /**
    Return the dielectric constant of water at the given temperature.

    t : double
        Temperature (K)

    This equation was fit to values given in the reference. For temperatures
    from 263.15 to 368.15 K values at 1 bar were used. For temperatures from
    368.15 to 443.15 K values at 10 bar were used.

    Reference:
    D. G. Archer and P. Wang, “The Dielectric Constant of Water and Debye‐Hückel
    Limiting Law Slopes,” J. Phys. Chem. Ref. Data, vol. 19, no. 2, pp. 371–411,
    Mar. 1990.
    */
    double dielc;
    if (t < 263.15) {
        throw ValueError("The current function for the dielectric constant for water is only valid for temperatures above 263.15 K.");
    }
    else if (t <= 368.15) {
        dielc = 7.6555618295E-04*t*t - 8.1783881423E-01*t + 2.5419616803E+02;
    }
    else if (t <= 443.15) {
        dielc = 0.0005003272124*t*t - 0.6285556029*t + 220.4467027;
    }
    else {
        throw ValueError("The current function for the dielectric constant for water is only valid for temperatures less than 443.15 K.");
    }
    return dielc;
}

double calc_water_sigma(double t) {
    return 3.8395 + 1.2828 * std::exp(-0.0074944 * t) - 1.3939 * std::exp(-0.00056029 * t);
}

add_args get_single_component(int i, const add_args &cppargs) {
    add_args args_single;
    args_single.born_model = cppargs.born_model;
    args_single.born_radius_model = cppargs.born_radius_model;
    args_single.born_diff_mode = cppargs.born_diff_mode;
    args_single.born_eps_mode = cppargs.born_eps_mode;
    args_single.DH_model = cppargs.DH_model;
    args_single.dielc_rule = cppargs.dielc_rule;
    args_single.dielc_diff_mode = cppargs.dielc_diff_mode;
    args_single.hc_dadx_diff_mode = cppargs.hc_dadx_diff_mode;
    args_single.disp_dadx_diff_mode = cppargs.disp_dadx_diff_mode;
    args_single.assoc_dadx_diff_mode = cppargs.assoc_dadx_diff_mode;
    args_single.polar_dadx_diff_mode = cppargs.polar_dadx_diff_mode;
    args_single.d_ion_mode = cppargs.d_ion_mode;
    args_single.mu_DH_diff_mode = cppargs.mu_DH_diff_mode;
    args_single.mu_DH_comp_dep_rel_perm = cppargs.mu_DH_comp_dep_rel_perm;
    args_single.mu_DH_include_sum_term = cppargs.mu_DH_include_sum_term;
    args_single.include_born_model = cppargs.include_born_model;
    args_single.d_born_mode = cppargs.d_born_mode;
    args_single.born_solvation_shell_model = cppargs.born_solvation_shell_model;
    args_single.born_dielectric_saturation = cppargs.born_dielectric_saturation;
    args_single.born_bulk_mode = cppargs.born_bulk_mode;
    args_single.mu_born_diff_mode = cppargs.mu_born_diff_mode;
    args_single.mu_born_comp_dep_rel_perm = cppargs.mu_born_comp_dep_rel_perm;
    args_single.mu_born_include_sum_term = cppargs.mu_born_include_sum_term;
    args_single.mu_born_comp_dep_delta_d = cppargs.mu_born_comp_dep_delta_d;
    args_single.mixed_rel_perm_water_index = (cppargs.mixed_rel_perm_water_index == i) ? 0 : -1;
    args_single.debug = cppargs.debug;
    args_single.m.push_back(cppargs.m[i]);
    args_single.s.push_back(cppargs.s[i]);
    args_single.e.push_back(cppargs.e[i]);
    if (!cppargs.e_assoc.empty()) {
        args_single.e_assoc.push_back(cppargs.e_assoc[i]);
        args_single.vol_a.push_back(cppargs.vol_a[i]);
    }
    if (!cppargs.dipm.empty()) {
        args_single.dipm.push_back(cppargs.dipm[i]);
        args_single.dip_num.push_back(cppargs.dip_num[i]);
    }
    if (!cppargs.z.empty()) {
        args_single.z.push_back(cppargs.z[i]);
        if (cppargs.dielc.size() > static_cast<size_t>(i)) {
            args_single.dielc.push_back(cppargs.dielc[i]);
        }
        if (cppargs.mw.size() > static_cast<size_t>(i)) {
            args_single.mw.push_back(cppargs.mw[i]);
        }
        if (cppargs.mixed_rel_perm_a.size() > static_cast<size_t>(i)) {
            args_single.mixed_rel_perm_a.push_back(cppargs.mixed_rel_perm_a[i]);
        }
        if (cppargs.mixed_rel_perm_b.size() > static_cast<size_t>(i)) {
            args_single.mixed_rel_perm_b.push_back(cppargs.mixed_rel_perm_b[i]);
        }
        if (cppargs.mixed_rel_perm_c.size() > static_cast<size_t>(i)) {
            args_single.mixed_rel_perm_c.push_back(cppargs.mixed_rel_perm_c[i]);
        }
        if (cppargs.mixed_rel_perm_mask.size() > static_cast<size_t>(i)) {
            args_single.mixed_rel_perm_mask.push_back(cppargs.mixed_rel_perm_mask[i]);
        }
        if (cppargs.d_born.size() > static_cast<size_t>(i)) {
            args_single.d_born.push_back(cppargs.d_born[i]);
        }
        if (cppargs.f_solv.size() > static_cast<size_t>(i)) {
            args_single.f_solv.push_back(cppargs.f_solv[i]);
        }
    }
    if (!cppargs.assoc_num.empty()) {
        args_single.assoc_num.push_back(cppargs.assoc_num[i]);

        if (args_single.assoc_num[0] > 0) {
            int nassoc = static_cast<int>(cppargs.assoc_num.size());
            int start = 0;
            for (int l = 0; l < static_cast<int>(cppargs.assoc_num.size()); l++) {
                if (l < i) {
                    start += 1;
                }
            }
            for (int j = 0; j < nassoc; j++) {
                for (int k = 0; k < args_single.assoc_num[0]; k++) {
                    args_single.assoc_matrix.push_back(cppargs.assoc_matrix[j*nassoc + start + k]);
                }
            }
        }
    }

    return args_single;
}

/*
----------------------------------------------------------------------------------------------------------------------
The code for the solvers was taken from CoolProp (https://github.com/CoolProp/CoolProp) and somewhat modified.
*/
/**

This function implements a 1-D bounded solver using the algorithm from Brent, R. P., Algorithms for Minimization Without Derivatives.
Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.

a and b must bound the solution of interest and f(a) and f(b) must have opposite signs.  If the function is continuous, there must be
at least one solution in the interval [a,b].

@param a The minimum bound for the solution of f=0
@param b The maximum bound for the solution of f=0
@param macheps The machine precision
@param tol_abs Tolerance (absolute)
@param maxiter Maximum number of steps allowed.  Will throw a SolutionError if the solution cannot be found
*/
double BrentRho(double t, double p, vector<double> x, int phase, const add_args &cppargs, double a, double b,
    double macheps, double tol_abs, int maxiter)
{
    int iter;
    double fa,fb,c,fc,m,tol,d,e,pp,q,s,r;
    fa = resid_rho(a, t, p, x, cppargs);
    fb = resid_rho(b, t, p, x, cppargs);

    // If one of the boundaries is to within tolerance, just stop
    if (std::abs(fb) < tol_abs) { return b;}
    if (std::isnan(fb)){
        throw ValueError("BrentRho's method f(b) is NAN for b");
    }
    if (std::abs(fa) < tol_abs) { return a;}
    if (std::isnan(fa)){
        throw ValueError("BrentRho's method f(a) is NAN for a");
    }
    if (fa*fb>0){
        throw ValueError("Inputs in BrentRho do not bracket the root");
    }

    c=a;
    fc=fa;
    iter=1;
    if (std::abs(fc)<std::abs(fb)){
        // Goto ext: from BrentRho ALGOL code
        a=b;
        b=c;
        c=a;
        fa=fb;
        fb=fc;
        fc=fa;
    }
    d=b-a;
    e=b-a;
    m=0.5*(c-b);
    tol=2*macheps*std::abs(b)+tol_abs;
    while (std::abs(m)>tol && fb!=0){
        // See if a bisection is forced
        if (std::abs(e)<tol || std::abs(fa) <= std::abs(fb)){
            m=0.5*(c-b);
            d=e=m;
        }
        else{
            s=fb/fa;
            if (a==c){
                //Linear interpolation
                pp=2*m*s;
                q=1-s;
            }
            else{
                //Inverse quadratic interpolation
                q=fa/fc;
                r=fb/fc;
                m=0.5*(c-b);
                pp=s*(2*m*q*(q-r)-(b-a)*(r-1));
                q=(q-1)*(r-1)*(s-1);
            }
            if (pp>0){
                q=-q;
            }
            else{
                pp=-pp;
            }
            s=e;
            e=d;
            m=0.5*(c-b);
            if (2*pp<3*m*q-std::abs(tol*q) || pp<std::abs(0.5*s*q)){
                d=pp/q;
            }
            else{
                m=0.5*(c-b);
                d=e=m;
            }
        }
        a=b;
        fa=fb;
        if (std::abs(d)>tol){
            b+=d;
        }
        else if (m>0){
            b+=tol;
        }
        else{
            b+=-tol;
        }
        fb=resid_rho(b, t, p, x, cppargs);
        if (std::isnan(fb)){
            throw ValueError("BrentRho's method f(t) is NAN for t");
        }
        if (std::abs(fb) < macheps){
            return b;
        }
        if (fb*fc>0){
            c=a;
            fc=fa;
            d=e=b-a;
        }
        if (std::abs(fc)<std::abs(fb)){
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        m=0.5*(c-b);
        tol=2*macheps*std::abs(b)+tol_abs;
        iter+=1;
        if (std::isnan(a)){
            throw ValueError("BrentRho's method a is NAN");}
        if (std::isnan(b)){
            throw ValueError("BrentRho's method b is NAN");}
        if (std::isnan(c)){
            throw ValueError("BrentRho's method c is NAN");}
        if (iter>maxiter){
            throw SolutionError("BrentRho's method reached maximum number of steps");}
        if (std::abs(fb)< 2*macheps*std::abs(b)){
            return b;
        }
    }
    return b;
}

double resid_rho(double rhomolar, double t, double p, vector<double> x, const add_args &cppargs){
    double peos = pcsaft_p_cpp(t, rhomolar, x, cppargs);
    double pressure_scale = std::max(std::abs(p), 1e-3);
    double cost = (peos-p)/pressure_scale;
    if (std::isfinite(cost)) {
        return cost;
    }
    else {
        return _HUGE;
    }
}

/**
In the secant function, a 1-D Newton-Raphson solver is implemented.  An initial guess for the solution is provided.

@param x0 The initial guess for the solution
@param xmax The upper bound for the solution
@param xmin The lower bound for the solution
@param dx The initial amount that is added to x in order to build the numerical derivative
@param tol The absolute value of the tolerance accepted for the objective function
@param maxiter Maximum number of iterations
@returns If no errors are found, the solution, otherwise the value _HUGE, the value for infinity
*/
double BoundedSecantInner(double kb0, double Q, vector<double> u, vector<double> x, const add_args &cppargs, double x0, double xmin,
    double xmax, double dx, double tol, int maxiter) {
    double x1=0,x2=0,x3=0,y1=0,y2=0,R,fval=999;
    int iter=1;
    if (std::abs(dx)==0){ throw ValueError("dx cannot be zero"); }
    while (std::abs(fval)>tol)
    {
        if (iter==1){
          x1 = x0;
          R = x1;
          x3 = R;
        }
        else if (iter==2){
          x2 = x0+dx;
          R = x2;
          x3 = R;
        }
        else {R=x2;}
        fval=resid_inner(R, kb0, Q, u, x, cppargs);

        if (iter==1){y1=fval;}
        else
        {
            if (std::isfinite(fval)) {
                y2 = fval;
            }
            else {
                y2 = 1e40;
            }
            x3=x2-y2/(y2-y1)*(x2-x1);
            // Check bounds, go half the way to the limit if limit is exceeded
            if (x3 < xmin)
            {
                x3 = (xmin + x2)/2;
            }
            if (x3 > xmax)
            {
                x3 = (xmax + x2)/2;
            }
            y1=y2; x1=x2; x2=x3;

        }
        if (iter>maxiter){
            throw SolutionError("BoundedSecant reached maximum number of iterations");
        }
        iter=iter+1;
    }
    return x3;
}

namespace {
vector<double> build_miac_gamma_vector(double t, double rho, const vector<double>& x, const add_args& cppargs)
{
    add_args args = cppargs;
    const int ncomp = static_cast<int>(x.size());
    if (args.z.empty() || std::all_of(args.z.begin(), args.z.end(), [](double v) { return std::abs(v) <= 1e-12; })) {
        throw ValueError("miac requires ionic species (non-zero z).");
    }
    if (args.mw.size() != x.size()) {
        throw ValueError("miac requires params['MW'] to be present and aligned with x.");
    }
    ChargeGroups groups = collect_charge_groups(args, x.size());
    if (groups.cations.empty() || groups.anions.empty()) {
        throw ValueError("miac requires at least one cation and one anion.");
    }
    if (groups.solvents.empty()) {
        throw ValueError("miac requires a neutral solvent reference.");
    }

    vector<double> fugcoef = pcsaft_fugcoef_cpp(t, rho, x, args);
    double p = pcsaft_p_cpp(t, rho, x, args);

    const double eps = 1e-12;
    vector<double> x_inf(ncomp, eps);
    vector<double> solvent_ref(groups.solvents.size(), 0.0);
    double solvent_sum = 0.0;
    for (size_t k = 0; k < groups.solvents.size(); ++k) {
        solvent_ref[k] = x[groups.solvents[k]];
        solvent_sum += solvent_ref[k];
    }
    if (solvent_sum <= 0.0) {
        throw ValueError("miac requires a positive solvent fraction.");
    }
    for (size_t k = 0; k < groups.solvents.size(); ++k) {
        x_inf[groups.solvents[k]] = solvent_ref[k] / solvent_sum;
    }
    double solvent_budget = std::max(1.0 - eps * static_cast<double>(ncomp - groups.solvents.size()), eps * static_cast<double>(groups.solvents.size()));
    for (size_t k = 0; k < groups.solvents.size(); ++k) {
        x_inf[groups.solvents[k]] *= solvent_budget;
    }
    double x_inf_sum = 0.0;
    for (double xi : x_inf) {
        x_inf_sum += xi;
    }
    for (double& xi : x_inf) {
        xi /= x_inf_sum;
    }

    double rho_inf = pcsaft_den_cpp(t, p, x_inf, 0, args);
    vector<double> fugcoef_inf = pcsaft_fugcoef_cpp(t, rho_inf, x_inf, args);
    vector<double> gamma_i(ncomp, 1.0);
    for (int i = 0; i < ncomp; ++i) {
        gamma_i[i] = fugcoef[i] / fugcoef_inf[i];
    }
    return gamma_i;
}

vector<double> build_gsolv_values(double t, double rho, const vector<double>& x, const add_args& cppargs)
{
    add_args args = cppargs;
    const int ncomp = static_cast<int>(x.size());
    if (args.z.empty() || std::all_of(args.z.begin(), args.z.end(), [](double v) { return std::abs(v) <= 1e-12; })) {
        throw ValueError("gsolv requires ionic species in params['z'].");
    }
    if (args.mw.size() != x.size()) {
        throw ValueError("gsolv requires params['MW'] to be present and aligned with x.");
    }
    ChargeGroups groups = collect_charge_groups(args, x.size());
    if (groups.cations.empty() && groups.anions.empty()) {
        throw ValueError("gsolv requires ionic species in params['z'].");
    }
    if (groups.solvents.empty()) {
        throw ValueError("gsolv requires at least one solvent species (z=0).");
    }

    vector<double> x_ref = x;
    vector<int> idx_ion = groups.cations;
    idx_ion.insert(idx_ion.end(), groups.anions.begin(), groups.anions.end());
    for (int i : idx_ion) {
        x_ref[i] = 0.0;
    }
    double solv_sum = 0.0;
    for (int i : groups.solvents) {
        solv_sum += x_ref[i];
    }
    if (solv_sum > 0.0) {
        for (int i : groups.solvents) {
            x_ref[i] /= solv_sum;
        }
    }
    else {
        double equal = 1.0 / static_cast<double>(groups.solvents.size());
        for (int i : groups.solvents) {
            x_ref[i] = equal;
        }
    }

    double p = pcsaft_p_cpp(t, rho, x_ref, args);
    int phase = (rho < 900.0) ? 1 : 0;
    vector<double> result(ncomp, 0.0);
    const double eps = 1e-12;
    for (int i : idx_ion) {
        vector<double> x_inf = x_ref;
        x_inf[i] = eps;
        double sum_inf = 0.0;
        for (double xi : x_inf) {
            sum_inf += xi;
        }
        for (double& xi : x_inf) {
            xi /= sum_inf;
        }
        double rho_inf = pcsaft_den_cpp(t, p, x_inf, phase, args);
        vector<double> lnfug_inf = pcsaft_lnfug_cpp(t, rho_inf, x_inf, args);
        result[i] = 8.31446261815324 * t * lnfug_inf[i];
    }
    return result;
}
} // namespace

PCSAFTMixtureNative::PCSAFTMixtureNative(const add_args& args)
    : args_(args)
{
}

const add_args& PCSAFTMixtureNative::args() const
{
    return args_;
}

size_t PCSAFTMixtureNative::ncomp() const
{
    return args_.m.size();
}

std::shared_ptr<PCSAFTStateNative> PCSAFTMixtureNative::state(double t, vector<double> x, int phase,
    bool has_p, double p, bool has_rho, double rho)
{
    return std::make_shared<PCSAFTStateNative>(shared_from_this(), t, std::move(x), phase, has_p, p, has_rho, rho);
}

PCSAFTStateNative::PCSAFTStateNative(std::shared_ptr<PCSAFTMixtureNative> mixture, double t, vector<double> x,
    int phase, bool has_p, double p, bool has_rho, double rho)
    : mixture_(std::move(mixture)), t_(t), x_(std::move(x)), phase_(phase),
      has_p_(has_p), has_rho_(has_rho), p_(p), rho_(rho),
      pressure_cached_(has_p), density_cached_(has_rho)
{
    if (!mixture_) {
        throw ValueError("PCSAFTStateNative requires a valid mixture.");
    }
    if (x_.size() != mixture_->ncomp()) {
        throw ValueError("State composition size does not match mixture size.");
    }
    if (phase_ != 0 && phase_ != 1) {
        throw ValueError("phase must be 0 (liquid) or 1 (vapor).");
    }
}

double PCSAFTStateNative::temperature() const
{
    return t_;
}

int PCSAFTStateNative::phase() const
{
    return phase_;
}

const vector<double>& PCSAFTStateNative::composition() const
{
    return x_;
}

double PCSAFTStateNative::pressure()
{
    if (pressure_cached_) {
        return p_;
    }
    if (!density_cached_) {
        throw ValueError("PCSAFTStateNative cannot compute pressure without density or pressure data.");
    }
    const add_args& args = mixture_->args();
    p_ = pcsaft_p_cpp(t_, rho_, x_, args);
    pressure_cached_ = true;
    return p_;
}

double PCSAFTStateNative::density()
{
    if (density_cached_) {
        return rho_;
    }
    if (!pressure_cached_) {
        throw ValueError("PCSAFTStateNative cannot compute density without pressure or density data.");
    }
    const add_args& args = mixture_->args();
    rho_ = pcsaft_den_cpp(t_, p_, x_, phase_, args);
    density_cached_ = true;
    return rho_;
}

double PCSAFTStateNative::Z()
{
    const add_args& args = mixture_->args();
    return pcsaft_Z_cpp(t_, density(), x_, args);
}

double PCSAFTStateNative::ares()
{
    const add_args& args = mixture_->args();
    return pcsaft_ares_cpp(t_, density(), x_, args);
}

double PCSAFTStateNative::dadt()
{
    const add_args& args = mixture_->args();
    return pcsaft_dadt_cpp(t_, density(), x_, args);
}

double PCSAFTStateNative::hres()
{
    const add_args& args = mixture_->args();
    return pcsaft_hres_cpp(t_, density(), x_, args);
}

double PCSAFTStateNative::sres()
{
    const add_args& args = mixture_->args();
    return pcsaft_sres_cpp(t_, density(), x_, args);
}

double PCSAFTStateNative::gres()
{
    const add_args& args = mixture_->args();
    return pcsaft_gres_cpp(t_, density(), x_, args);
}

vector<double> PCSAFTStateNative::lnfugcoef()
{
    const add_args& args = mixture_->args();
    return pcsaft_lnfug_cpp(t_, density(), x_, args);
}

vector<double> PCSAFTStateNative::fugcoef()
{
    const add_args& args = mixture_->args();
    return pcsaft_fugcoef_cpp(t_, density(), x_, args);
}

vector<double> PCSAFTStateNative::lnfugcoef_terms()
{
    const add_args& args = mixture_->args();
    return pcsaft_lnfug_terms_cpp(t_, density(), x_, args);
}

vector<double> PCSAFTStateNative::dielc_eval()
{
    const add_args& args = mixture_->args();
    vector<double> out;
    out.push_back(pcsaft_dielc_eps_cpp(x_, args));
    vector<double> deps = pcsaft_dielc_diff_cpp(x_, args);
    out.insert(out.end(), deps.begin(), deps.end());
    return out;
}

double PCSAFTStateNative::osmoticC()
{
    const add_args& args = mixture_->args();
    vector<double> xx = x_;

    int indx_solvent = -1;
    if (args.mw.size() == xx.size()) {
        for (size_t i = 0; i < xx.size(); ++i) {
            if (args.z.size() == xx.size() && std::abs(args.z[i]) < 1e-12) {
                indx_solvent = static_cast<int>(i);
                break;
            }
        }
    }
    if (indx_solvent < 0) {
        throw ValueError("PCSAFTStateNative::osmoticC requires a neutral solvent species.");
    }
    double mw_solvent = args.mw[indx_solvent];
    if (mw_solvent > 1.0) {
        mw_solvent /= 1000.0;
    }
    if (mw_solvent <= 0.0) {
        throw ValueError("Solvent molecular weight must be positive.");
    }
    if (xx[indx_solvent] <= 0.0) {
        throw ValueError("Solvent mole fraction must be positive.");
    }

    vector<double> molality(xx.size(), 0.0);
    double molality_sum = 0.0;
    for (size_t i = 0; i < xx.size(); ++i) {
        molality[i] = xx[i] / (xx[indx_solvent] * mw_solvent);
        if (static_cast<int>(i) != indx_solvent) {
            molality_sum += molality[i];
        }
    }
    if (molality_sum <= 0.0) {
        throw ValueError("Total molality is zero; osmotic coefficient is undefined.");
    }

    vector<double> x0(xx.size(), 0.0);
    x0[indx_solvent] = 1.0;

    double rho = density();
    double p = pcsaft_p_cpp(t_, rho, xx, args);
    int ph = (rho < 900.0) ? 1 : 0;
    double rho0 = pcsaft_den_cpp(t_, p, x0, ph, args);
    vector<double> fugcoef = pcsaft_fugcoef_cpp(t_, rho, xx, args);
    vector<double> fugcoef0 = pcsaft_fugcoef_cpp(t_, rho0, x0, args);
    double gamma = fugcoef[indx_solvent] / fugcoef0[indx_solvent];
    double osmC = -std::log(xx[indx_solvent] * gamma) / (mw_solvent * molality_sum);
    return osmC;
}

vector<double> PCSAFTStateNative::miac_m()
{
    const add_args& args_ref = mixture_->args();
    vector<double> gamma_i = build_miac_gamma_vector(t_, density(), x_, args_ref);
    ChargeGroups groups = collect_charge_groups(args_ref, x_.size());
    if (groups.cations.empty() || groups.anions.empty()) {
        throw ValueError("miac_m requires at least one cation and one anion.");
    }
    if (groups.solvents.empty()) {
        throw ValueError("miac_m requires a neutral solvent reference.");
    }
    if (args_ref.mw.size() != x_.size()) {
        throw ValueError("miac_m requires params['MW'] to be present and aligned with x.");
    }

    double mass_solvent = 0.0;
    for (int i : groups.solvents) {
        mass_solvent += x_[i] * args_ref.mw[i];
    }
    if (mass_solvent <= 0.0) {
        throw ValueError("Solvent mass is zero; check solvent mole fraction and MW.");
    }

    vector<double> mass_neutral(groups.solvents.size(), 0.0);
    for (size_t k = 0; k < groups.solvents.size(); ++k) {
        mass_neutral[k] = x_[groups.solvents[k]] * args_ref.mw[groups.solvents[k]];
    }
    double mass_neutral_sum = std::accumulate(mass_neutral.begin(), mass_neutral.end(), 0.0);
    if (mass_neutral_sum <= 0.0) {
        throw ValueError("Solvent mass is zero; check solvent mole fraction and MW.");
    }
    double w_sf_sum = 0.0;
    double M_solvent_mix_inv = 0.0;
    for (size_t k = 0; k < groups.solvents.size(); ++k) {
        double w_sf = mass_neutral[k] / mass_neutral_sum;
        w_sf_sum += w_sf;
        M_solvent_mix_inv += w_sf / args_ref.mw[groups.solvents[k]];
    }
    (void)w_sf_sum;
    if (M_solvent_mix_inv <= 0.0) {
        throw ValueError("Solvent molecular weight mixture is invalid.");
    }
    double M_solvent_mix = 1.0 / M_solvent_mix_inv;

    vector<double> out;
    out.reserve(static_cast<size_t>(groups.cations.size() * groups.anions.size()));
    for (int ic : groups.cations) {
        for (int ia : groups.anions) {
            int zc = static_cast<int>(std::round(std::abs(args_ref.z[ic])));
            int za = static_cast<int>(std::round(std::abs(args_ref.z[ia])));
            int g = gcd_int(zc, za);
            int nu_cat = za / g;
            int nu_an = zc / g;
            double n_salt = 0.5 * (x_[ic] / static_cast<double>(nu_cat) + x_[ia] / static_cast<double>(nu_an));
            double m_salt = n_salt / mass_solvent;
            double sum_nu = static_cast<double>(nu_cat + nu_an);
            double ln_gamma_pm = (static_cast<double>(nu_cat) * std::log(gamma_i[ic]) + static_cast<double>(nu_an) * std::log(gamma_i[ia])) / sum_nu;
            double gamma_pm_x = std::exp(ln_gamma_pm);
            double gamma_pm_m = gamma_pm_x / (1.0 + M_solvent_mix * m_salt * sum_nu);
            out.push_back(gamma_pm_m);
        }
    }
    return out;
}

vector<double> PCSAFTStateNative::miac()
{
    const add_args& args_ref = mixture_->args();
    vector<double> gamma_i = build_miac_gamma_vector(t_, density(), x_, args_ref);

    ChargeGroups groups = collect_charge_groups(args_ref, x_.size());
    if (groups.cations.empty() || groups.anions.empty()) {
        throw ValueError("miac requires at least one cation and one anion.");
    }
    if (groups.solvents.empty()) {
        throw ValueError("miac requires a neutral solvent reference.");
    }

    vector<double> out;
    out.reserve(static_cast<size_t>(groups.cations.size() * groups.anions.size()));
    for (int ic : groups.cations) {
        for (int ia : groups.anions) {
            int zc = static_cast<int>(std::round(std::abs(args_ref.z[ic])));
            int za = static_cast<int>(std::round(std::abs(args_ref.z[ia])));
            int g = gcd_int(zc, za);
            int nu_cat = za / g;
            int nu_an = zc / g;
            double sum_nu = static_cast<double>(nu_cat + nu_an);
            double ln_gamma_pm = (static_cast<double>(nu_cat) * std::log(gamma_i[ic]) + static_cast<double>(nu_an) * std::log(gamma_i[ia])) / sum_nu;
            out.push_back(std::exp(ln_gamma_pm));
        }
    }
    return out;
}

vector<double> PCSAFTStateNative::gsolv()
{
    const add_args& args = mixture_->args();
    return build_gsolv_values(t_, density(), x_, args);
}

FlashResultNative PCSAFTStateNative::flashTQ(double q, bool has_p_guess, double p_guess)
{
    const add_args& args = mixture_->args();
    vector<double> payload;
    if (has_p_guess) {
        payload = flashTQ_cpp(t_, q, x_, args, p_guess);
    } else {
        payload = flashTQ_cpp(t_, q, x_, args);
    }
    FlashResultNative result;
    result.value = payload.empty() ? 0.0 : payload[0];
    if (payload.size() > 1) {
        result.xl.assign(payload.begin() + 1, payload.begin() + 1 + static_cast<std::ptrdiff_t>(x_.size()));
        result.xv.assign(payload.begin() + 1 + static_cast<std::ptrdiff_t>(x_.size()), payload.end());
    }
    return result;
}

FlashResultNative PCSAFTStateNative::flashPQ(double p, double q, bool has_t_guess, double t_guess)
{
    const add_args& args = mixture_->args();
    vector<double> payload;
    if (has_t_guess) {
        payload = flashPQ_cpp(p, q, x_, args, t_guess);
    } else {
        payload = flashPQ_cpp(p, q, x_, args);
    }
    FlashResultNative result;
    result.value = payload.empty() ? 0.0 : payload[0];
    if (payload.size() > 1) {
        result.xl.assign(payload.begin() + 1, payload.begin() + 1 + static_cast<std::ptrdiff_t>(x_.size()));
        result.xv.assign(payload.begin() + 1 + static_cast<std::ptrdiff_t>(x_.size()), payload.end());
    }
    return result;
}

VaporizationResultNative PCSAFTStateNative::Hvap(bool has_p_guess, double p_guess)
{
    const add_args& args = mixture_->args();
    vector<double> payload;
    if (has_p_guess) {
        payload = flashTQ_cpp(t_, 0.0, x_, args, p_guess);
    } else {
        payload = flashTQ_cpp(t_, 0.0, x_, args);
    }
    double pvap = payload.empty() ? 0.0 : payload[0];
    double rho_l = pcsaft_den_cpp(t_, pvap, x_, 0, args);
    double hres_l = pcsaft_hres_cpp(t_, rho_l, x_, args);
    double rho_v = pcsaft_den_cpp(t_, pvap, x_, 1, args);
    double hres_v = pcsaft_hres_cpp(t_, rho_v, x_, args);
    VaporizationResultNative result;
    result.value = hres_v - hres_l;
    result.pressure = pvap;
    return result;
}

