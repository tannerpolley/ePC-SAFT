#include "epcsaft_core_internal.h"

using namespace thermo_detail;

ScalarContributionTerms make_scalar_terms(double hc, double disp, double polar, double assoc, double ion, double born, double total) {
    ScalarContributionTerms out;
    out.hc = hc;
    out.disp = disp;
    out.polar = polar;
    out.assoc = assoc;
    out.ion = ion;
    out.born = born;
    out.total = total;
    return out;
}

VectorContributionTerms make_vector_terms(
    const vector<double> &hc,
    const vector<double> &disp,
    const vector<double> &polar,
    const vector<double> &assoc,
    const vector<double> &ion,
    const vector<double> &born,
    const vector<double> &total
) {
    VectorContributionTerms out;
    out.hc = hc;
    out.disp = disp;
    out.polar = polar;
    out.assoc = assoc;
    out.ion = ion;
    out.born = born;
    out.total = total;
    return out;
}

vector<double> exp_vector(const vector<double> &values) {
    vector<double> out(values.size(), 0.0);
    for (int i = 0; i < static_cast<int>(values.size()); ++i) {
        out[i] = std::exp(values[i]);
    }
    return out;
}

vector<double> solve_association_site_fractions_cpp(
    const vector<double> &delta_ij,
    double den,
    const vector<double> &x_assoc
) {
    int num_sites = static_cast<int>(x_assoc.size());
    vector<double> XA(num_sites, 0.0);
    for (int i = 0; i < num_sites; ++i) {
        XA[i] = (-1.0 + std::sqrt(1.0 + 8.0 * den * delta_ij[i * num_sites + i]))
            / (4.0 * den * delta_ij[i * num_sites + i]);
        if (!std::isfinite(XA[i])) {
            XA[i] = 0.02;
        }
    }

    int ctr = 0;
    double dif = 1000.0;
    vector<double> XA_old = XA;
    while ((ctr < 100) && (dif > 1e-15)) {
        ctr += 1;
        XA = ::association_site_fractions_cpp(XA_old, delta_ij, den, x_assoc);
        dif = 0.0;
        for (int i = 0; i < num_sites; ++i) {
            dif += std::abs(XA[i] - XA_old[i]);
        }
        for (int i = 0; i < num_sites; ++i) {
            XA_old[i] = (XA[i] + XA_old[i]) / 2.0;
        }
    }
    return XA;
}

double pair_sigma_cpp(size_t idx, int i, int j, const add_args &cppargs) {
    double sigma = 0.5 * (cppargs.s[i] + cppargs.s[j]);
    if (!cppargs.l_ij.empty()) {
        sigma *= (1.0 - cppargs.l_ij[idx]);
    }
    return sigma;
}

double pair_epsilon_cpp(size_t idx, int i, int j, const add_args &cppargs) {
    if (!cppargs.z.empty() && cppargs.z[i] * cppargs.z[j] > 0.0) {
        return 0.0;
    }
    double epsilon = std::sqrt(cppargs.e[i] * cppargs.e[j]);
    if (!cppargs.k_ij.empty()) {
        epsilon *= (1.0 - cppargs.k_ij[idx]);
    }
    return epsilon;
}

double pair_diameter_cpp(double d_i, double d_j) {
    return d_i * d_j / (d_i + d_j);
}

double hs_contact_value_cpp(double pair_diameter, double zeta2, double zeta3) {
    return 1.0 / (1.0 - zeta3)
        + pair_diameter * 3.0 * zeta2 / std::pow(1.0 - zeta3, 2.0)
        + std::pow(pair_diameter, 2.0) * 2.0 * zeta2 * zeta2 / std::pow(1.0 - zeta3, 3.0);
}

double hs_contact_density_derivative_cpp(double pair_diameter, double zeta2, double zeta3) {
    return zeta3 / std::pow(1.0 - zeta3, 2.0)
        + pair_diameter * (3.0 * zeta2 / std::pow(1.0 - zeta3, 2.0) + 6.0 * zeta2 * zeta3 / std::pow(1.0 - zeta3, 3.0))
        + std::pow(pair_diameter, 2.0) * (4.0 * zeta2 * zeta2 / std::pow(1.0 - zeta3, 3.0) + 6.0 * zeta2 * zeta2 * zeta3 / std::pow(1.0 - zeta3, 4.0));
}

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

ThermoCommonState thermo_common_state_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs, bool include_dt) {
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
        if (include_dt) {
            state.dd_dt[i] = -0.36 * cppargs.s[i] * cppargs.e[i] * std::exp(-3.0 * cppargs.e[i] / t) / (t * t);
        }
        if (!cppargs.z.empty() && std::abs(cppargs.z[i]) > 1e-12) {
            int mode = cppargs.d_ion_mode;
            double sigma_i = cppargs.s[i];
            if (sigma_i <= 0.0) {
                throw ValueError("DH/ion diameter requires positive ionic sigma_i.");
            }
            if (mode == 0) {
                state.d[i] = sigma_i;
                if (include_dt) {
                    state.dd_dt[i] = 0.0;
                }
            }
            else if (mode == 1) {
                state.d[i] = sigma_i * (1.0 - 0.12);
                if (include_dt) {
                    state.dd_dt[i] = 0.0;
                }
            }
            else if (mode == 2) {
                double expo = std::exp(-3.0 * cppargs.e[i] / t);
                state.d[i] = sigma_i * (1.0 - 0.12 * expo);
                if (include_dt) {
                    state.dd_dt[i] = -0.36 * sigma_i * cppargs.e[i] * expo / (t * t);
                }
            }
            else {
                throw ValueError("Unknown d_ion_mode. Supported values are 0, 1, 2.");
            }
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
        for (int j = 0; j < ncomp; ++j) {
            idx += 1;
            state.s_ij[idx] = pair_sigma_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.e_ij[idx] = pair_epsilon_cpp(static_cast<size_t>(idx), i, j, cppargs);
            state.m2es3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * state.e_ij[idx] / t * std::pow(state.s_ij[idx], 3);
            state.m2e2s3 += x[i] * x[j] * cppargs.m[i] * cppargs.m[j] * std::pow(state.e_ij[idx] / t, 2) * std::pow(state.s_ij[idx], 3);

            double pair_diameter = pair_diameter_cpp(state.d[i], state.d[j]);
            state.ghs[idx] = hs_contact_value_cpp(pair_diameter, state.zeta[2], state.zeta[3]);
            state.denghs[idx] = hs_contact_density_derivative_cpp(pair_diameter, state.zeta[2], state.zeta[3]);
        }
    }

    return state;
}

DispersionPolynomialState dispersion_polynomials_cpp(double m_avg, double eta) {
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

vector<double> dipole_coefficients_cpp(const std::array<double, 5> &c0, const std::array<double, 5> &c1, const std::array<double, 5> &c2, double m) {
    vector<double> coeffs(5, 0.0);
    double a = (m - 1.0) / m;
    double b = (m - 2.0) / m;
    for (size_t i = 0; i < coeffs.size(); ++i) {
        coeffs[i] = c0[i] + a * c1[i] + a * b * c2[i];
    }
    return coeffs;
}

double dh_kappa_cpp(double den, double t, double eps, double q2_sum) {
    return std::sqrt(den * E_CHRG * E_CHRG / kb / t / (eps * perm_vac) * q2_sum);
}

double dh_chi_cpp(double kappa, double diameter) {
    double ka = kappa * diameter;
    return 3.0 / std::pow(ka, 3.0) * (1.5 + std::log(1.0 + ka) - 2.0 * (1.0 + ka) + 0.5 * std::pow(1.0 + ka, 2.0));
}

AssociationSetup association_setup_cpp(
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

void build_charge_metadata_cpp(
    const add_args& args,
    bool& has_ionic,
    vector<int>& cation_indices,
    vector<int>& anion_indices,
    vector<int>& solvent_indices,
    vector<int>& pair_cation_indices,
    vector<int>& pair_anion_indices,
    vector<int>& pair_nu_cation,
    vector<int>& pair_nu_anion
) {
    has_ionic = false;
    cation_indices.clear();
    anion_indices.clear();
    solvent_indices.clear();
    pair_cation_indices.clear();
    pair_anion_indices.clear();
    pair_nu_cation.clear();
    pair_nu_anion.clear();

    if (args.z.empty()) {
        return;
    }
    ChargeGroups groups = collect_charge_groups(args, args.z.size());
    cation_indices = groups.cations;
    anion_indices = groups.anions;
    solvent_indices = groups.solvents;
    has_ionic = (!cation_indices.empty() || !anion_indices.empty());

    for (int ic : cation_indices) {
        for (int ia : anion_indices) {
            int zc = static_cast<int>(std::round(std::abs(args.z[ic])));
            int za = static_cast<int>(std::round(std::abs(args.z[ia])));
            int g = gcd_int(zc, za);
            pair_cation_indices.push_back(ic);
            pair_anion_indices.push_back(ia);
            pair_nu_cation.push_back(za / g);
            pair_nu_anion.push_back(zc / g);
        }
    }
}

