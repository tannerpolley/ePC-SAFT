#include "epcsaft_contrib_internal.h"

using thermo_detail::AssociationIntermediateState;
using thermo_detail::AssociationSetup;
using thermo_detail::ContributionDadxResult;
using thermo_detail::HardChainState;
using thermo_detail::MixtureState;
using thermo_detail::parameter_setup_detail::association_volume_cpp;
using thermo_detail::parameter_setup_detail::pair_diameter_cpp;

namespace assoc_detail {

// EqID: x_assoc_site
static vector<double> association_site_fractions_cpp(vector<double> XA_guess, vector<double> delta_ij, double den, vector<double> x) {
    int num_sites = static_cast<int>(XA_guess.size());
    vector<double> XA = XA_guess;

    int idxij = -1;
    for (int i = 0; i < num_sites; ++i) {
        double summ = 0.0;
        for (int j = 0; j < num_sites; ++j) {
            idxij += 1;
            summ += den * x[j] * XA_guess[j] * delta_ij[idxij];
        }
        XA[i] = 1.0 / (1.0 + summ);
    }

    return XA;
}

// EqID: dx_assoc_drho
static vector<double> association_site_fraction_dt_cpp(vector<double> delta_ij, double den, vector<double> XA, vector<double> ddelta_dt, vector<double> x) {
    int num_sites = static_cast<int>(XA.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_sites, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites, num_sites);

    int ij = 0;
    for (int i = 0; i < num_sites; ++i) {
        double summ = 0.0;
        for (int j = 0; j < num_sites; ++j) {
            B(i) -= x[j] * XA[j] * ddelta_dt[ij];
            A(i, j) = x[j] * delta_ij[ij];
            summ += x[j] * XA[j] * delta_ij[ij];
            ij += 1;
        }
        A(i, i) = std::pow(1.0 + den * summ, 2.0) / den;
    }

    Eigen::MatrixXd solution = A.lu().solve(B);
    vector<double> dXA_dt(num_sites);
    for (int i = 0; i < num_sites; ++i) {
        dXA_dt[i] = solution(i);
    }
    return dXA_dt;
}

// EqID: dx_assoc_dxk
static vector<double> association_site_fraction_dx_cpp(vector<int> assoc_num, vector<double> delta_ij, double den, vector<double> XA, vector<double> ddelta_dx, vector<double> x) {
    int num_sites = static_cast<int>(XA.size());
    int ncomp = static_cast<int>(assoc_num.size());
    Eigen::MatrixXd B(num_sites * ncomp, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites * ncomp, num_sites * ncomp);

    int idx1 = 0;
    int ij = 0;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < num_sites; ++j) {
            double sum1 = 0.0;
            for (int k = 0; k < num_sites; ++k) {
                sum1 += den * x[k] * (XA[k] * ddelta_dx[i * num_sites * num_sites + j * num_sites + k]);
                A(ij, i * num_sites + k) = XA[j] * XA[j] * den * x[k] * delta_ij[j * num_sites + k];
            }

            double sum2 = 0.0;
            for (int l = 0; l < assoc_num[i]; ++l) {
                sum2 += XA[idx1 + l] * delta_ij[idx1 * num_sites + l * num_sites + j];
            }

            A(ij, ij) = A(ij, ij) + 1.0;
            B(ij) = -XA[j] * XA[j] * (sum1 + sum2);
            ij += 1;
        }
        idx1 += assoc_num[i];
    }

    Eigen::MatrixXd solution = A.lu().solve(B);
    vector<double> dXA_dx(num_sites * ncomp);
    for (int i = 0; i < num_sites * ncomp; ++i) {
        dXA_dx[i] = solution(i);
    }
    return dXA_dx;
}

static vector<double> solve_association_site_fractions_cpp(const vector<double> &delta_ij, double den, const vector<double> &x_assoc) {
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
        XA = association_site_fractions_cpp(XA_old, delta_ij, den, x_assoc);
        dif = 0.0;
        for (int i = 0; i < num_sites; ++i) {
            dif += std::abs(XA[i] - XA_old[i]);
        }
        for (int i = 0; i < num_sites; ++i) {
            XA_old[i] = 0.5 * (XA[i] + XA_old[i]);
        }
    }
    return XA;
}

}  // namespace assoc_detail

// EqID: rho_j_assoc
// EqID: delta_assoc
AssociationSetup association_setup_cpp(const vector<double> &x, const add_args &cppargs, const vector<double> &s_ij, const vector<double> &ghs, double t) {
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
        for (int j = 0; j < num_sites; ++j) {
            int comp_j = setup.site_component_index[j];
            if (cppargs.assoc_matrix[idxa] != 0) {
                double eABij = 0.5 * (cppargs.e_assoc[comp_i] + cppargs.e_assoc[comp_j]);
                double volABij = association_volume_cpp(comp_i, comp_j, ncomp, s_ij, cppargs);
                setup.delta_ij[idxa] = ghs[comp_i * ncomp + comp_j] * (std::exp(eABij / t) - 1.0)
                    * std::pow(s_ij[comp_i * ncomp + comp_j], 3.0) * volABij;
            }
            ++idxa;
        }
    }

    return setup;
}

namespace assoc_detail {

// EqID: ddelta_assoc_drho
static vector<double> association_site_fraction_density_terms_cpp(
    const vector<double> &delta_ij,
    double den,
    const vector<double> &XA,
    const vector<double> &ddelta_weighted,
    const vector<double> &x_assoc
) {
    int num_sites = static_cast<int>(XA.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_sites, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites, num_sites);

    int ij = 0;
    for (int i = 0; i < num_sites; ++i) {
        double summ = 0.0;
        for (int j = 0; j < num_sites; ++j) {
            B(i) -= x_assoc[j] * XA[j] * ddelta_weighted[ij];
            A(i, j) = x_assoc[j] * delta_ij[ij];
            summ += x_assoc[j] * XA[j] * delta_ij[ij];
            ++ij;
        }
        B(i) -= summ;
        A(i, i) = std::pow(1.0 + den * summ, 2.0) / den;
    }

    Eigen::MatrixXd solution = A.lu().solve(B);
    vector<double> dXA_weighted(num_sites, 0.0);
    for (int i = 0; i < num_sites; ++i) {
        dXA_weighted[i] = solution(i);
    }
    return dXA_weighted;
}

// EqID: ddelta_assoc_dxk
static vector<double> association_site_fraction_composition_terms_cpp(
    const vector<double> &delta_ij,
    double den,
    const vector<double> &XA,
    const vector<double> &ddelta_dx,
    const vector<int> &site_component_index,
    const vector<double> &x_assoc,
    int ncomp
) {
    int num_sites = static_cast<int>(XA.size());
    vector<double> dXA_dx(ncomp * num_sites, 0.0);

    for (int k = 0; k < ncomp; ++k) {
        Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_sites, 1);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites, num_sites);

        int ij = 0;
        for (int i = 0; i < num_sites; ++i) {
            double direct_sum = 0.0;
            double delta_sum = 0.0;
            for (int j = 0; j < num_sites; ++j) {
                if (site_component_index[j] == k) {
                    direct_sum += XA[j] * delta_ij[ij];
                }
                delta_sum += x_assoc[j] * XA[j] * ddelta_dx[k * num_sites * num_sites + ij];
                A(i, j) = x_assoc[j] * delta_ij[ij];
                ++ij;
            }
            B(i) = -(direct_sum + delta_sum);
            A(i, i) += 1.0 / (den * XA[i] * XA[i]);
        }

        Eigen::MatrixXd solution = A.lu().solve(B);
        for (int i = 0; i < num_sites; ++i) {
            dXA_dx[k * num_sites + i] = solution(i);
        }
    }

    return dXA_dx;
}

}  // namespace assoc_detail

AssociationIntermediateState association_intermediate_state_cpp(
    const MixtureState &thermo,
    const HardChainState &hc_state,
    double t,
    const vector<double> &x,
    const add_args &cppargs,
    bool include_dt,
    bool include_dx,
    const vector<double> *dghs_dt
) {
    AssociationIntermediateState state;
    if (cppargs.e_assoc.empty()) {
        return state;
    }
    state.active = true;
    state.setup = association_setup_cpp(x, cppargs, thermo.s_ij, hc_state.ghs, t);
    const vector<int> &site_component_index = state.setup.site_component_index;
    const vector<double> &x_assoc = state.setup.x_assoc;
    const vector<double> &delta_ij = state.setup.delta_ij;
    int ncomp = static_cast<int>(x.size());
    int num_sites = static_cast<int>(site_component_index.size());

    state.XA = assoc_detail::solve_association_site_fractions_cpp(delta_ij, thermo.den, x_assoc);

    if (include_dt) {
        vector<double> ddelta_dt(num_sites * num_sites, 0.0);
        if (dghs_dt == nullptr) {
            throw ValueError("Association temperature derivatives require externally provided hard-sphere contact time derivatives.");
        }

        for (int i = 0; i < num_sites; ++i) {
            for (int j = 0; j < num_sites; ++j) {
                if (cppargs.assoc_matrix[i * num_sites + j] != 0) {
                    double eABij = 0.5 * (cppargs.e_assoc[site_component_index[i]] + cppargs.e_assoc[site_component_index[j]]);
                    double volABij = association_volume_cpp(site_component_index[i], site_component_index[j], ncomp, thermo.s_ij, cppargs);
                    double pair_diameter = pair_diameter_cpp(
                        thermo.d[site_component_index[i]],
                        thermo.d[site_component_index[j]]
                    );
                    ddelta_dt[i * num_sites + j] = std::pow(pair_diameter, 3) * volABij * (
                        -eABij / std::pow(t, 2) * std::exp(eABij / t) * hc_state.ghs[site_component_index[i] * ncomp + site_component_index[j]]
                        + (*dghs_dt)[site_component_index[i] * ncomp + site_component_index[j]] * (std::exp(eABij / t) - 1.0)
                    );
                }
            }
        }

        state.dXA_dt = assoc_detail::association_site_fraction_dt_cpp(delta_ij, thermo.den, state.XA, ddelta_dt, x_assoc);
    }

    if (include_dx) {
        vector<double> ddelta_dx(num_sites * num_sites * ncomp, 0.0);
        int idx_ddelta = 0;
        for (int k = 0; k < ncomp; ++k) {
            for (int i = 0; i < num_sites; ++i) {
                for (int j = 0; j < num_sites; ++j) {
                    if (cppargs.assoc_matrix[i * num_sites + j] != 0) {
                        double pair_diameter = pair_diameter_cpp(thermo.d[site_component_index[i]], thermo.d[site_component_index[j]]);
                        double dzeta2_dx = PI / 6.0 * thermo.den * cppargs.m[k] * std::pow(thermo.d[k], 2);
                        double dzeta3_dx = PI / 6.0 * thermo.den * cppargs.m[k] * std::pow(thermo.d[k], 3);
                        double dghsd_dx = hs_contact_composition_derivative_cpp(
                            pair_diameter,
                            hc_state.zeta[2],
                            hc_state.zeta[3],
                            dzeta2_dx,
                            dzeta3_dx
                        );
                        double eABij = 0.5 * (cppargs.e_assoc[site_component_index[i]] + cppargs.e_assoc[site_component_index[j]]);
                        double volABij = association_volume_cpp(site_component_index[i], site_component_index[j], ncomp, thermo.s_ij, cppargs);
                        ddelta_dx[idx_ddelta] = dghsd_dx * (std::exp(eABij / t) - 1.0)
                            * std::pow(thermo.s_ij[site_component_index[i] * ncomp + site_component_index[j]], 3) * volABij;
                    }
                    ++idx_ddelta;
                }
            }
        }

        state.dXA_dx = assoc_detail::association_site_fraction_dx_cpp(cppargs.assoc_num, delta_ij, thermo.den, state.XA, ddelta_dx, x_assoc);
    }

    return state;
}

// EqID: assoc_ares_dadrho
double dadrho_assoc_cpp(
    const MixtureState &thermo,
    const HardChainState &hc_state,
    const AssociationIntermediateState &assoc_state,
    const vector<double> &x,
    const add_args &cppargs,
    double t
) {
    if (!assoc_state.active) {
        return 0.0;
    }

    const vector<int> &site_component_index = assoc_state.setup.site_component_index;
    const vector<double> &x_assoc = assoc_state.setup.x_assoc;
    const vector<double> &delta_ij = assoc_state.setup.delta_ij;
    int ncomp = static_cast<int>(x.size());
    int num_sites = static_cast<int>(site_component_index.size());
    vector<double> ddelta_weighted(num_sites * num_sites, 0.0);

    int ij = 0;
    for (int i = 0; i < num_sites; ++i) {
        int comp_i = site_component_index[i];
        for (int j = 0; j < num_sites; ++j) {
            int comp_j = site_component_index[j];
            if (cppargs.assoc_matrix[ij] != 0) {
                double pair_diameter = pair_diameter_cpp(thermo.d[comp_i], thermo.d[comp_j]);
                double eABij = 0.5 * (cppargs.e_assoc[comp_i] + cppargs.e_assoc[comp_j]);
                double volABij = association_volume_cpp(comp_i, comp_j, ncomp, thermo.s_ij, cppargs);
                ddelta_weighted[ij] = hs_contact_density_derivative_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3])
                    * (std::exp(eABij / t) - 1.0)
                    * std::pow(thermo.s_ij[comp_i * ncomp + comp_j], 3.0)
                    * volABij;
            }
            ++ij;
        }
    }

    vector<double> dXA_weighted = assoc_detail::association_site_fraction_density_terms_cpp(
        delta_ij, thermo.den, assoc_state.XA, ddelta_weighted, x_assoc
    );

    double value = 0.0;
    for (int i = 0; i < num_sites; ++i) {
        int component_index = site_component_index[i];
        value += x[component_index] * (1.0 / assoc_state.XA[i] - 0.5) * dXA_weighted[i];
    }
    return value;
}

// EqID: assoc_ares_dT
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

// EqID: assoc_ares_dxk
ContributionDadxResult dadx_assoc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const AssociationIntermediateState &assoc_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);
    if (!assoc_state.active) {
        return result;
    }

    int num_sites = static_cast<int>(assoc_state.setup.site_component_index.size());
    vector<double> ddelta_dx(num_sites * num_sites * ncomp, 0.0);
    int idx_ddelta = 0;
    for (int k = 0; k < ncomp; ++k) {
        for (int i = 0; i < num_sites; ++i) {
            int comp_i = assoc_state.setup.site_component_index[i];
            for (int j = 0; j < num_sites; ++j) {
                int comp_j = assoc_state.setup.site_component_index[j];
                if (cppargs.assoc_matrix[i * num_sites + j] != 0) {
                    double pair_diameter = pair_diameter_cpp(thermo.d[comp_i], thermo.d[comp_j]);
                    double dzeta2_dx = PI / 6.0 * thermo.den * cppargs.m[k] * std::pow(thermo.d[k], 2);
                    double dzeta3_dx = PI / 6.0 * thermo.den * cppargs.m[k] * std::pow(thermo.d[k], 3);
                    double dghsd_dx = hs_contact_composition_derivative_cpp(
                        pair_diameter,
                        hc_state.zeta[2],
                        hc_state.zeta[3],
                        dzeta2_dx,
                        dzeta3_dx
                    );
                    double eABij = 0.5 * (cppargs.e_assoc[comp_i] + cppargs.e_assoc[comp_j]);
                    double volABij = association_volume_cpp(comp_i, comp_j, ncomp, thermo.s_ij, cppargs);
                    ddelta_dx[idx_ddelta] = dghsd_dx * (std::exp(eABij / t) - 1.0)
                        * std::pow(thermo.s_ij[comp_i * ncomp + comp_j], 3.0) * volABij;
                }
                ++idx_ddelta;
            }
        }
    }

    vector<double> dXA_dx = assoc_detail::association_site_fraction_composition_terms_cpp(
        assoc_state.setup.delta_ij,
        thermo.den,
        assoc_state.XA,
        ddelta_dx,
        assoc_state.setup.site_component_index,
        assoc_state.setup.x_assoc,
        ncomp
    );

    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < num_sites; ++j) {
            result.dadx[i] += x[assoc_state.setup.site_component_index[j]] * dXA_dx[i * num_sites + j] * (1.0 / assoc_state.XA[j] - 0.5);
        }
    }

    for (int i = 0; i < num_sites; ++i) {
        int component_index = assoc_state.setup.site_component_index[i];
        result.dadx[component_index] += std::log(assoc_state.XA[i]) - 0.5 * assoc_state.XA[i] + 0.5;
        result.ares += x[component_index] * (std::log(assoc_state.XA[i]) - 0.5 * assoc_state.XA[i] + 0.5);
    }

    if (cppargs.assoc_dadx_diff_mode == 1) {
        throw ValueError("unsupported: association composition derivative backend is not enabled.");
    } else if (cppargs.assoc_dadx_diff_mode == 2) {
        result.dadx = contribution_dadx_cppad_cpp(AresContributionKind::ASSOC, t, rho, x, cppargs);
    }

    for (int i = 0; i < ncomp; ++i) {
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    return result;
}
