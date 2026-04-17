#include "epcsaft_contrib_internal.h"

using namespace thermo_detail;

double ion_born_radius_cpp(int i, double t, const add_args &cppargs) {
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
        return sigma_i * (1.0 - 0.12);
    }
    if (mode == 2) {
        return sigma_i * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
    }
    if (mode == 3) {
        if (cppargs.d_born.size() <= static_cast<size_t>(i) || cppargs.d_born[i] <= 0.0) {
            throw ValueError("d_Born_mode=fitted_param requires positive ionic params['d_born'] values.");
        }
        return cppargs.d_born[i];
    }
    throw ValueError("Unknown d_Born_mode. Supported values are 0, 1, 2, 3.");
}

double ion_born_radius_cpp_dt(int i, double t, const add_args &cppargs) {
    if (!is_ion_species(cppargs, i)) {
        return 0.0;
    }
    if (cppargs.d_born_mode == 2) {
        double sigma_i = cppargs.s[i];
        double expo = std::exp(-3.0 * cppargs.e[i] / t);
        return -0.36 * sigma_i * cppargs.e[i] * expo / (t * t);
    }
    return 0.0;
}

BornSSMDSData born_shell_data_cpp(vector<double> x, const add_args &cppargs, double t, double eps_r, double eps_r_ion) {
    int ncomp = static_cast<int>(x.size());
    const bool use_ssm = (cppargs.born_solvation_shell_model != 0);
    const bool use_ds = (cppargs.born_dielectric_saturation != 0);

    BornSSMDSData data;
    data.d_born.assign(ncomp, 1.0);
    data.D.assign(ncomp, 1.0);
    data.ddelta_prefac.assign(ncomp, 0.0);
    data.f_k.assign(ncomp, 1.0);
    data.bracket.assign(ncomp, 0.0);

    double f_mix = 0.0;
    for (int i = 0; i < ncomp; i++) {
        bool is_ion = is_ion_species(cppargs, i);
        double fi = 1.0;
        if (!is_ion && cppargs.f_solv.size() > static_cast<size_t>(i)) {
            fi = cppargs.f_solv[i];
        }
        data.f_k[i] = fi;
        f_mix += x[i] * fi;

        if (is_ion) {
            data.d_born[i] = ion_born_radius_cpp(i, t, cppargs);
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
            data.ddelta_prefac[i] = data.d_born[i] / std::abs(cppargs.z[i]);
        }
    }

    for (int i = 0; i < ncomp; i++) {
        bool is_ion = std::abs(cppargs.z[i]) > 1e-12;
        if (!is_ion) {
            data.D[i] = data.d_born[i];
            continue;
        }

        double delta_di = use_ssm ? ((f_mix - 1.0) * data.ddelta_prefac[i]) : 0.0;
        data.D[i] = data.d_born[i] + delta_di;
        if (data.D[i] <= 0.0) {
            throw ValueError("Born model generated a non-positive d_born + Delta d.");
        }

        double z2 = cppargs.z[i] * cppargs.z[i];
        double invD = 1.0 / data.D[i];
        double gap = (1.0 / data.d_born[i] - invD);
        double base_term = (1.0 - 1.0 / eps_r) * invD;
        double ds_term = use_ds ? ((1.0 - 1.0 / eps_r_ion) * gap) : 0.0;

        data.bracket[i] = base_term + ds_term;
        data.sum_bracket += x[i] * z2 * data.bracket[i];
        data.sum_invD += x[i] * z2 * invD;
        data.sum_gap += x[i] * z2 * gap;
        if (use_ssm) {
            data.sum_dpref_over_D2 += x[i] * z2 * data.ddelta_prefac[i] * invD * invD;
        }
    }
    return data;
}

double reference_solvent_dielectric_constant_cpp(const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
        return dielectric_constant_rule_cpp(cppargs.dielc_rule, x, cppargs);
    }
    double x_sol = 0.0;
    double eps_sol_num = 0.0;
    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) <= 1e-12) {
            x_sol += x[i];
            eps_sol_num += x[i] * cppargs.dielc[i];
        }
    }
    if (x_sol <= 0.0) {
        return dielectric_constant_rule_cpp(cppargs.dielc_rule, x, cppargs);
    }
    return eps_sol_num / x_sol;
}

vector<double> reference_solvent_dielectric_derivative_cpp(const vector<double> &x, const add_args &cppargs) {
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
            eps_sol_num += x[i] * cppargs.dielc[i];
        }
    }
    if (x_sol <= 0.0) {
        return deps;
    }
    double inv_xsol2 = 1.0 / (x_sol * x_sol);
    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) <= 1e-12) {
            deps[i] = (cppargs.dielc[i] * x_sol - eps_sol_num) * inv_xsol2;
        }
    }
    return deps;
}

AutoDual reference_solvent_dielectric_constant_ad_cpp(const vector<AutoDual> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    if (cppargs.z.size() != static_cast<size_t>(ncomp)) {
        vector<double> xd(ncomp, 0.0);
        for (int i = 0; i < ncomp; ++i) {
            xd[i] = scalar_value(x[i]);
        }
        return make_autodiff_scalar(dielectric_constant_rule_cpp(cppargs.dielc_rule, xd, cppargs), 0.0);
    }
    AutoDual x_sol = make_autodiff_scalar(0.0, 0.0);
    AutoDual eps_sol_num = make_autodiff_scalar(0.0, 0.0);
    for (int i = 0; i < ncomp; i++) {
        if (std::abs(cppargs.z[i]) <= 1e-12) {
            x_sol += x[i];
            eps_sol_num += x[i] * cppargs.dielc[i];
        }
    }
    if (scalar_value(x_sol) <= 0.0) {
        vector<double> xd(ncomp, 0.0);
        for (int i = 0; i < ncomp; ++i) {
            xd[i] = scalar_value(x[i]);
        }
        return make_autodiff_scalar(dielectric_constant_rule_cpp(cppargs.dielc_rule, xd, cppargs), 0.0);
    }
    return eps_sol_num / x_sol;
}

BornIntermediateState born_intermediate_state_cpp(
    double t,
    const vector<double> &x,
    const add_args &cppargs,
    bool include_dt,
    bool include_dx
) {
    BornIntermediateState state;
    if (cppargs.z.empty() || cppargs.born_model == 0) {
        return state;
    }
    state.model = cppargs.born_model;

    state.eps_value = dielectric_constant_rule_cpp(cppargs.dielc_rule, x, cppargs);
    state.deps_dx = dielectric_derivative_rule_cpp(cppargs.dielc_rule, x, cppargs);
    if (cppargs.born_eps_mode == 1) {
        state.eps_value = reference_solvent_dielectric_constant_cpp(x, cppargs);
        state.deps_dx = reference_solvent_dielectric_derivative_cpp(x, cppargs);
    }

    if (cppargs.born_model == 1) {
        for (int i = 0; i < static_cast<int>(x.size()); ++i) {
            if (is_ion_species(cppargs, i)) {
                double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                state.charge_radius_sum += x[i] * cppargs.z[i] * cppargs.z[i] / d_born_i;
                if (include_dt) {
                    double d_born_dt = ion_born_radius_cpp_dt(i, t, cppargs);
                    state.charge_radius_sum_dt += x[i] * cppargs.z[i] * cppargs.z[i] * (-d_born_dt) / (d_born_i * d_born_i);
                }
            }
        }
        return state;
    }

    if (cppargs.born_model == 2) {
        const double eps_r_ion = 8.0;
        state.shell = born_shell_data_cpp(x, cppargs, t, state.eps_value, eps_r_ion);
        return state;
    }

    throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
}

// EqID: born_ares_dadrho
double dadrho_born_cpp() {
    return 0.0;
}

// EqID: born_ares_dT
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

ContributionDadxResult dadx_born_cpp(const BornIntermediateState &born_state, double t, double rho, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    ContributionDadxResult result;
    result.dadx.assign(ncomp, 0.0);
    if (cppargs.z.empty()) {
        return result;
    }

    if (born_state.model == 1) {
        double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        result.ares = -Kborn * (1.0 - 1.0 / born_state.eps_value) * born_state.charge_radius_sum;

        if (cppargs.born_diff_mode == 1) {
            result.dadx = contribution_dadx_fd_cpp(AresContributionKind::BORN, t, rho, x, cppargs, result.ares);
        } else if (cppargs.born_diff_mode == 4) {
            result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::BORN, t, rho, x, cppargs);
        } else {
            for (int i = 0; i < ncomp; ++i) {
                double ion_part = 0.0;
                if (is_ion_species(cppargs, i)) {
                    double d_born_i = ion_born_radius_cpp(i, t, cppargs);
                    ion_part = (1.0 - 1.0 / born_state.eps_value) * cppargs.z[i] * cppargs.z[i] / d_born_i;
                }
                double eps_part = 0.0;
                if (cppargs.born_diff_mode == 2) {
                    eps_part = born_state.deps_dx[i] / (born_state.eps_value * born_state.eps_value);
                } else if (cppargs.born_diff_mode == 3) {
                    eps_part = 0.0;
                } else {
                    eps_part = born_state.charge_radius_sum * born_state.deps_dx[i] / (born_state.eps_value * born_state.eps_value);
                }
                result.dadx[i] = -Kborn * (ion_part + eps_part);
            }
        }
    } else if (born_state.model == 2) {
        const double eps_r_ion = 8.0;
        const double Kborn = E_CHRG * E_CHRG / (4.0 * PI * kb * t * perm_vac);
        result.ares = -Kborn * born_state.shell.sum_bracket;

        if (cppargs.born_diff_mode == 1) {
            result.dadx = contribution_dadx_fd_cpp(AresContributionKind::BORN, t, rho, x, cppargs, result.ares);
        } else if (cppargs.born_diff_mode == 4) {
            result.dadx = contribution_dadx_autodiff_cpp(AresContributionKind::BORN, t, rho, x, cppargs);
        } else {
            const double inv_eps2 = 1.0 / (born_state.eps_value * born_state.eps_value);
            const double shell_coeff = 1.0 / eps_r_ion - 1.0 / born_state.eps_value;
            const bool use_deps = (cppargs.mu_born_comp_dep_rel_perm != 0);
            const bool use_shell_chain = (cppargs.mu_born_comp_dep_delta_d != 0);
            const double deps_multiplier = (cppargs.mu_born_include_sum_term != 0) ? born_state.shell.sum_gap : 1.0;
            for (int k = 0; k < ncomp; ++k) {
                double direct_part = 0.0;
                if (std::abs(cppargs.z[k]) > 1e-12) {
                    direct_part = cppargs.z[k] * cppargs.z[k] * born_state.shell.bracket[k];
                }
                double deps_part = use_deps ? deps_multiplier * born_state.deps_dx[k] * inv_eps2 : 0.0;
                double ddelta_part = use_shell_chain ? shell_coeff * born_state.shell.sum_dpref_over_D2 * born_state.shell.f_k[k] : 0.0;
                result.dadx[k] = -Kborn * (direct_part + deps_part + ddelta_part);
            }
        }
    } else if (born_state.model != 0) {
        throw ValueError("Unknown born_model. Supported values are 0, 1, 2.");
    }

    for (int i = 0; i < ncomp; ++i) {
        result.sum_x_dadx += x[i] * result.dadx[i];
    }
    return result;
}
