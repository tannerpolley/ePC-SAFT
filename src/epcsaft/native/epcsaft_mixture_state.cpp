#include "epcsaft_core_internal.h"
#include "contributions/epcsaft_contrib_internal.h"

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

MixtureState mixture_state_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs, bool include_dt) {
    MixtureState state;
    int ncomp = static_cast<int>(x.size());
    state.d.assign(ncomp, 0.0);
    if (include_dt) {
        state.dd_dt.assign(ncomp, 0.0);
    }
    state.e_ij.assign(ncomp * ncomp, 0.0);
    state.s_ij.assign(ncomp * ncomp, 0.0);
    state.den = rho * N_AV / 1.0e30;

    for (int i = 0; i < ncomp; ++i) {
        state.d[i] = cppargs.s[i] * (1.0 - 0.12 * std::exp(-3.0 * cppargs.e[i] / t));
        if (include_dt) {
            state.dd_dt[i] = -0.36 * cppargs.s[i] * cppargs.e[i] * std::exp(-3.0 * cppargs.e[i] / t) / (t * t);
        }
        if (!cppargs.z.empty() && std::abs(cppargs.z[i]) > 1e-12) {
            state.d[i] = ion_diameter_cpp(i, t, cppargs);
            if (include_dt) {
                state.dd_dt[i] = ion_diameter_cpp_dt(i, t, cppargs);
            }
        }
    }

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
        }
    }

    return state;
}

