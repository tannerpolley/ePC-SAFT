#include "epcsaft_core_internal.h"

using namespace thermo_detail;

namespace {

// EqID: mu_res_from_ares
vector<double> mu_contribution_cpp(
    double ares_term,
    double z_term,
    const vector<double> &dadx_term,
    double sum_x_term
) {
    vector<double> out(dadx_term.size(), 0.0);
    for (int i = 0; i < static_cast<int>(dadx_term.size()); ++i) {
        out[i] = ares_term + z_term + dadx_term[i] - sum_x_term;
    }
    return out;
}

vector<double> mu_hc_cpp(const CompositionContributionResult &composition) {
    return mu_contribution_cpp(composition.ares.hc, composition.z_raw.hc, composition.dadx.hc, composition.sum_x_dadx.hc);
}

vector<double> mu_disp_cpp(const CompositionContributionResult &composition) {
    return mu_contribution_cpp(composition.ares.disp, composition.z_raw.disp, composition.dadx.disp, composition.sum_x_dadx.disp);
}

vector<double> mu_assoc_cpp(const CompositionContributionResult &composition) {
    return mu_contribution_cpp(composition.ares.assoc, composition.z_raw.assoc, composition.dadx.assoc, composition.sum_x_dadx.assoc);
}

vector<double> mu_ion_cpp(const CompositionContributionResult &composition) {
    return mu_contribution_cpp(composition.ares.ion, composition.z_raw.ion, composition.dadx.ion, composition.sum_x_dadx.ion);
}

vector<double> mu_born_cpp(const CompositionContributionResult &composition) {
    return mu_contribution_cpp(composition.ares.born, composition.z_raw.born, composition.dadx.born, composition.sum_x_dadx.born);
}

vector<double> mu_total_cpp(
    const vector<double> &mu_hc,
    const vector<double> &mu_disp,
    const vector<double> &mu_assoc,
    const vector<double> &mu_ion,
    const vector<double> &mu_born
) {
    vector<double> total(mu_hc.size(), 0.0);
    for (int i = 0; i < static_cast<int>(total.size()); ++i) {
        total[i] = mu_hc[i] + mu_disp[i] + mu_assoc[i] + mu_ion[i] + mu_born[i];
    }
    return total;
}

}  // namespace

// EqID: mu_res
// EqID: mu_res_sum
ResidualChemicalPotentialResult residual_chemical_potential_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    CompositionContributionResult composition = composition_derivative_residual_helmholtz_result_cpp(t, rho, std::move(x), cppargs);

    vector<double> mu_hc = mu_hc_cpp(composition);
    vector<double> mu_disp = mu_disp_cpp(composition);
    vector<double> mu_assoc = mu_assoc_cpp(composition);
    vector<double> mu_ion = mu_ion_cpp(composition);
    vector<double> mu_born = mu_born_cpp(composition);
    vector<double> mu_total = mu_total_cpp(mu_hc, mu_disp, mu_assoc, mu_ion, mu_born);

    ResidualChemicalPotentialResult result;
    result.mu = make_vector_terms(mu_hc, mu_disp, mu_assoc, mu_ion, mu_born, mu_total);
    result.composition = std::move(composition);
    return result;
}

vector<double> mures_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return residual_chemical_potential_result_cpp(t, rho, std::move(x), cppargs).mu.total;
}
