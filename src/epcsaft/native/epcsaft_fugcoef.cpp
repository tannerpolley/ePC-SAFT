#include "epcsaft_core_internal.h"

using namespace thermo_detail;

// EqID: lnphi_alpha_near_ideal
double stable_logz_over_zminus1(double Z) {
    double dz = Z - 1.0;
    if (std::abs(dz) < 1e-8) {
        double dz2 = dz * dz;
        double dz3 = dz2 * dz;
        return 1.0 - 0.5 * dz + dz2 / 3.0 - 0.25 * dz3;
    }
    return std::log(Z) / dz;
}

// EqID: lnphi_total_sum
// EqID: lnphi_alpha
VectorContributionTerms ln_fugacity_terms_from_mu_z_cpp(
    const VectorContributionTerms &mu_terms,
    const ScalarContributionTerms &z_raw_terms
) {
    vector<double> z_terms = {
        z_raw_terms.hc,
        z_raw_terms.disp,
        z_raw_terms.polar,
        z_raw_terms.assoc,
        z_raw_terms.ion,
        z_raw_terms.born,
    };
    double z_scale = z_term_scale_cpp(z_terms, z_raw_terms.total);
    double z_weight = stable_logz_over_zminus1(1.0 + z_raw_terms.total);
    double z_correction_scale = z_scale * z_weight;

    auto subtract_z = [z_correction_scale](const vector<double> &mu, double z_value) {
        vector<double> out(mu.size(), 0.0);
        for (int i = 0; i < static_cast<int>(mu.size()); ++i) {
            out[i] = mu[i] - z_value * z_correction_scale;
        }
        return out;
    };

    vector<double> lnfug_hc = subtract_z(mu_terms.hc, z_raw_terms.hc);
    vector<double> lnfug_disp = subtract_z(mu_terms.disp, z_raw_terms.disp);
    vector<double> lnfug_polar = subtract_z(mu_terms.polar, z_raw_terms.polar);
    vector<double> lnfug_assoc = subtract_z(mu_terms.assoc, z_raw_terms.assoc);
    vector<double> lnfug_ion = subtract_z(mu_terms.ion, z_raw_terms.ion);
    vector<double> lnfug_born = subtract_z(mu_terms.born, z_raw_terms.born);
    vector<double> lnfug_total(mu_terms.total.size(), 0.0);
    for (int i = 0; i < static_cast<int>(lnfug_total.size()); ++i) {
        lnfug_total[i] = lnfug_hc[i] + lnfug_disp[i] + lnfug_polar[i] + lnfug_assoc[i] + lnfug_ion[i] + lnfug_born[i];
    }
    return make_vector_terms(lnfug_hc, lnfug_disp, lnfug_polar, lnfug_assoc, lnfug_ion, lnfug_born, lnfug_total);
}

// EqID: lnphi_total
// EqID: lnphi_total_sum
FugacityContributionResult fugacity_coefficient_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    ResidualChemicalPotentialResult mu_result = residual_chemical_potential_result_cpp(t, rho, std::move(x), cppargs);
    FugacityContributionResult result;
    result.mu = mu_result.mu;
    result.composition = mu_result.composition;
    result.lnfugcoef = ln_fugacity_terms_from_mu_z_cpp(mu_result.mu, mu_result.composition.z_raw);
    return result;
}

vector<double> lnfug_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return fugacity_coefficient_result_cpp(t, rho, std::move(x), cppargs).lnfugcoef.total;
}

vector<double> fugcoef_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return exp_vector(lnfug_cpp(t, rho, std::move(x), cppargs));
}
