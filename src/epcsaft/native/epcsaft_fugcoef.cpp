#include "epcsaft_core_internal.h"

namespace fugcoef_detail {

static vector<double> exp_vector(const vector<double> &values) {
    vector<double> out(values.size(), 0.0);
    for (int i = 0; i < static_cast<int>(values.size()); ++i) {
        out[i] = std::exp(values[i]);
    }
    return out;
}

// EqID: lnphi_alpha_near_ideal
static double stable_logz_over_zminus1(double Z) {
    double dz = Z - 1.0;
    if (std::abs(dz) < 1e-8) {
        double dz2 = dz * dz;
        double dz3 = dz2 * dz;
        return 1.0 - 0.5 * dz + dz2 / 3.0 - 0.25 * dz3;
    }
    return std::log(Z) / dz;
}

static double lnfug_correction_scale_cpp(const ScalarContributionTerms &z_raw_terms) {
    vector<double> z_terms = {
        z_raw_terms.hc,
        z_raw_terms.disp,
        z_raw_terms.assoc,
        z_raw_terms.ion,
        z_raw_terms.born,
    };
    double z_scale = z_term_scale_cpp(z_terms, z_raw_terms.total);
    double z_weight = stable_logz_over_zminus1(1.0 + z_raw_terms.total);
    return z_scale * z_weight;
}

// EqID: lnphi_alpha
static vector<double> lnfug_contribution_cpp(
    const vector<double> &mu_term,
    double z_value,
    double z_correction_scale
) {
    vector<double> out(mu_term.size(), 0.0);
    for (int i = 0; i < static_cast<int>(mu_term.size()); ++i) {
        out[i] = mu_term[i] - z_value * z_correction_scale;
    }
    return out;
}

static vector<double> lnfug_hc_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.hc, z_raw_terms.hc, z_correction_scale);
}

static vector<double> lnfug_disp_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.disp, z_raw_terms.disp, z_correction_scale);
}

static vector<double> lnfug_assoc_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.assoc, z_raw_terms.assoc, z_correction_scale);
}

static vector<double> lnfug_ion_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.ion, z_raw_terms.ion, z_correction_scale);
}

static vector<double> lnfug_born_cpp(const VectorContributionTerms &mu_terms, const ScalarContributionTerms &z_raw_terms, double z_correction_scale) {
    return lnfug_contribution_cpp(mu_terms.born, z_raw_terms.born, z_correction_scale);
}

static vector<double> lnfug_total_cpp(
    const vector<double> &lnfug_hc,
    const vector<double> &lnfug_disp,
    const vector<double> &lnfug_assoc,
    const vector<double> &lnfug_ion,
    const vector<double> &lnfug_born
) {
    vector<double> total(lnfug_hc.size(), 0.0);
    for (int i = 0; i < static_cast<int>(total.size()); ++i) {
        total[i] = lnfug_hc[i] + lnfug_disp[i] + lnfug_assoc[i] + lnfug_ion[i] + lnfug_born[i];
    }
    return total;
}

}  // namespace fugcoef_detail

// EqID: lnphi_total
// EqID: lnphi_total_sum
FugacityContributionResult fugacity_coefficient_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    ResidualChemicalPotentialResult mu_result = residual_chemical_potential_result_cpp(t, rho, std::move(x), cppargs);
    double z_correction_scale = fugcoef_detail::lnfug_correction_scale_cpp(mu_result.composition.z_raw);

    vector<double> lnfug_hc = fugcoef_detail::lnfug_hc_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_disp = fugcoef_detail::lnfug_disp_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_assoc = fugcoef_detail::lnfug_assoc_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_ion = fugcoef_detail::lnfug_ion_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_born = fugcoef_detail::lnfug_born_cpp(mu_result.mu, mu_result.composition.z_raw, z_correction_scale);
    vector<double> lnfug_total = fugcoef_detail::lnfug_total_cpp(lnfug_hc, lnfug_disp, lnfug_assoc, lnfug_ion, lnfug_born);

    FugacityContributionResult result;
    result.mu = mu_result.mu;
    result.composition = mu_result.composition;
    result.lnfugcoef = make_vector_terms(lnfug_hc, lnfug_disp, lnfug_assoc, lnfug_ion, lnfug_born, lnfug_total);
    return result;
}

vector<double> lnfug_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return fugacity_coefficient_result_cpp(t, rho, std::move(x), cppargs).lnfugcoef.total;
}

vector<double> fugcoef_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    return fugcoef_detail::exp_vector(lnfug_cpp(t, rho, std::move(x), cppargs));
}
