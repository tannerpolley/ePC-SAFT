#include "epcsaft_core_internal.h"

using namespace thermo_detail;

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
