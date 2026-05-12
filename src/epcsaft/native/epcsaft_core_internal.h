#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

#include "Eigen/Dense"

#include "epcsaft_electrolyte.h"
#include "epcsaft_autodiff_internal.h"
#include "autodiff/ad_derivative_result.h"

namespace thermo_detail {

enum class AresContributionKind {
    HC,
    DISP,
    ASSOC,
    ION,
    BORN
};

struct AresContributions {
    double hc = 0.0;
    double disp = 0.0;
    double assoc = 0.0;
    double ion = 0.0;
    double born = 0.0;
};

struct DadrhoResult {
    ScalarContributionTerms raw;
    ScalarContributionTerms terms;
};

struct MixtureState {
    vector<double> d;
    vector<double> dd_dt;
    vector<double> e_ij;
    vector<double> s_ij;
    double den = 0.0;
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

struct ChargeGroups {
    vector<int> cations;
    vector<int> anions;
    vector<int> solvents;
};

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

struct BornSSMDSData {
    vector<double> d_born;
    vector<double> D;
    vector<double> ddelta_prefac;
    vector<double> f_k;
    vector<double> bracket;
    double sum_bracket = 0.0;
    double sum_invD = 0.0;
    double sum_gap = 0.0;
    double sum_dpref_over_D2 = 0.0;
    double sum_bracket_dt = 0.0;
};

struct DielectricState {
    double eps = 0.0;
    vector<double> deps_dx;
};

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
template <size_t N>
double polynomial_value_cpp(const std::array<double, N> &coeffs, double x) {
    double value = 0.0;
    for (size_t i = 0; i < N; ++i) {
        value += coeffs[i] * std::pow(x, static_cast<int>(i));
    }
    return value;
}

template <size_t N>
double polynomial_derivative_cpp(const std::array<double, N> &coeffs, double x) {
    double value = 0.0;
    for (size_t i = 1; i < N; ++i) {
        value += coeffs[i] * static_cast<double>(i) * std::pow(x, static_cast<int>(i - 1));
    }
    return value;
}

template <size_t N>
double eta_weighted_derivative_cpp(const std::array<double, N> &coeffs, double x) {
    double value = 0.0;
    for (size_t i = 0; i < N; ++i) {
        value += coeffs[i] * static_cast<double>(i + 1) * std::pow(x, static_cast<int>(i));
    }
    return value;
}

namespace parameter_setup_detail {

inline double pair_sigma_cpp(size_t idx, int i, int j, const add_args &cppargs) {
    double sigma = 0.5 * (cppargs.s[i] + cppargs.s[j]);
    if (!cppargs.l_ij.empty()) {
        sigma *= (1.0 - cppargs.l_ij[idx]);
    }
    return sigma;
}

inline double pair_epsilon_cpp(size_t idx, int i, int j, const add_args &cppargs) {
    if (!cppargs.z.empty() && cppargs.z[i] * cppargs.z[j] > 0.0) {
        return 0.0;
    }
    double epsilon = std::sqrt(cppargs.e[i] * cppargs.e[j]);
    if (!cppargs.k_ij.empty()) {
        epsilon *= (1.0 - cppargs.k_ij[idx]);
    }
    return epsilon;
}

double pair_diameter_cpp(double d_i, double d_j);
double association_volume_cpp(int comp_i, int comp_j, int ncomp, const vector<double> &s_ij, const add_args &cppargs);
double ion_diameter_cpp(int i, double t, const add_args &cppargs);
double ion_diameter_cpp_dt(int i, double t, const add_args &cppargs);
double ion_born_radius_cpp(int i, double t, const add_args &cppargs);
double ion_born_radius_cpp_dt(int i, double t, const add_args &cppargs);

}  // namespace parameter_setup_detail

}  // namespace thermo_detail

using thermo_detail::AresContributionKind;
using thermo_detail::AresContributions;
using thermo_detail::AssociationSetup;
using thermo_detail::BornSSMDSData;
using thermo_detail::ChargeGroups;
using thermo_detail::DensityBracket;
using thermo_detail::DadrhoResult;
using thermo_detail::DensityRootCandidate;
using thermo_detail::DensityScanPoint;
using thermo_detail::DielectricState;
using thermo_detail::DispersionPolynomialState;
using thermo_detail::MixtureState;

ScalarContributionTerms make_scalar_terms(double hc, double disp, double assoc, double ion, double born, double total);
VectorContributionTerms make_vector_terms(
    const vector<double> &hc,
    const vector<double> &disp,
    const vector<double> &assoc,
    const vector<double> &ion,
    const vector<double> &born,
    const vector<double> &total
);
MixtureState mixture_state_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs, bool include_dt);
DispersionPolynomialState dispersion_polynomials_cpp(double m_avg, double eta);
double dh_kappa_cpp(double den, double t, double eps, double q2_sum);
double dh_chi_cpp(double kappa, double diameter);
AssociationSetup association_setup_cpp(const vector<double> &x, const add_args &cppargs, const vector<double> &s_ij, const vector<double> &ghs, double t);
ChargeGroups collect_charge_groups(const add_args& args, size_t ncomp);
double z_term_scale_cpp(const vector<double> &z_term, double target_total);
ScalarContributionTerms normalized_dadrho_terms_cpp(const ScalarContributionTerms &raw_terms);
ScalarContributionTerms compressibility_terms_from_dadrho_cpp(const DadrhoResult &result);
vector<double> density_scan_grid_cpp();
DensityScanPoint density_scan_point_cpp(double nu, double t, int ncomp, const vector<double> &x, double p, const add_args &cppargs);
std::string density_failure_message_cpp(const std::string &outcome, double t, double p, const vector<double> &x, int phase);
vector<DensityBracket> density_brackets_cpp(const vector<DensityScanPoint> &points);
void refine_density_brackets_cpp(
    const DensityBracket &coarse,
    double t,
    int ncomp,
    const vector<double> &x,
    double p,
    const add_args &cppargs,
    vector<DensityBracket> &refined_brackets
);
bool density_root_valid_cpp(double t, double p, const vector<double> &x, const add_args &cppargs, double rho, DensityRootCandidate *candidate);
DensitySolveResult density_solve_report_cpp(double t, double p, vector<double> x, int phase, const add_args &cppargs);
bool density_root_from_seed_cpp(
    double t,
    double p,
    const vector<double> &x,
    int phase,
    const add_args &cppargs,
    double rho_seed,
    DensityRootCandidate *candidate,
    double *rho_root_out
);
inline bool is_ion_species(const add_args &cppargs, int i) { return std::abs(cppargs.z[i]) > 1e-12; }
double reference_solvent_dielectric_constant_cpp(const vector<double> &x, const add_args &cppargs);
vector<double> reference_solvent_dielectric_derivative_cpp(const vector<double> &x, const add_args &cppargs);
BornSSMDSData born_shell_data_cpp(vector<double> x, const add_args &cppargs, double t, double eps_r, double eps_r_ion);
double ares_contribution_value_cpp(const AresContributions &terms, AresContributionKind kind);
AresContributions ares_contributions_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs);
epcsaft::native::autodiff::ADDerivativeResult cppad_eos_contribution_derivatives_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs);
epcsaft::native::autodiff::ADDerivativeResult cppad_pressure_density_derivative_cpp(double t, double rho, const vector<double> &x, const add_args &cppargs);
double dielectric_constant_rule_cpp(int rule, const vector<double> &x, const add_args &cppargs);
vector<double> dielectric_derivative_rule_cpp(int rule, const vector<double> &x, const add_args &cppargs);
AutoDual reference_solvent_dielectric_constant_ad_cpp(const vector<AutoDual> &x, const add_args &cppargs);
AutoDual dielectric_constant_rule_autodiff_cpp(int rule, const vector<AutoDual> &x, const add_args &cppargs);
vector<double> dielectric_derivative_rule_ad_cpp(int rule, const vector<double> &x, const add_args &cppargs);
DielectricState dielectric_state_cpp(const vector<double> &x, const add_args &cppargs);
DadrhoResult dadrho_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double hs_contact_density_derivative_cpp(double pair_diameter, double zeta2, double zeta3);
double hs_contact_time_derivative_cpp(
    double pair_diameter,
    double pair_diameter_dt,
    double zeta2,
    double zeta3,
    double dzeta2_dt,
    double dzeta3_dt
);
double hs_contact_composition_derivative_cpp(
    double pair_diameter,
    double zeta2,
    double zeta3,
    double dzeta2_dx,
    double dzeta3_dx
);
CompressibilityFactorResult compressibility_factor_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
ScalarContributionTerms residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
CompositionContributionResult composition_derivative_residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
ResidualChemicalPotentialResult residual_chemical_potential_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
FugacityContributionResult fugacity_coefficient_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
BornSSMDSDerivativeResult born_ssmds_liquid_derivatives_cpp(double t, double rho, int phase, vector<double> x, const add_args &cppargs);
ScalarContributionTerms temperature_derivative_residual_helmholtz_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
ActivityCoefficientNative activity_coefficient_values_cpp(
    ePCSAFTMixtureNative* mixture,
    double t,
    double rho,
    double p,
    int phase,
    const vector<double>& x,
    const add_args& args,
    const vector<int>& cation_indices,
    const vector<int>& anion_indices,
    const vector<int>& solvent_indices,
    const vector<int>& pair_cation_indices,
    const vector<int>& pair_anion_indices,
    const vector<int>& pair_nu_cation,
    const vector<int>& pair_nu_anion,
    bool include_aux,
    bool has_solvent_override,
    int solvent_override_index
);
