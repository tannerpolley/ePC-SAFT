#include <vector>
#include <string>
#include <memory>
#include <limits>
#include <map>

using std::vector;

const static double kb = 1.380648465952442093e-23; // Boltzmann constant, J K^-1
const static double PI = 3.141592653589793;
const static double N_AV = 6.022140857e23; // Avagadro's number
const static double E_CHRG = 1.6021766208e-19; // elementary charge, units of coulomb
const static double perm_vac = 8.854187817e-22; //permittivity in vacuum, C V^-1 Angstrom^-1

#ifndef DBL_EPSILON
    #define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif

const static double HUGE_DBL = std::numeric_limits<double>::infinity();

struct add_args {
    vector<double> m;
    vector<double> s;
    vector<double> e;
    vector<double> k_ij;
    vector<double> e_assoc;
    vector<double> vol_a;
    vector<double> z;
    vector<double> dielc;
    vector<double> mw;
    vector<double> mixed_rel_perm_a;
    vector<double> mixed_rel_perm_b;
    vector<double> mixed_rel_perm_c;
    vector<int> mixed_rel_perm_mask;
    int mixed_rel_perm_water_index;
    int dielc_rule;
    int dielc_diff_mode;
    int hc_dadx_diff_mode;
    int disp_dadx_diff_mode;
    int assoc_dadx_diff_mode;
    int d_ion_mode;
    int mu_DH_diff_mode;
    int mu_DH_comp_dep_rel_perm;
    int mu_DH_include_sum_term;
    int include_born_model;
    int d_born_mode;
    int born_solvation_shell_model;
    int born_dielectric_saturation;
    int born_bulk_mode;
    int mu_born_diff_mode;
    int mu_born_comp_dep_rel_perm;
    int mu_born_include_sum_term;
    int mu_born_comp_dep_delta_d;
    vector<double> d_born;
    vector<double> f_solv;
    int born_model;
    int born_radius_model;
    int born_diff_mode;
    int born_eps_mode;
    int DH_model;
    vector<int> assoc_num;
    vector<int> assoc_matrix;
    vector<double> k_hb;
    vector<double> l_ij;
};

struct ActivityCoefficientNative {
    vector<double> component_activity_coefficients;
    vector<double> mean_ionic_activity_coefficients_mole_fraction;
    vector<double> mean_ionic_activity_coefficients_molality;
    vector<double> solvation_free_energy;
    vector<double> pair_molality;
    vector<double> pair_conversion_factor;
    vector<int> cation_indices;
    vector<int> anion_indices;
    vector<int> solvent_indices;
    vector<int> pair_cation_indices;
    vector<int> pair_anion_indices;
    vector<int> pair_nu_cation;
    vector<int> pair_nu_anion;
    int solvent_index;
    double osmotic_coefficient;
};

struct ScalarContributionTerms {
    double hc = 0.0;
    double disp = 0.0;
    double assoc = 0.0;
    double ion = 0.0;
    double born = 0.0;
    double total = 0.0;
};

struct CompressibilityFactorResult {
    ScalarContributionTerms raw;
    ScalarContributionTerms terms;
};

struct VectorContributionTerms {
    vector<double> hc;
    vector<double> disp;
    vector<double> assoc;
    vector<double> ion;
    vector<double> born;
    vector<double> total;
};

struct CompositionContributionResult {
    VectorContributionTerms dadx;
    ScalarContributionTerms ares;
    ScalarContributionTerms sum_x_dadx;
    ScalarContributionTerms z_raw;
    ScalarContributionTerms z;
    std::map<std::string, std::string> derivative_backend;
    bool derivative_available = true;
    std::string not_available_reason;
};

struct ResidualChemicalPotentialResult {
    VectorContributionTerms mu;
    CompositionContributionResult composition;
};

struct FugacityContributionResult {
    VectorContributionTerms mu;
    VectorContributionTerms lnfugcoef;
    CompositionContributionResult composition;
};

struct BornSSMDSDerivativeResult {
    bool supported = false;
    std::string backend = "not_available";
    std::string message;
    vector<double> a_born_d_d_born;
    vector<double> a_born_d_f_solv;
    vector<double> mu_res_d_d_born_row_major;
    vector<double> mu_res_d_f_solv_row_major;
    vector<double> lnfug_d_d_born_row_major;
    vector<double> lnfug_d_f_solv_row_major;
    vector<double> lngamma_d_d_born_row_major;
    vector<double> lngamma_d_f_solv_row_major;
    int ncomp = 0;
};

struct ReferenceStateKey {
    double t = 0.0;
    double p = 0.0;
    int phase = 0;
    vector<double> x_ref;
};

struct ReferenceStateValue {
    double rho = 0.0;
    vector<double> fugcoef;
};

struct ReferenceStateCacheEntry {
    ReferenceStateKey key;
    ReferenceStateValue value;
};

struct DensityCandidateDiagnostics {
    double rho_sort = 0.0;
    double rho = 0.0;
    double gres = 1.0e300;
    double rel_resid = 1.0e300;
    double abs_p_error = 1.0e300;
    double dpdrho = 1.0e300;
    bool valid = false;
};

struct DensitySolveDiagnostics {
    std::string phase_label = "density";
    std::string phase_kind = "liq";
    std::string rejection_reason = "";
    double t = 0.0;
    double p = 0.0;
    std::vector<double> composition;
    int scan_point_count = 0;
    int finite_point_count = 0;
    int coarse_bracket_count = 0;
    int refined_bracket_count = 0;
    int candidate_root_count = 0;
    DensityCandidateDiagnostics best_near_root;
    std::vector<DensityCandidateDiagnostics> candidate_roots;
    bool fallback_used = false;
    std::string fallback_rejected_reason = "";
    std::string warm_start_source = "scan";
    std::string validity_gate = "not_evaluated";
};

struct DensitySolveResult {
    double rho = 0.0;
    bool valid = false;
    DensitySolveDiagnostics diagnostics;
};

struct PureNeutralRegressionDensityRecord {
    double t = 0.0;
    double p = 0.0;
    double rho_exp = 0.0;
    int phase = 0;
};

struct PureNeutralRegressionVLERecord {
    double t = 0.0;
    double p = 0.0;
};

struct PureNeutralRegressionDebugResult {
    double objective = 0.0;
    vector<double> gradient;
    vector<double> residuals;
    vector<double> jacobian_row_major;
    int jacobian_rows = 0;
    int jacobian_cols = 0;
    bool jacobian_available = true;
    std::string jacobian_backend = "autodiff";
    bool jacobian_fallback_used = false;
    std::string jacobian_fallback_reason;
    vector<double> hessian_row_major;
    int hessian_rows = 0;
    int hessian_cols = 0;
    bool hessian_available = false;
    std::string hessian_backend = "not_implemented";
    bool hessian_fallback_used = false;
    std::string hessian_fallback_reason = "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
    vector<double> density_raw_residuals;
    vector<double> pure_vle_raw_residuals;
    int residual_evaluations = 0;
    int density_solves = 0;
    int fused_state_evaluations = 0;
    double callback_wall_time_s = 0.0;
};

struct PureNeutralRegressionResult {
    vector<double> x;
    double cost = HUGE_DBL;
    double residual_norm = HUGE_DBL;
    double density_metric = HUGE_DBL;
    double pure_vle_metric = HUGE_DBL;
    double initial_cost = HUGE_DBL;
    double initial_density_metric = HUGE_DBL;
    double initial_pure_vle_metric = HUGE_DBL;
    bool success = false;
    int status = 0;
    int nfev = 0;
    int iterations = 0;
    int starts_tried = 0;
    int objective_evaluations = 0;
    int gradient_evaluations = 0;
    int residual_evaluations = 0;
    int density_solves = 0;
    int fused_state_evaluations = 0;
    double callback_wall_time_s = 0.0;
    double solve_wall_time_s = 0.0;
    std::string message;
    std::string backend;
    std::string optimizer_backend;
    std::string derivative_backend;
    double gradient_norm = 0.0;
    double step_norm = 0.0;
    bool jacobian_available = true;
    std::string jacobian_backend = "autodiff";
    bool jacobian_fallback_used = false;
    std::string jacobian_fallback_reason;
    bool hessian_available = false;
    std::string hessian_backend = "not_implemented";
    bool hessian_fallback_used = false;
    std::string hessian_fallback_reason = "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
};

struct GenericRegressionRecord {
    std::string term_name;
    int term = 0;
    double t = 0.0;
    double p = 0.0;
    int phase = 0;
    vector<double> x;
    vector<double> y;
    double target = 0.0;
    int target_index = -1;
    int target_index_2 = -1;
    int density_kind = 0;
    int activity_basis = 0;
    int solvent_index = -1;
    double scale = 1.0;
};

struct GenericRegressionDebugResult {
    double cost = HUGE_DBL;
    double residual_norm = HUGE_DBL;
    vector<double> residuals;
    vector<double> jacobian_row_major;
    int jacobian_rows = 0;
    int jacobian_cols = 0;
    bool jacobian_available = true;
    std::string jacobian_backend = "not_available";
    bool jacobian_fallback_used = true;
    std::string jacobian_fallback_reason = "Generic regression sensitivities are not implemented for all residual state calls yet.";
    std::string not_available_reason = "not_available";
    vector<double> hessian_row_major;
    int hessian_rows = 0;
    int hessian_cols = 0;
    bool hessian_available = false;
    std::string hessian_backend = "not_implemented";
    bool hessian_fallback_used = false;
    std::string hessian_fallback_reason = "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
    std::map<std::string, double> metrics_by_term;
};

struct GenericRegressionResult {
    vector<double> x;
    double cost = HUGE_DBL;
    double residual_norm = HUGE_DBL;
    double initial_cost = HUGE_DBL;
    double initial_residual_norm = HUGE_DBL;
    std::map<std::string, double> metrics_by_term;
    bool success = false;
    int status = 0;
    int nfev = 0;
    int iterations = 0;
    int starts_tried = 0;
    bool jacobian_available = true;
    std::string jacobian_backend = "not_available";
    bool jacobian_fallback_used = true;
    std::string jacobian_fallback_reason = "Generic regression sensitivities are not implemented for all residual state calls yet.";
    std::string not_available_reason = "not_available";
    bool hessian_available = false;
    std::string hessian_backend = "not_implemented";
    bool hessian_fallback_used = false;
    std::string hessian_fallback_reason = "Hessian support is a skeleton for future IPOPT-compatible optimizer integration.";
    std::string message;
    std::string backend;
};

struct NeutralBinaryKijPhaseDerivatives {
    double pressure = 0.0;
    double rho = 0.0;
    double z = 0.0;
    double dpdrho = 0.0;
    double dpdk = 0.0;
    double drhodk = 0.0;
    vector<double> mu_res;
    vector<double> dmu_res_dk_fixed_rho;
    vector<double> lnphi;
    vector<double> dlnphi_drho;
    vector<double> dlnphi_dk_fixed_rho;
    vector<double> dlnphi_dk_total;
    std::string backend = "cppad_implicit";
};

struct PhaseStateCompositionSensitivityResult {
    bool supported = false;
    std::string backend = "not_available";
    std::string density_backend = "not_available";
    std::string message;
    double temperature = 0.0;
    double pressure = 0.0;
    double density = 0.0;
    double pressure_density_derivative = 0.0;
    int rows = 0;
    int cols = 0;
    vector<double> composition;
    vector<double> ln_fugacity;
    vector<double> density_composition_derivative;
    vector<double> pressure_composition_fixed_density_derivative;
    vector<double> ln_fugacity_density_derivative;
    vector<double> fixed_density_jacobian_row_major;
    vector<double> jacobian_row_major;
};

class ePCSAFTMixtureNative;

class ePCSAFTStateNative {
public:
    ePCSAFTStateNative(std::shared_ptr<ePCSAFTMixtureNative> mixture, double t, vector<double> x,
        int phase, bool has_p, double p, bool has_rho, double rho,
        bool has_rho_guess = false, double rho_guess = 0.0);

    double temperature() const;
    int phase() const;
    const vector<double>& composition() const;

    double pressure();
    double density();
    double compressibility_factor();
    CompressibilityFactorResult compressibility_factor_result();
    double residual_helmholtz();
    ScalarContributionTerms residual_helmholtz_result();
    double temperature_derivative_residual_helmholtz();
    ScalarContributionTerms temperature_derivative_residual_helmholtz_result();
    double residual_enthalpy();
    double residual_entropy();
    double residual_gibbs();
    vector<double> residual_chemical_potential();
    ResidualChemicalPotentialResult residual_chemical_potential_result();
    CompositionContributionResult composition_derivative_residual_helmholtz_result();
    vector<double> ln_fugacity_coefficient();
    vector<double> fugacity_coefficient();
    FugacityContributionResult fugacity_coefficient_result();
    BornSSMDSDerivativeResult born_ssmds_liquid_derivatives();
    vector<double> relative_permittivity();
    double osmotic_coefficient();
    vector<double> solvation_free_energy();
    ActivityCoefficientNative activity_coefficient_native(
        bool include_aux = true,
        bool has_solvent_override = false,
        int solvent_override_index = -1
    );

private:
    std::shared_ptr<ePCSAFTMixtureNative> mixture_;
    double t_;
    vector<double> x_;
    int phase_;
    bool has_p_;
    bool has_rho_;
    double p_;
    double rho_;
    bool pressure_cached_;
    bool density_cached_;
    bool activity_coefficient_cached_;
    ActivityCoefficientNative activity_coefficient_cache_;
};

class ePCSAFTMixtureNative : public std::enable_shared_from_this<ePCSAFTMixtureNative> {
public:
    explicit ePCSAFTMixtureNative(const add_args& args);
    const add_args& args() const;
    std::shared_ptr<ePCSAFTStateNative> state(double t, vector<double> x, int phase,
        bool has_p, double p, bool has_rho, double rho, bool has_rho_guess = false, double rho_guess = 0.0);
    size_t ncomp() const;
    bool has_ionic() const;
    const vector<int>& cation_indices() const;
    const vector<int>& anion_indices() const;
    const vector<int>& solvent_indices() const;
    const vector<int>& pair_cation_indices() const;
    const vector<int>& pair_anion_indices() const;
    const vector<int>& pair_nu_cation() const;
    const vector<int>& pair_nu_anion() const;
    double solve_density(double t, double p, const vector<double>& x, int phase);
    double solve_density_with_guess(double t, double p, const vector<double>& x, int phase, double rho_guess);
    double solve_density_scoped(double t, double p, const vector<double>& x, int phase, const std::string& scope);
    const DensitySolveDiagnostics& last_density_diagnostics() const;
    bool lookup_reference_state(const ReferenceStateKey& key, ReferenceStateValue* out);
    void store_reference_state(const ReferenceStateKey& key, const ReferenceStateValue& value);
    void clear_runtime_caches();
    void reset_runtime_cache_stats();
    size_t reference_state_cache_hits() const;
    size_t reference_state_cache_misses() const;
    size_t density_warm_start_hits() const;
    size_t density_warm_start_fallbacks() const;

private:
    add_args args_;
    bool has_ionic_;
    vector<int> cation_indices_;
    vector<int> anion_indices_;
    vector<int> solvent_indices_;
    vector<int> pair_cation_indices_;
    vector<int> pair_anion_indices_;
    vector<int> pair_nu_cation_;
    vector<int> pair_nu_anion_;
    vector<ReferenceStateCacheEntry> reference_state_cache_;
    double liquid_density_seed_ = 0.0;
    double vapor_density_seed_ = 0.0;
    bool liquid_density_seed_valid_ = false;
    bool vapor_density_seed_valid_ = false;
    std::map<std::string, double> scoped_density_seeds_;
    DensitySolveDiagnostics last_density_diagnostics_;
    size_t reference_state_cache_hits_ = 0;
    size_t reference_state_cache_misses_ = 0;
    size_t density_warm_start_hits_ = 0;
    size_t density_warm_start_fallbacks_ = 0;
};

double Z_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
vector<double> mures_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
vector<double> lnfug_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
vector<double> fugcoef_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double p_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double den_cpp(double t, double p, vector<double> x, int phase, const add_args &cppargs);
double ares_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
NeutralBinaryKijPhaseDerivatives neutral_binary_kij_phase_derivatives_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs,
    int k_index
);
NeutralBinaryKijPhaseDerivatives neutral_binary_pair_parameter_phase_derivatives_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs,
    int parameter_index,
    const std::string &parameter_name
);
NeutralBinaryKijPhaseDerivatives generic_component_parameter_phase_derivatives_cpp(
    double t,
    double rho,
    const vector<double> &x,
    const add_args &cppargs,
    int target_kind,
    int target_index
);
double dadt_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double hres_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double sres_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double gres_cpp(double t, double rho, vector<double> x, const add_args &cppargs);

// functions used to solve for XA and its derivatives
vector<double> XA_find(vector<double> XA_guess, vector<double> delta_ij, double den,
    vector<double> x);
vector<double> dXAdx_find(vector<int> assoc_num, vector<double> delta_ij,
    double den, vector<double> XA, vector<double> ddelta_dx, vector<double> x);
vector<double> dXAdt_find(vector<double> delta_ij, double den,
    vector<double> XA, vector<double> ddelta_dt, vector<double> x);

double dielectric_eps_cpp(vector<double> x, const add_args &cppargs);
vector<double> dielectric_diff_cpp(vector<double> x, const add_args &cppargs);
PureNeutralRegressionResult fit_pure_neutral_ceres_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart
);
PureNeutralRegressionDebugResult evaluate_pure_neutral_objective_debug_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x
);
GenericRegressionDebugResult evaluate_generic_regression_debug_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &x
);
GenericRegressionResult fit_generic_least_squares_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart,
    int max_nfev
);
GenericRegressionResult fit_generic_ceres_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart,
    int max_nfev
);

class ValueError: public std::exception
{
public:
    ValueError(const std::string &err) throw() : m_err(err) {}
    ~ValueError() throw() {};
    virtual const char* what() const throw() { return m_err.c_str(); }
private:
    std::string m_err;
};

class SolutionError: public std::exception
{
public:
    SolutionError(const std::string &err) throw() : m_err(err) {}
    ~SolutionError() throw() {};
    virtual const char* what() const throw() { return m_err.c_str(); }
private:
    std::string m_err;
};

// density and composition helpers
double density_root_residual_cpp(double rhomolar, double t, double p, vector<double> x, const add_args &cppargs);
double density_brent_cpp(double t, double p, vector<double> x, const add_args &cppargs, double a, double b,
    double macheps, double tol_abs, int maxiter);
double reduced_density_to_molar(double nu, double t, int ncomp, vector<double> x, const add_args &cppargs);
vector<double> association_site_fractions_cpp(vector<double> XA_guess, vector<double> delta_ij, double den,
    vector<double> x);
vector<double> association_site_fraction_dt_cpp(vector<double> delta_ij, double den,
    vector<double> XA, vector<double> ddelta_dt, vector<double> x);
vector<double> association_site_fraction_dx_cpp(vector<int> assoc_num, vector<double> delta_ij,
    double den, vector<double> XA, vector<double> ddelta_dx, vector<double> x);


