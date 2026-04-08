#include <vector>
#include <string>
#include <memory>

using std::vector;

const static double kb = 1.380648465952442093e-23; // Boltzmann constant, J K^-1
const static double PI = 3.141592653589793;
const static double N_AV = 6.022140857e23; // Avagadro's number
const static double E_CHRG = 1.6021766208e-19; // elementary charge, units of coulomb
const static double perm_vac = 8.854187817e-22; //permittivity in vacuum, C V^-1 Angstrom^-1

#ifndef DBL_EPSILON
    #include <limits>
    #define DBL_EPSILON std::numeric_limits<double>::epsilon()
#endif

struct add_args {
    vector<double> m;
    vector<double> s;
    vector<double> e;
    vector<double> k_ij;
    vector<double> e_assoc;
    vector<double> vol_a;
    vector<double> dipm;
    vector<double> dip_num;
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
    int polar_dadx_diff_mode;
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
    int debug;
    vector<int> assoc_num;
    vector<int> assoc_matrix;
    vector<double> k_hb;
    vector<double> l_ij;
};

struct FlashResultNative {
    double value;
    vector<double> xl;
    vector<double> xv;
};

struct VaporizationResultNative {
    double value;
    double pressure;
};

struct ActivityCoeffNative {
    vector<double> gamma_components;
    vector<double> gamma_mean_ionic_x;
    vector<double> gamma_mean_ionic_m;
    vector<double> gsolv;
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
    double osmotic_c;
};

class PCSAFTMixtureNative;

class PCSAFTStateNative {
public:
    PCSAFTStateNative(std::shared_ptr<PCSAFTMixtureNative> mixture, double t, vector<double> x,
        int phase, bool has_p, double p, bool has_rho, double rho);

    double temperature() const;
    int phase() const;
    const vector<double>& composition() const;

    double pressure();
    double density();
    double Z();
    double a_res();
    double ares();
    double dadt();
    double h_res();
    double hres();
    double s_res();
    double sres();
    double g_res();
    double gres();
    vector<double> mu_res();
    vector<double> gamma();
    vector<double> lnfugcoef();
    vector<double> fugcoef();
    vector<double> lnfugcoef_terms();
    vector<double> dielectric_eval();
    vector<double> dielc_eval();
    double osmoticC();
    vector<double> miac_m();
    vector<double> miac();
    vector<double> gsolv();
    ActivityCoeffNative actcoeff(bool has_solvent_override = false, int solvent_override_index = -1);
    FlashResultNative flashTQ(double q, bool has_p_guess = false, double p_guess = 0.0);
    FlashResultNative flashPQ(double p, double q, bool has_t_guess = false, double t_guess = 0.0);
    VaporizationResultNative Hvap(bool has_p_guess = false, double p_guess = 0.0);

private:
    std::shared_ptr<PCSAFTMixtureNative> mixture_;
    double t_;
    vector<double> x_;
    int phase_;
    bool has_p_;
    bool has_rho_;
    double p_;
    double rho_;
    bool pressure_cached_;
    bool density_cached_;
    bool actcoeff_cached_;
    ActivityCoeffNative actcoeff_cache_;
};

class PCSAFTMixtureNative : public std::enable_shared_from_this<PCSAFTMixtureNative> {
public:
    explicit PCSAFTMixtureNative(const add_args& args);
    const add_args& args() const;
    std::shared_ptr<PCSAFTStateNative> state(double t, vector<double> x, int phase,
        bool has_p, double p, bool has_rho, double rho);
    size_t ncomp() const;
    bool has_ionic() const;
    const vector<int>& cation_indices() const;
    const vector<int>& anion_indices() const;
    const vector<int>& solvent_indices() const;
    const vector<int>& pair_cation_indices() const;
    const vector<int>& pair_anion_indices() const;
    const vector<int>& pair_nu_cation() const;
    const vector<int>& pair_nu_anion() const;

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
};

double Z_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
vector<double> mures_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
vector<double> lnfug_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
vector<double> lnfug_terms_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
vector<double> fugcoef_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double p_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double den_cpp(double t, double p, vector<double> x, int phase, const add_args &cppargs);
double ares_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double dadt_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double hres_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double sres_cpp(double t, double rho, vector<double> x, const add_args &cppargs);
double gres_cpp(double t, double rho, vector<double> x, const add_args &cppargs);

vector<double> flashTQ_cpp(double t, double Q, vector<double> x, const add_args &cppargs);
vector<double> flashTQ_cpp(double t, double Q, vector<double> x, const add_args &cppargs, double p_guess); // used if a guess value is given
vector<double> flashPQ_cpp(double p, double Q, vector<double> x, const add_args &cppargs);
vector<double> flashPQ_cpp(double p, double Q, vector<double> x, const add_args &cppargs, double t_guess); // used if a guess value is given

// functions used to solve for XA and its derivatives
vector<double> XA_find(vector<double> XA_guess, vector<double> delta_ij, double den,
    vector<double> x);
vector<double> dXAdx_find(vector<int> assoc_num, vector<double> delta_ij,
    double den, vector<double> XA, vector<double> ddelta_dx, vector<double> x);
vector<double> dXAdt_find(vector<double> delta_ij, double den,
    vector<double> XA, vector<double> ddelta_dt, vector<double> x);

// helper functions
inline bool IsNotZero (double x) {return x != 0.0;}
double reduced_to_molar(double nu, double t, int ncomp, vector<double> x, const add_args &cppargs);
double calc_water_sigma(double t);
inline double calc_sigma(double t, double (*function)(double)){return function(t);} // this can allow us to accept a custom function for a temperature dependent sigma
 add_args single_component_args_cpp(int i, const add_args &cppargs);
double dielectric_eps_cpp(vector<double> x, const add_args &cppargs);
vector<double> dielectric_diff_cpp(vector<double> x, const add_args &cppargs);
double dielc_eps_cpp(vector<double> x, const add_args &cppargs);
vector<double> dielc_diff_cpp(vector<double> x, const add_args &cppargs);

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

// functions used in flash calculations and root finding
vector<double> flash_pq_cpp(double t_guess, double p, double Q, vector<double> x, add_args cppargs);
vector<double> flash_tq_cpp(double p_guess, double t, double Q, vector<double> x, add_args cppargs);
double density_root_residual_cpp(double rhomolar, double t, double p, vector<double> x, const add_args &cppargs);
double density_brent_cpp(double t, double p, vector<double> x, int phase, const add_args &cppargs, double a, double b,
    double macheps, double tol_abs, int maxiter);
namespace thermo_detail {
double vle_temperature_guess_cpp(double p, double Q, vector<double> x, add_args cppargs);
double vle_pressure_guess_cpp(double t, double Q, vector<double> x, const add_args &cppargs);
}
double secant_residual_cpp(double R, double kb0, double Q, vector<double> u, vector<double> x, const add_args &cppargs);
double bounded_secant_cpp(double kb0, double Q, vector<double> u, vector<double> x, const add_args &cppargs, double x0, double xmin,
    double xmax, double dx, double tol, int maxiter);
double reduced_density_to_molar(double nu, double t, int ncomp, vector<double> x, const add_args &cppargs);
vector<double> association_site_fractions_cpp(vector<double> XA_guess, vector<double> delta_ij, double den,
    vector<double> x);
vector<double> association_site_fraction_dt_cpp(vector<double> delta_ij, double den,
    vector<double> XA, vector<double> ddelta_dt, vector<double> x);
vector<double> association_site_fraction_dx_cpp(vector<int> assoc_num, vector<double> delta_ij,
    double den, vector<double> XA, vector<double> ddelta_dx, vector<double> x);

