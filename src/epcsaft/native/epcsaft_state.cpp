#include "epcsaft_core_internal.h"

using namespace thermo_detail;

namespace {

int gcd_int(int a, int b) {
    a = std::abs(a);
    b = std::abs(b);
    while (b != 0) {
        int t = a % b;
        a = b;
        b = t;
    }
    return a == 0 ? 1 : a;
}

void build_charge_metadata_cpp(
    const add_args& args,
    bool& has_ionic,
    vector<int>& cation_indices,
    vector<int>& anion_indices,
    vector<int>& solvent_indices,
    vector<int>& pair_cation_indices,
    vector<int>& pair_anion_indices,
    vector<int>& pair_nu_cation,
    vector<int>& pair_nu_anion
) {
    has_ionic = false;
    cation_indices.clear();
    anion_indices.clear();
    solvent_indices.clear();
    pair_cation_indices.clear();
    pair_anion_indices.clear();
    pair_nu_cation.clear();
    pair_nu_anion.clear();

    if (args.z.empty()) {
        return;
    }
    ChargeGroups groups = collect_charge_groups(args, args.z.size());
    cation_indices = groups.cations;
    anion_indices = groups.anions;
    solvent_indices = groups.solvents;
    has_ionic = (!cation_indices.empty() || !anion_indices.empty());

    for (int ic : cation_indices) {
        for (int ia : anion_indices) {
            int zc = static_cast<int>(std::round(std::abs(args.z[ic])));
            int za = static_cast<int>(std::round(std::abs(args.z[ia])));
            int g = gcd_int(zc, za);
            pair_cation_indices.push_back(ic);
            pair_anion_indices.push_back(ia);
            pair_nu_cation.push_back(za / g);
            pair_nu_anion.push_back(zc / g);
        }
    }
}

}  // namespace

ChargeGroups collect_charge_groups(const add_args& args, size_t ncomp) {
    ChargeGroups groups;
    groups.cations.reserve(ncomp);
    groups.anions.reserve(ncomp);
    groups.solvents.reserve(ncomp);
    for (size_t i = 0; i < ncomp; ++i) {
        if (i >= args.z.size()) {
            throw ValueError("Composition and charge vectors must be aligned.");
        }
        if (std::abs(args.z[i]) < 1e-12) {
            groups.solvents.push_back(static_cast<int>(i));
        }
        else if (args.z[i] > 0.0) {
            groups.cations.push_back(static_cast<int>(i));
        }
        else {
            groups.anions.push_back(static_cast<int>(i));
        }
    }
    return groups;
}

ePCSAFTMixtureNative::ePCSAFTMixtureNative(const add_args& args)
    : args_(args), has_ionic_(false)
{
    build_charge_metadata_cpp(
        args_,
        has_ionic_,
        cation_indices_,
        anion_indices_,
        solvent_indices_,
        pair_cation_indices_,
        pair_anion_indices_,
        pair_nu_cation_,
        pair_nu_anion_
    );
}

const add_args& ePCSAFTMixtureNative::args() const
{
    return args_;
}

size_t ePCSAFTMixtureNative::ncomp() const
{
    return args_.m.size();
}

bool ePCSAFTMixtureNative::has_ionic() const
{
    return has_ionic_;
}

const vector<int>& ePCSAFTMixtureNative::cation_indices() const
{
    return cation_indices_;
}

const vector<int>& ePCSAFTMixtureNative::anion_indices() const
{
    return anion_indices_;
}

const vector<int>& ePCSAFTMixtureNative::solvent_indices() const
{
    return solvent_indices_;
}

const vector<int>& ePCSAFTMixtureNative::pair_cation_indices() const
{
    return pair_cation_indices_;
}

const vector<int>& ePCSAFTMixtureNative::pair_anion_indices() const
{
    return pair_anion_indices_;
}

const vector<int>& ePCSAFTMixtureNative::pair_nu_cation() const
{
    return pair_nu_cation_;
}

const vector<int>& ePCSAFTMixtureNative::pair_nu_anion() const
{
    return pair_nu_anion_;
}

std::shared_ptr<ePCSAFTStateNative> ePCSAFTMixtureNative::state(double t, vector<double> x, int phase,
    bool has_p, double p, bool has_rho, double rho)
{
    return std::make_shared<ePCSAFTStateNative>(shared_from_this(), t, std::move(x), phase, has_p, p, has_rho, rho);
}

ePCSAFTStateNative::ePCSAFTStateNative(std::shared_ptr<ePCSAFTMixtureNative> mixture, double t, vector<double> x,
    int phase, bool has_p, double p, bool has_rho, double rho)
    : mixture_(std::move(mixture)), t_(t), x_(std::move(x)), phase_(phase),
      has_p_(has_p), has_rho_(has_rho), p_(p), rho_(rho),
      pressure_cached_(has_p), density_cached_(has_rho), activity_coefficient_cached_(false)
{
    if (!mixture_) {
        throw ValueError("ePCSAFTStateNative requires a valid mixture.");
    }
    if (x_.size() != mixture_->ncomp()) {
        throw ValueError("State composition size does not match mixture size.");
    }
    if (phase_ != 0 && phase_ != 1) {
        throw ValueError("phase must be 0 (liquid) or 1 (vapor).");
    }
    if (pressure_cached_ && !density_cached_) {
        const add_args& args = mixture_->args();
        rho_ = den_cpp(t_, p_, x_, phase_, args);
        has_rho_ = true;
        density_cached_ = true;
    }
}

double ePCSAFTStateNative::temperature() const
{
    return t_;
}

int ePCSAFTStateNative::phase() const
{
    return phase_;
}

const vector<double>& ePCSAFTStateNative::composition() const
{
    return x_;
}

double ePCSAFTStateNative::pressure()
{
    if (pressure_cached_) {
        return p_;
    }
    if (!density_cached_) {
        throw ValueError("ePCSAFTStateNative cannot compute pressure without density or pressure data.");
    }
    const add_args& args = mixture_->args();
    p_ = p_cpp(t_, rho_, x_, args);
    pressure_cached_ = true;
    return p_;
}

double ePCSAFTStateNative::density()
{
    if (density_cached_) {
        return rho_;
    }
    if (!pressure_cached_) {
        throw ValueError("ePCSAFTStateNative cannot compute density without pressure or density data.");
    }
    const add_args& args = mixture_->args();
    rho_ = den_cpp(t_, p_, x_, phase_, args);
    density_cached_ = true;
    return rho_;
}

double ePCSAFTStateNative::compressibility_factor()
{
    const add_args& args = mixture_->args();
    return Z_cpp(t_, density(), x_, args);
}

CompressibilityFactorResult ePCSAFTStateNative::compressibility_factor_result()
{
    const add_args& args = mixture_->args();
    return compressibility_factor_result_cpp(t_, density(), x_, args);
}

double ePCSAFTStateNative::residual_helmholtz()
{
    const add_args& args = mixture_->args();
    return ares_cpp(t_, density(), x_, args);
}

ScalarContributionTerms ePCSAFTStateNative::residual_helmholtz_result()
{
    const add_args& args = mixture_->args();
    return residual_helmholtz_result_cpp(t_, density(), x_, args);
}

double ePCSAFTStateNative::temperature_derivative_residual_helmholtz()
{
    const add_args& args = mixture_->args();
    return dadt_cpp(t_, density(), x_, args);
}

ScalarContributionTerms ePCSAFTStateNative::temperature_derivative_residual_helmholtz_result()
{
    const add_args& args = mixture_->args();
    return temperature_derivative_residual_helmholtz_result_cpp(t_, density(), x_, args);
}

double ePCSAFTStateNative::residual_enthalpy()
{
    const add_args& args = mixture_->args();
    return hres_cpp(t_, density(), x_, args);
}

double ePCSAFTStateNative::residual_entropy()
{
    const add_args& args = mixture_->args();
    return sres_cpp(t_, density(), x_, args);
}

double ePCSAFTStateNative::residual_gibbs()
{
    const add_args& args = mixture_->args();
    return gres_cpp(t_, density(), x_, args);
}

vector<double> ePCSAFTStateNative::residual_chemical_potential()
{
    const add_args& args = mixture_->args();
    return mures_cpp(t_, density(), x_, args);
}

ResidualChemicalPotentialResult ePCSAFTStateNative::residual_chemical_potential_result()
{
    const add_args& args = mixture_->args();
    return residual_chemical_potential_result_cpp(t_, density(), x_, args);
}

CompositionContributionResult ePCSAFTStateNative::composition_derivative_residual_helmholtz_result()
{
    const add_args& args = mixture_->args();
    return composition_derivative_residual_helmholtz_result_cpp(t_, density(), x_, args);
}

vector<double> ePCSAFTStateNative::ln_fugacity_coefficient()
{
    const add_args& args = mixture_->args();
    return lnfug_cpp(t_, density(), x_, args);
}

vector<double> ePCSAFTStateNative::fugacity_coefficient()
{
    const add_args& args = mixture_->args();
    return fugcoef_cpp(t_, density(), x_, args);
}

FugacityContributionResult ePCSAFTStateNative::fugacity_coefficient_result()
{
    const add_args& args = mixture_->args();
    return fugacity_coefficient_result_cpp(t_, density(), x_, args);
}

vector<double> ePCSAFTStateNative::relative_permittivity()
{
    const add_args& args = mixture_->args();
    vector<double> out;
    out.push_back(dielectric_eps_cpp(x_, args));
    vector<double> deps = dielectric_diff_cpp(x_, args);
    out.insert(out.end(), deps.begin(), deps.end());
    return out;
}

double dielectric_eps_cpp(vector<double> x, const add_args &cppargs) {
    return dielectric_state_cpp(x, cppargs).eps;
}

vector<double> dielectric_diff_cpp(vector<double> x, const add_args &cppargs) {
    return dielectric_state_cpp(x, cppargs).deps_dx;
}

double dielc_eps_cpp(vector<double> x, const add_args &cppargs) {
    return dielectric_eps_cpp(std::move(x), cppargs);
}

vector<double> dielc_diff_cpp(vector<double> x, const add_args &cppargs) {
    return dielectric_diff_cpp(std::move(x), cppargs);
}

double ePCSAFTStateNative::osmotic_coefficient()
{
    return activity_coefficient_native(false, -1).osmotic_coefficient;
}

vector<double> ePCSAFTStateNative::solvation_free_energy()
{
    return activity_coefficient_native(false, -1).solvation_free_energy;
}

ActivityCoefficientNative ePCSAFTStateNative::activity_coefficient_native(bool has_solvent_override, int solvent_override_index)
{
    if (!mixture_->has_ionic()) {
        throw ValueError("activity_coefficient requires ionic species (non-zero z).");
    }
    if (!has_solvent_override && activity_coefficient_cached_) {
        return activity_coefficient_cache_;
    }
    const add_args& args = mixture_->args();
    double rho = density();
    double p = pressure();
    ActivityCoefficientNative out = activity_coefficient_values_cpp(
        t_, rho, p, phase_, x_, args,
        mixture_->cation_indices(),
        mixture_->anion_indices(),
        mixture_->solvent_indices(),
        mixture_->pair_cation_indices(),
        mixture_->pair_anion_indices(),
        mixture_->pair_nu_cation(),
        mixture_->pair_nu_anion(),
        true,
        has_solvent_override,
        solvent_override_index
    );
    if (!has_solvent_override) {
        activity_coefficient_cache_ = out;
        activity_coefficient_cached_ = true;
    }
    return out;
}
