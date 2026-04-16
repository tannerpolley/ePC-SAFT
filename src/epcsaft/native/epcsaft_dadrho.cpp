#include "epcsaft_core_internal.h"
#include "contributions/epcsaft_contrib_internal.h"

using namespace thermo_detail;

namespace {

// EqID: hs_ares_dadrho
double dadrho_hs_cpp(const HardChainState &hc_state) {
    const auto &zeta = hc_state.zeta;
    return zeta[3] / (1.0 - zeta[3])
        + 3.0 * zeta[1] * zeta[2] / zeta[0] / std::pow(1.0 - zeta[3], 2.0)
        + (3.0 * std::pow(zeta[2], 3.0) - zeta[3] * std::pow(zeta[2], 3.0)) / zeta[0] / std::pow(1.0 - zeta[3], 3.0);
}

}  // namespace

// EqID: ghs_contact_dadrho
double hs_contact_density_derivative_cpp(double pair_diameter, double zeta2, double zeta3) {
    return zeta3 / std::pow(1.0 - zeta3, 2.0)
        + pair_diameter * (3.0 * zeta2 / std::pow(1.0 - zeta3, 2.0) + 6.0 * zeta2 * zeta3 / std::pow(1.0 - zeta3, 3.0))
        + std::pow(pair_diameter, 2.0) * (4.0 * zeta2 * zeta2 / std::pow(1.0 - zeta3, 3.0) + 6.0 * zeta2 * zeta2 * zeta3 / std::pow(1.0 - zeta3, 4.0));
}

namespace {

// EqID: dadrho_hc
// EqID: hc_ares_dadrho
double dadrho_hc_cpp(const MixtureState &thermo, const HardChainState &hc_state, const vector<double> &x, const add_args &cppargs) {
    int ncomp = static_cast<int>(x.size());
    double summ = 0.0;
    for (int i = 0; i < ncomp; ++i) {
        double pair_diameter = pair_diameter_cpp(thermo.d[i], thermo.d[i]);
        double dghs_drho = hs_contact_density_derivative_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3]);
        summ += x[i] * (cppargs.m[i] - 1.0) / hc_state.ghs[i * ncomp + i] * dghs_drho;
    }
    return thermo.m_avg * dadrho_hs_cpp(hc_state) - summ;
}

// EqID: dadrho_disp
// EqID: disp_ares_dadrho
double dadrho_disp_cpp(const MixtureState &thermo, const HardChainState &hc_state, const DispersionPolynomialState &dispersion) {
    return -2.0 * PI * thermo.den * dispersion.dEtaI1_deta * thermo.m2es3
        - PI * thermo.den * thermo.m_avg * (dispersion.C1 * dispersion.dEtaI2_deta + dispersion.C2 * hc_state.eta * dispersion.I2) * thermo.m2e2s3;
}

vector<double> association_site_fraction_density_terms_cpp(
    const vector<double> &delta_ij,
    double den,
    const vector<double> &XA,
    const vector<double> &ddelta_weighted,
    const vector<double> &x_assoc
) {
    int num_sites = static_cast<int>(XA.size());
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(num_sites, 1);
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(num_sites, num_sites);

    int ij = 0;
    for (int i = 0; i < num_sites; ++i) {
        double summ = 0.0;
        for (int j = 0; j < num_sites; ++j) {
            B(i) -= x_assoc[j] * XA[j] * ddelta_weighted[ij];
            A(i, j) = x_assoc[j] * delta_ij[ij];
            summ += x_assoc[j] * XA[j] * delta_ij[ij];
            ++ij;
        }
        B(i) -= summ;
        A(i, i) = std::pow(1.0 + den * summ, 2.0) / den;
    }

    Eigen::MatrixXd solution = A.lu().solve(B);
    vector<double> dXA_weighted(num_sites, 0.0);
    for (int i = 0; i < num_sites; ++i) {
        dXA_weighted[i] = solution(i);
    }
    return dXA_weighted;
}

// EqID: dadrho_assoc
// EqID: assoc_ares_dadrho
double dadrho_assoc_cpp(
    const MixtureState &thermo,
    const HardChainState &hc_state,
    const AssociationIntermediateState &assoc_state,
    const vector<double> &x,
    const add_args &cppargs,
    double t
) {
    if (!assoc_state.active) {
        return 0.0;
    }

    const vector<int> &site_component_index = assoc_state.setup.site_component_index;
    const vector<double> &x_assoc = assoc_state.setup.x_assoc;
    const vector<double> &delta_ij = assoc_state.setup.delta_ij;
    int ncomp = static_cast<int>(x.size());
    int num_sites = static_cast<int>(site_component_index.size());
    vector<double> ddelta_weighted(num_sites * num_sites, 0.0);

    int ij = 0;
    for (int i = 0; i < num_sites; ++i) {
        int comp_i = site_component_index[i];
        for (int j = 0; j < num_sites; ++j) {
            int comp_j = site_component_index[j];
            if (cppargs.assoc_matrix[ij] != 0) {
                double pair_diameter = pair_diameter_cpp(thermo.d[comp_i], thermo.d[comp_j]);
                double eABij = 0.5 * (cppargs.e_assoc[comp_i] + cppargs.e_assoc[comp_j]);
                double volABij = association_volume_cpp(comp_i, comp_j, ncomp, thermo.s_ij, cppargs);
                ddelta_weighted[ij] = hs_contact_density_derivative_cpp(pair_diameter, hc_state.zeta[2], hc_state.zeta[3])
                    * (std::exp(eABij / t) - 1.0)
                    * std::pow(thermo.s_ij[comp_i * ncomp + comp_j], 3.0)
                    * volABij;
            }
            ++ij;
        }
    }

    vector<double> dXA_weighted = association_site_fraction_density_terms_cpp(delta_ij, thermo.den, assoc_state.XA, ddelta_weighted, x_assoc);

    double value = 0.0;
    for (int i = 0; i < num_sites; ++i) {
        int component_index = site_component_index[i];
        value += x[component_index] * (1.0 / assoc_state.XA[i] - 0.5) * dXA_weighted[i];
    }
    return value;
}

// EqID: dadrho_dh
// EqID: dadrho_dh_explicit
double dadrho_ion_cpp(double t, const IonIntermediateState &ion_state) {
    if (!ion_state.active) {
        return 0.0;
    }
    return -ion_state.kappa / 24.0 / PI / kb / t / (ion_state.dielectric.eps * perm_vac) * ion_state.sigma_sum * E_CHRG * E_CHRG;
}

// EqID: dadrho_born
// EqID: born_ares_dadrho
double dadrho_born_cpp() {
    return 0.0;
}

}  // namespace

DadrhoResult dadrho_result_cpp(double t, double rho, vector<double> x, const add_args &cppargs) {
    MixtureState thermo = mixture_state_cpp(t, rho, x, cppargs, false);
    HardChainState hc_state = hard_chain_state_cpp(thermo, x, cppargs);
    DispersionPolynomialState dispersion = dispersion_polynomials_cpp(thermo.m_avg, hc_state.eta);
    AssociationIntermediateState assoc_state = association_intermediate_state_cpp(thermo, hc_state, t, x, cppargs, false, false);
    IonIntermediateState ion_state = ion_intermediate_state_cpp(thermo, t, x, cppargs, false);

    double hc = dadrho_hc_cpp(thermo, hc_state, x, cppargs);
    double disp = dadrho_disp_cpp(thermo, hc_state, dispersion);
    double assoc = dadrho_assoc_cpp(thermo, hc_state, assoc_state, x, cppargs, t);
    double ion = dadrho_ion_cpp(t, ion_state);
    double born = dadrho_born_cpp();
    double total = hc + disp + assoc + ion + born;

    ScalarContributionTerms raw_terms = make_scalar_terms(hc, disp, assoc, ion, born, total);
    DadrhoResult result;
    result.raw = raw_terms;
    result.terms = normalized_dadrho_terms_cpp(raw_terms);
    return result;
}
