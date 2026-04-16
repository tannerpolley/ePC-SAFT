#pragma once

#include "epcsaft_core_internal.h"

namespace thermo_detail {

struct HardChainState {
    vector<double> zeta;
    vector<double> ghs;
    double eta = 0.0;
};

struct AutodiffHardChainState {
    vector<AutoDual> zeta;
    vector<AutoDual> ghs;
    AutoDual eta = 0.0;
};

struct ContributionDadxResult {
    vector<double> dadx;
    double ares = 0.0;
    double z = 0.0;
    double sum_x_dadx = 0.0;
};

struct AssociationIntermediateState {
    bool active = false;
    AssociationSetup setup;
    vector<double> XA;
    vector<double> dXA_dt;
    vector<double> dXA_dx;
};

struct IonIntermediateState {
    bool active = false;
    DielectricState dielectric;
    double charge_square_sum = 0.0;
    double kappa = 0.0;
    vector<double> chi;
    vector<double> sigma_k;
    double chi_sum = 0.0;
    double sigma_sum = 0.0;
    vector<double> dkappa_dx;
    vector<double> dchi_sum_dx;
};

struct BornIntermediateState {
    int model = 0;
    double eps_value = 0.0;
    vector<double> deps_dx;
    double charge_radius_sum = 0.0;
    double charge_radius_sum_dt = 0.0;
    BornSSMDSData shell;
};

}  // namespace thermo_detail

using thermo_detail::HardChainState;
using thermo_detail::AutodiffHardChainState;
using thermo_detail::AssociationIntermediateState;
using thermo_detail::BornIntermediateState;
using thermo_detail::ContributionDadxResult;
using thermo_detail::IonIntermediateState;

HardChainState hard_chain_state_cpp(const MixtureState &thermo, const vector<double> &x, const add_args &cppargs);
AutodiffHardChainState hard_chain_state_autodiff_cpp(double den, const vector<double> &d, const vector<AutoDual> &x, const add_args &cppargs);
double pair_diameter_cpp(double d_i, double d_j);
double hs_contact_value_cpp(double pair_diameter, double zeta2, double zeta3);
double association_volume_cpp(int comp_i, int comp_j, int ncomp, const vector<double> &s_ij, const add_args &cppargs);

AssociationIntermediateState association_intermediate_state_cpp(
    const MixtureState &thermo,
    const HardChainState &hc_state,
    double t,
    const vector<double> &x,
    const add_args &cppargs,
    bool include_dt,
    bool include_dx,
    const vector<double> *dghs_dt = nullptr
);

IonIntermediateState ion_intermediate_state_cpp(
    const MixtureState &thermo,
    double t,
    const vector<double> &x,
    const add_args &cppargs,
    bool include_dx
);

BornIntermediateState born_intermediate_state_cpp(
    double t,
    const vector<double> &x,
    const add_args &cppargs,
    bool include_dt,
    bool include_dx
);
