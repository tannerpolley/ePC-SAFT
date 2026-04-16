#include "epcsaft_contrib_internal.h"

using namespace thermo_detail;

namespace {

template <typename Scalar>
Scalar hs_contact_value_scalar_cpp(double pair_diameter, const Scalar &zeta2, const Scalar &zeta3) {
    const Scalar one = scalar_constant<Scalar>(1.0);
    return one / (one - zeta3)
        + pair_diameter * 3.0 * zeta2 / scalar_pow(one - zeta3, 2)
        + scalar_pow(pair_diameter, 2.0) * 2.0 * zeta2 * zeta2 / scalar_pow(one - zeta3, 3);
}

template <typename Scalar>
struct HardChainStateScalar {
    vector<Scalar> zeta;
    vector<Scalar> ghs;
    Scalar eta = scalar_constant<Scalar>(0.0);
};

template <typename Scalar>
HardChainStateScalar<Scalar> hard_chain_state_scalar_cpp(double den, const vector<double> &d, const vector<Scalar> &x, const add_args &cppargs) {
    HardChainStateScalar<Scalar> state;
    int ncomp = static_cast<int>(x.size());
    state.zeta.assign(4, scalar_constant<Scalar>(0.0));
    state.ghs.assign(ncomp * ncomp, scalar_constant<Scalar>(0.0));

    for (int k = 0; k < 4; ++k) {
        Scalar summ = scalar_constant<Scalar>(0.0);
        for (int j = 0; j < ncomp; ++j) {
            summ += x[j] * cppargs.m[j] * scalar_pow(d[j], k);
        }
        state.zeta[k] = PI / 6.0 * den * summ;
    }

    state.eta = state.zeta[3];

    int idx = -1;
    for (int i = 0; i < ncomp; ++i) {
        for (int j = 0; j < ncomp; ++j) {
            ++idx;
            double pair_diameter = pair_diameter_cpp(d[i], d[j]);
            state.ghs[idx] = hs_contact_value_scalar_cpp(pair_diameter, state.zeta[2], state.zeta[3]);
        }
    }

    return state;
}

}  // namespace

double pair_diameter_cpp(double d_i, double d_j) {
    return d_i * d_j / (d_i + d_j);
}

double hs_contact_value_cpp(double pair_diameter, double zeta2, double zeta3) {
    return hs_contact_value_scalar_cpp(pair_diameter, zeta2, zeta3);
}

HardChainState hard_chain_state_cpp(const MixtureState &thermo, const vector<double> &x, const add_args &cppargs) {
    auto scalar_state = hard_chain_state_scalar_cpp(thermo.den, thermo.d, x, cppargs);
    HardChainState state;
    state.zeta = std::move(scalar_state.zeta);
    state.ghs = std::move(scalar_state.ghs);
    state.eta = scalar_state.eta;
    return state;
}

AutodiffHardChainState hard_chain_state_autodiff_cpp(double den, const vector<double> &d, const vector<AutoDual> &x, const add_args &cppargs) {
    auto scalar_state = hard_chain_state_scalar_cpp(den, d, x, cppargs);
    AutodiffHardChainState state;
    state.zeta = std::move(scalar_state.zeta);
    state.ghs = std::move(scalar_state.ghs);
    state.eta = scalar_state.eta;
    return state;
}
