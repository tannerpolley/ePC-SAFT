#pragma once

#include <cmath>
#include <string>
#include <vector>

#ifdef EPCSAFT_HAS_CPPAD
#include <cppad/cppad.hpp>
#endif

namespace epcsaft::autodiff {

#ifdef EPCSAFT_HAS_CPPAD
using CppADScalar = CppAD::AD<double>;
#endif

inline bool cppad_compiled() {
#ifdef EPCSAFT_HAS_CPPAD
    return true;
#else
    return false;
#endif
}

template <class Scalar>
inline Scalar scaled_residual(const Scalar& predicted, const Scalar& observed, const Scalar& scale) {
    return (predicted - observed) * scale;
}

template <class Scalar>
inline Scalar pressure_log_residual(const Scalar& predicted_pressure, const Scalar& observed_pressure) {
    using std::log;
    return log(predicted_pressure) - log(observed_pressure);
}

template <class Scalar>
inline Scalar reaction_log_residual(const Scalar& ln_q, const Scalar& ln_k) {
    return ln_q - ln_k;
}

}  // namespace epcsaft::autodiff
