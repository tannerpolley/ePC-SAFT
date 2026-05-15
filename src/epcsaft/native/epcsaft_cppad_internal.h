#pragma once

#include <cmath>
#include <string>

#include "cppad/cppad_scalar.h"

#ifdef EPCSAFT_HAS_CPPAD
using epcsaft::native::cppad_support::CppADScalar;
using epcsaft::native::cppad_support::make_cppad_scalar;
#endif

template <typename Scalar>
inline Scalar scalar_constant(double value) {
    return static_cast<Scalar>(value);
}

#ifdef EPCSAFT_HAS_CPPAD
template <>
inline CppADScalar scalar_constant<CppADScalar>(double value) {
    return make_cppad_scalar(value);
}
#endif

inline double scalar_value(double x) {
    return x;
}

#ifdef EPCSAFT_HAS_CPPAD
inline double scalar_value(const CppADScalar &x) {
    return CppAD::Value(x);
}
#endif

inline double scalar_derivative(double) {
    return 0.0;
}

inline double scalar_log(double x) {
    return std::log(x);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppADScalar scalar_log(const CppADScalar &x) {
    return CppAD::log(x);
}
#endif

inline double scalar_exp(double x) {
    return std::exp(x);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppADScalar scalar_exp(const CppADScalar &x) {
    return CppAD::exp(x);
}
#endif

inline double scalar_sqrt(double x) {
    return std::sqrt(x);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppADScalar scalar_sqrt(const CppADScalar &x) {
    return CppAD::sqrt(x);
}
#endif

inline double scalar_pow(double x, int exponent) {
    return std::pow(x, exponent);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppADScalar scalar_pow(const CppADScalar &x, int exponent) {
    return CppAD::pow(x, static_cast<double>(exponent));
}
#endif

inline double scalar_pow(double x, double exponent) {
    return std::pow(x, exponent);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppADScalar scalar_pow(const CppADScalar &x, double exponent) {
    return CppAD::pow(x, exponent);
}
#endif
