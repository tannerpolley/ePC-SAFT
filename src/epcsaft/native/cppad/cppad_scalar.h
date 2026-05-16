#pragma once

#include <cmath>
#include <string>

#ifdef EPCSAFT_HAS_CPPAD
#include <cppad/cppad.hpp>
#endif

namespace epcsaft::native::cppad_support {

inline bool cppad_compiled() {
#ifdef EPCSAFT_HAS_CPPAD
    return true;
#else
    return false;
#endif
}

inline std::string cppad_build_status() {
#ifdef EPCSAFT_CPPAD_STATUS
    return EPCSAFT_CPPAD_STATUS;
#else
    return "not_configured";
#endif
}

inline std::string cppad_disabled_status() {
    return "cppad_disabled";
}

inline double scalar_value(double x) {
    return x;
}

inline double scalar_log(double x) {
    return std::log(x);
}

inline double scalar_exp(double x) {
    return std::exp(x);
}

inline double scalar_pow(double x, double exponent) {
    return std::pow(x, exponent);
}

inline double scalar_pow(double x, int exponent) {
    return std::pow(x, exponent);
}

#ifdef EPCSAFT_HAS_CPPAD
using CppADScalar = CppAD::AD<double>;

inline CppADScalar make_cppad_scalar(double value) {
    return CppADScalar(value);
}

inline double scalar_value(const CppADScalar& x) {
    return CppAD::Value(x);
}

inline CppADScalar scalar_log(const CppADScalar& x) {
    return CppAD::log(x);
}

inline CppADScalar scalar_exp(const CppADScalar& x) {
    return CppAD::exp(x);
}

inline CppADScalar scalar_pow(const CppADScalar& x, double exponent) {
    return CppAD::pow(x, exponent);
}

inline CppADScalar scalar_pow(const CppADScalar& x, int exponent) {
    return CppAD::pow(x, static_cast<double>(exponent));
}
#endif

}  // namespace epcsaft::native::cppad_support
