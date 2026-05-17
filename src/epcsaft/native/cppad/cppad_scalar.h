#pragma once

#include <cmath>
#include <string>

#ifndef EPCSAFT_HAS_CPPAD
#error "EPCSAFT_HAS_CPPAD must be defined; CppAD is a required native dependency."
#endif

#include <cppad/cppad.hpp>

namespace epcsaft::native::cppad_support {

inline bool cppad_compiled() {
    return true;
}

inline std::string cppad_build_status() {
#ifdef EPCSAFT_CPPAD_STATUS
    return EPCSAFT_CPPAD_STATUS;
#else
    return "enabled_available";
#endif
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

}  // namespace epcsaft::native::cppad_support
