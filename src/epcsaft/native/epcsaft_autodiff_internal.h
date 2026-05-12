#pragma once

#include <cmath>
#include <string>

#ifdef EPCSAFT_HAS_CPPAD
#include <cppad/cppad.hpp>
#endif

#include <unsupported/Eigen/AutoDiff>

using AutoDualDerivative = Eigen::Matrix<double, 1, 1>;
using AutoDual = Eigen::AutoDiffScalar<AutoDualDerivative>;

enum class DerivativeBackend {
    Analytic = 0,
    BackendUnavailable = 1,
    Autodiff = 2,
    Auto = 3,
};

inline DerivativeBackend derivative_backend_from_int(int mode) {
    switch (mode) {
        case 0:
            return DerivativeBackend::Analytic;
        case 1:
            return DerivativeBackend::BackendUnavailable;
        case 2:
            return DerivativeBackend::Autodiff;
        case 3:
            return DerivativeBackend::Auto;
        default:
            return DerivativeBackend::Auto;
    }
}

inline std::string derivative_backend_name(DerivativeBackend backend) {
    switch (backend) {
        case DerivativeBackend::Analytic:
            return "analytic";
        case DerivativeBackend::BackendUnavailable:
            return "backend_unavailable";
        case DerivativeBackend::Autodiff:
            return "autodiff";
        case DerivativeBackend::Auto:
            return "auto";
    }
    return "auto";
}

inline AutoDual make_autodiff_scalar(double value, double derivative = 0.0) {
    AutoDual x;
    x.value() = value;
    x.derivatives() = AutoDualDerivative::Constant(derivative);
    return x;
}

template <typename Scalar>
inline Scalar scalar_constant(double value) {
    return static_cast<Scalar>(value);
}

template <>
inline AutoDual scalar_constant<AutoDual>(double value) {
    return make_autodiff_scalar(value, 0.0);
}

inline double scalar_value(double x) {
    return x;
}

inline double scalar_value(const AutoDual &x) {
    return x.value();
}

#ifdef EPCSAFT_HAS_CPPAD
inline double scalar_value(const CppAD::AD<double> &x) {
    return CppAD::Value(x);
}
#endif

inline double scalar_derivative(double) {
    return 0.0;
}

inline double scalar_derivative(const AutoDual &x) {
    return x.derivatives().size() == 0 ? 0.0 : x.derivatives()[0];
}

inline double scalar_log(double x) {
    return std::log(x);
}

inline AutoDual scalar_log(const AutoDual &x) {
    using std::log;
    return log(x);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppAD::AD<double> scalar_log(const CppAD::AD<double> &x) {
    return CppAD::log(x);
}
#endif

inline double scalar_exp(double x) {
    return std::exp(x);
}

inline AutoDual scalar_exp(const AutoDual &x) {
    using std::exp;
    return exp(x);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppAD::AD<double> scalar_exp(const CppAD::AD<double> &x) {
    return CppAD::exp(x);
}
#endif

inline double scalar_sqrt(double x) {
    return std::sqrt(x);
}

inline AutoDual scalar_sqrt(const AutoDual &x) {
    using std::sqrt;
    return sqrt(x);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppAD::AD<double> scalar_sqrt(const CppAD::AD<double> &x) {
    return CppAD::sqrt(x);
}
#endif

inline double scalar_pow(double x, int exponent) {
    return std::pow(x, exponent);
}

inline AutoDual scalar_pow(const AutoDual &x, int exponent) {
    using std::pow;
    return pow(x, static_cast<double>(exponent));
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppAD::AD<double> scalar_pow(const CppAD::AD<double> &x, int exponent) {
    return CppAD::pow(x, static_cast<double>(exponent));
}
#endif

inline double scalar_pow(double x, double exponent) {
    return std::pow(x, exponent);
}

inline AutoDual scalar_pow(const AutoDual &x, double exponent) {
    using std::pow;
    return pow(x, exponent);
}

#ifdef EPCSAFT_HAS_CPPAD
inline CppAD::AD<double> scalar_pow(const CppAD::AD<double> &x, double exponent) {
    return CppAD::pow(x, exponent);
}
#endif
