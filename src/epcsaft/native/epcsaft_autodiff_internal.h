#pragma once

#include "third_party/autodiff/dual.hpp"

using AutoDual = autodiff::dual;

inline double scalar_value(double x) {
    return x;
}

inline double scalar_value(const AutoDual &x) {
    return autodiff::value(x);
}

inline double scalar_derivative(double) {
    return 0.0;
}

inline double scalar_derivative(const AutoDual &x) {
    return autodiff::derivative(x);
}

inline double scalar_log(double x) {
    return std::log(x);
}

inline AutoDual scalar_log(const AutoDual &x) {
    return autodiff::log(x);
}

inline double scalar_exp(double x) {
    return std::exp(x);
}

inline AutoDual scalar_exp(const AutoDual &x) {
    return autodiff::exp(x);
}

inline double scalar_sqrt(double x) {
    return std::sqrt(x);
}

inline AutoDual scalar_sqrt(const AutoDual &x) {
    return autodiff::sqrt(x);
}

inline double scalar_pow(double x, int exponent) {
    return std::pow(x, exponent);
}

inline AutoDual scalar_pow(const AutoDual &x, int exponent) {
    return autodiff::pow(x, exponent);
}

inline double scalar_pow(double x, double exponent) {
    return std::pow(x, exponent);
}

inline AutoDual scalar_pow(const AutoDual &x, double exponent) {
    return autodiff::pow(x, exponent);
}
