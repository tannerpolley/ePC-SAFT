#pragma once

#include <cmath>

namespace autodiff {

struct dual {
    double val = 0.0;
    double grad = 0.0;

    dual() = default;
    dual(double value, double derivative = 0.0) : val(value), grad(derivative) {}
};

inline dual operator+(const dual &lhs, const dual &rhs) {
    return dual(lhs.val + rhs.val, lhs.grad + rhs.grad);
}

inline dual operator+(const dual &lhs, double rhs) {
    return dual(lhs.val + rhs, lhs.grad);
}

inline dual operator+(double lhs, const dual &rhs) {
    return dual(lhs + rhs.val, rhs.grad);
}

inline dual operator-(const dual &lhs, const dual &rhs) {
    return dual(lhs.val - rhs.val, lhs.grad - rhs.grad);
}

inline dual operator-(const dual &lhs, double rhs) {
    return dual(lhs.val - rhs, lhs.grad);
}

inline dual operator-(double lhs, const dual &rhs) {
    return dual(lhs - rhs.val, -rhs.grad);
}

inline dual operator-(const dual &value) {
    return dual(-value.val, -value.grad);
}

inline dual operator*(const dual &lhs, const dual &rhs) {
    return dual(lhs.val * rhs.val, lhs.grad * rhs.val + lhs.val * rhs.grad);
}

inline dual operator*(const dual &lhs, double rhs) {
    return dual(lhs.val * rhs, lhs.grad * rhs);
}

inline dual operator*(double lhs, const dual &rhs) {
    return dual(lhs * rhs.val, lhs * rhs.grad);
}

inline dual operator/(const dual &lhs, const dual &rhs) {
    return dual(
        lhs.val / rhs.val,
        (lhs.grad * rhs.val - lhs.val * rhs.grad) / (rhs.val * rhs.val)
    );
}

inline dual operator/(const dual &lhs, double rhs) {
    return dual(lhs.val / rhs, lhs.grad / rhs);
}

inline dual operator/(double lhs, const dual &rhs) {
    return dual(lhs / rhs.val, (-lhs * rhs.grad) / (rhs.val * rhs.val));
}

inline dual &operator+=(dual &lhs, const dual &rhs) {
    lhs = lhs + rhs;
    return lhs;
}

inline dual &operator-=(dual &lhs, const dual &rhs) {
    lhs = lhs - rhs;
    return lhs;
}

inline dual &operator*=(dual &lhs, const dual &rhs) {
    lhs = lhs * rhs;
    return lhs;
}

inline dual &operator/=(dual &lhs, const dual &rhs) {
    lhs = lhs / rhs;
    return lhs;
}

inline dual exp(const dual &x) {
    double ev = std::exp(x.val);
    return dual(ev, ev * x.grad);
}

inline dual log(const dual &x) {
    return dual(std::log(x.val), x.grad / x.val);
}

inline dual sqrt(const dual &x) {
    double sv = std::sqrt(x.val);
    return dual(sv, x.grad / (2.0 * sv));
}

inline dual pow(const dual &x, int exponent) {
    if (exponent == 0) {
        return dual(1.0, 0.0);
    }
    double value = std::pow(x.val, exponent);
    return dual(value, exponent * std::pow(x.val, exponent - 1) * x.grad);
}

inline dual pow(const dual &x, double exponent) {
    double value = std::pow(x.val, exponent);
    return dual(value, exponent * std::pow(x.val, exponent - 1.0) * x.grad);
}

inline dual abs(const dual &x) {
    return (x.val >= 0.0) ? dual(x.val, x.grad) : dual(-x.val, -x.grad);
}

inline double value(const dual &x) {
    return x.val;
}

inline double derivative(const dual &x) {
    return x.grad;
}

}  // namespace autodiff
