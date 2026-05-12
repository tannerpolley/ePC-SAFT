#include "implicit_sensitivity.h"

#include <Eigen/Dense>

#include <cmath>

NativeImplicitSensitivityResult solve_native_implicit_sensitivity(
    const std::vector<double>& residual_jacobian_u_row_major,
    int residual_rows,
    int state_cols,
    const std::vector<double>& residual_jacobian_theta_row_major,
    int parameter_cols
) {
    NativeImplicitSensitivityResult result;
    result.rows = state_cols;
    result.cols = parameter_cols;
    if (residual_rows <= 0 || state_cols <= 0 || parameter_cols <= 0) {
        result.status = "invalid_input";
        result.message = "implicit sensitivity dimensions must be positive.";
        return result;
    }
    if (residual_rows < state_cols) {
        result.status = "unsupported_derivative";
        result.message = "implicit sensitivity requires at least as many residual equations as state variables.";
        return result;
    }
    if (residual_jacobian_u_row_major.size() != static_cast<std::size_t>(residual_rows * state_cols)
        || residual_jacobian_theta_row_major.size() != static_cast<std::size_t>(residual_rows * parameter_cols)) {
        result.status = "invalid_input";
        result.message = "implicit sensitivity Jacobian payload sizes do not match dimensions.";
        return result;
    }
    Eigen::MatrixXd r_u(residual_rows, state_cols);
    Eigen::MatrixXd r_theta(residual_rows, parameter_cols);
    for (int row = 0; row < residual_rows; ++row) {
        for (int col = 0; col < state_cols; ++col) {
            const double value = residual_jacobian_u_row_major[static_cast<std::size_t>(row * state_cols + col)];
            if (!std::isfinite(value)) {
                result.status = "nonfinite_objective";
                result.message = "implicit sensitivity R_u contains a non-finite value.";
                return result;
            }
            r_u(row, col) = value;
        }
        for (int col = 0; col < parameter_cols; ++col) {
            const double value = residual_jacobian_theta_row_major[static_cast<std::size_t>(row * parameter_cols + col)];
            if (!std::isfinite(value)) {
                result.status = "nonfinite_objective";
                result.message = "implicit sensitivity R_theta contains a non-finite value.";
                return result;
            }
            r_theta(row, col) = value;
        }
    }
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(r_u);
    if (qr.rank() < state_cols) {
        result.status = "singular_jacobian";
        result.message = "implicit sensitivity inner residual Jacobian is singular.";
        return result;
    }
    const Eigen::MatrixXd u_theta = qr.solve(-r_theta);
    result.sensitivities_row_major.reserve(static_cast<std::size_t>(state_cols * parameter_cols));
    for (int row = 0; row < u_theta.rows(); ++row) {
        for (int col = 0; col < u_theta.cols(); ++col) {
            result.sensitivities_row_major.push_back(u_theta(row, col));
        }
    }
    result.residual_jacobian_condition_proxy = static_cast<double>(qr.rank());
    result.success = true;
    result.status = "converged";
    result.message = "implicit sensitivity solved from converged inner residual Jacobians.";
    return result;
}


