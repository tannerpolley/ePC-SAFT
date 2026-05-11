#pragma once

#include <string>
#include <vector>

struct NativeImplicitSensitivityResult {
    bool success = false;
    std::string status = "invalid_input";
    std::string message;
    std::vector<double> sensitivities_row_major;
    int rows = 0;
    int cols = 0;
    double residual_jacobian_condition_proxy = 0.0;
};

NativeImplicitSensitivityResult solve_native_implicit_sensitivity(
    const std::vector<double>& residual_jacobian_u_row_major,
    int residual_rows,
    int state_cols,
    const std::vector<double>& residual_jacobian_theta_row_major,
    int parameter_cols
);
