#include "ad_derivative_checks.h"

#include "ad_scalar.h"

#include <algorithm>
#include <cmath>

namespace epcsaft::autodiff {

NativeAutodiffDerivativeCheckResult native_autodiff_derivative_checks() {
    NativeAutodiffDerivativeCheckResult result;
    result.cppad_compiled = cppad_compiled();
    result.checked_residuals = {
        "scaled_residual",
        "pressure_log_residual",
        "reaction_log_residual",
    };

#ifndef EPCSAFT_HAS_CPPAD
    result.status = "unsupported_derivative";
    result.derivative_backend = "cppad_unavailable";
    return result;
#else
    std::vector<CppADScalar> independent(1);
    independent[0] = 2.0;
    CppAD::Independent(independent);

    const CppADScalar x = independent[0];
    std::vector<CppADScalar> dependent(3);
    dependent[0] = scaled_residual(x, CppADScalar(1.25), CppADScalar(4.0));
    dependent[1] = pressure_log_residual(x, CppADScalar(1.5));
    dependent[2] = reaction_log_residual(CppAD::log(x), CppADScalar(0.25));

    CppAD::ADFun<double> tape(independent, dependent);
    const std::vector<double> jacobian = tape.Jacobian(std::vector<double>{2.0});

    result.cppad_used = true;
    result.status = "ok";
    result.derivative_backend = "cppad";
    result.derivative_by_residual["scaled_residual"] = jacobian[0];
    result.derivative_by_residual["pressure_log_residual"] = jacobian[1];
    result.derivative_by_residual["reaction_log_residual"] = jacobian[2];

    const double expected_scaled = 4.0;
    const double expected_pressure_log = 0.5;
    const double expected_reaction_log = 0.5;
    result.max_abs_error = std::max(
        {std::abs(jacobian[0] - expected_scaled),
         std::abs(jacobian[1] - expected_pressure_log),
         std::abs(jacobian[2] - expected_reaction_log)}
    );
    return result;
#endif
}

}  // namespace epcsaft::autodiff


