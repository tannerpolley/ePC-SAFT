#pragma once

#include <map>
#include <string>
#include <vector>

namespace epcsaft::autodiff {

struct NativeAutodiffDerivativeCheckResult {
    bool cppad_compiled = false;
    bool cppad_used = false;
    bool unsupported_derivative_used = false;
    std::string status = "unsupported_derivative";
    std::string derivative_backend = "cppad_unavailable";
    std::vector<std::string> checked_residuals;
    std::map<std::string, double> derivative_by_residual;
    double max_abs_error = 0.0;
};

NativeAutodiffDerivativeCheckResult native_autodiff_derivative_checks();

}  // namespace epcsaft::autodiff



