#include "ad_derivative_checks.h"

#ifdef EPCSAFT_HAS_CPPAD
#include <cppad/cppad.hpp>
#endif

namespace epcsaft::native::autodiff {

ADDerivativeResult cppad_square_smoke_derivative(double x) {
#ifdef EPCSAFT_HAS_CPPAD
    std::vector<CppADScalar> ax(1);
    ax[0] = x;
    CppAD::Independent(ax);

    std::vector<CppADScalar> ay(1);
    ay[0] = scalar_pow(ax[0], 2);

    CppAD::ADFun<double> function(ax, ay);
    std::vector<double> point{x};
    std::vector<double> value = function.Forward(0, point);
    std::vector<double> jacobian = function.Jacobian(point);

    ADDerivativeResult result;
    result.supported = true;
    result.backend = "cppad";
    result.message = "CppAD smoke derivative available";
    result.value = value;
    result.jacobian_row_major = jacobian;
    result.rows = 1;
    result.cols = 1;
    return result;
#else
    (void)x;
    ADDerivativeResult result;
    result.supported = false;
    result.backend = "backend_unavailable";
    result.message = "CppAD support is disabled in this native build";
    return result;
#endif
}

}  // namespace epcsaft::native::autodiff
