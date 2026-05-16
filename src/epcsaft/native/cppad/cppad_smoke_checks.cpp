#include "cppad_smoke_checks.h"

#ifdef EPCSAFT_HAS_CPPAD
#include <cppad/cppad.hpp>
#endif

namespace epcsaft::native::cppad_support {

CppADDerivativeResult cppad_square_smoke_derivative(double x) {
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

    CppADDerivativeResult result;
    result.supported = true;
    result.backend = "cppad";
    result.message = "CppAD smoke derivative available";
    result.value = value;
    result.jacobian_row_major = jacobian;
    result.outputs = {"x_squared"};
    result.variables = {"x"};
    result.rows = 1;
    result.cols = 1;
    return result;
#else
    (void)x;
    CppADDerivativeResult result;
    result.supported = false;
    result.backend = cppad_disabled_status();
    result.message = "CppAD support is disabled in this native build";
    return result;
#endif
}

}  // namespace epcsaft::native::cppad_support
