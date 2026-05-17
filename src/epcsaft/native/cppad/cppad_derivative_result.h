#pragma once

#include <string>
#include <vector>

namespace epcsaft::native::cppad_support {

struct CppADDerivativeResult {
    bool supported = false;
    std::string backend = "cppad";
    std::string message;
    std::vector<double> value;
    std::vector<double> jacobian_row_major;
    std::vector<std::string> outputs;
    std::vector<std::string> variables;
    int rows = 0;
    int cols = 0;
};

}  // namespace epcsaft::native::cppad_support
