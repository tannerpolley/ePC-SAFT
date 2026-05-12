#pragma once

#include <string>
#include <vector>

namespace epcsaft::native::autodiff {

struct ADDerivativeResult {
    bool supported = false;
    std::string backend = "backend_unavailable";
    std::string message;
    std::vector<double> value;
    std::vector<double> jacobian_row_major;
    int rows = 0;
    int cols = 0;
};

}  // namespace epcsaft::native::autodiff
