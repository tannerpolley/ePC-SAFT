#pragma once

#include <string>

namespace epcsaft::native::regression {

struct RegressionProblemContract {
    std::string optimizer_backend = "ceres";
    std::string derivative_backend = "legacy_eigen_forward";
    bool python_objective_used = false;
};

}  // namespace epcsaft::native::regression

