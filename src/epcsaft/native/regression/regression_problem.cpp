#include "regression_problem.h"

namespace epcsaft::native::regression {

bool regression_problem_uses_python_objective_loop() {
    return RegressionProblemContract{}.python_objective_used;
}

}  // namespace epcsaft::native::regression
