#include "nlp_problem.h"

#include "epcsaft_electrolyte.h"

#include <cstddef>
#include <cmath>
#include <sstream>

namespace epcsaft::native::equilibrium_nlp {

namespace {

void require_size(const std::vector<double>& values, std::size_t expected, const std::string& label) {
    if (values.size() == expected) {
        return;
    }
    std::ostringstream msg;
    msg << label << " has size " << values.size() << " but expected " << expected << ".";
    throw ValueError(msg.str());
}

void require_optional_size(const std::vector<double>& values, std::size_t expected, const std::string& label) {
    if (values.empty() || values.size() == expected) {
        return;
    }
    std::ostringstream msg;
    msg << label << " has size " << values.size() << " but expected " << expected << " when provided.";
    throw ValueError(msg.str());
}

void require_positive_finite(double value, const std::string& label) {
    if (std::isfinite(value) && value > 0.0) {
        return;
    }
    throw ValueError(label + " must be positive and finite.");
}

void require_positive_finite(const std::vector<double>& values, const std::string& label) {
    for (double value : values) {
        require_positive_finite(value, label);
    }
}

}  // namespace

void validate_nlp_problem_shape(const NlpProblem& problem) {
    const int variables = problem.variable_count();
    const int constraints = problem.constraint_count();
    const int jacobian_nonzeros = problem.jacobian_nonzero_count();
    const int hessian_nonzeros = problem.hessian_nonzero_count();
    if (variables <= 0) {
        throw ValueError("NLP problem must expose at least one variable.");
    }
    if (constraints < 0) {
        throw ValueError("NLP problem constraint count cannot be negative.");
    }
    if (jacobian_nonzeros < 0) {
        throw ValueError("NLP problem Jacobian nonzero count cannot be negative.");
    }
    if (hessian_nonzeros < 0) {
        throw ValueError("NLP problem Hessian nonzero count cannot be negative.");
    }

    const NlpBounds problem_bounds = problem.bounds();
    require_size(problem_bounds.variable_lower, static_cast<std::size_t>(variables), "NLP variable lower bounds");
    require_size(problem_bounds.variable_upper, static_cast<std::size_t>(variables), "NLP variable upper bounds");
    require_size(problem_bounds.constraint_lower, static_cast<std::size_t>(constraints), "NLP constraint lower bounds");
    require_size(problem_bounds.constraint_upper, static_cast<std::size_t>(constraints), "NLP constraint upper bounds");
    require_size(problem.initial_point(), static_cast<std::size_t>(variables), "NLP initial point");

    const NlpJacobianStructure structure = problem.jacobian_structure();
    if (structure.rows.size() != static_cast<std::size_t>(jacobian_nonzeros)
        || structure.cols.size() != static_cast<std::size_t>(jacobian_nonzeros)) {
        throw ValueError("NLP Jacobian structure must match jacobian_nonzero_count.");
    }
    for (std::size_t index = 0; index < structure.rows.size(); ++index) {
        if (structure.rows[index] < 0 || structure.rows[index] >= constraints) {
            throw ValueError("NLP Jacobian row index is out of range.");
        }
        if (structure.cols[index] < 0 || structure.cols[index] >= variables) {
            throw ValueError("NLP Jacobian column index is out of range.");
        }
    }
    if (problem.has_exact_hessian()) {
        if (hessian_nonzeros <= 0) {
            throw ValueError("NLP exact Hessian support requires positive hessian_nonzero_count.");
        }
        const NlpHessianStructure hessian = problem.hessian_structure();
        if (hessian.rows.size() != static_cast<std::size_t>(hessian_nonzeros)
            || hessian.cols.size() != static_cast<std::size_t>(hessian_nonzeros)) {
            throw ValueError("NLP Hessian structure must match hessian_nonzero_count.");
        }
        for (std::size_t index = 0; index < hessian.rows.size(); ++index) {
            if (hessian.rows[index] < 0 || hessian.rows[index] >= variables) {
                throw ValueError("NLP Hessian row index is out of range.");
            }
            if (hessian.cols[index] < 0 || hessian.cols[index] >= variables) {
                throw ValueError("NLP Hessian column index is out of range.");
            }
        }
    } else if (hessian_nonzeros != 0) {
        throw ValueError("NLP Hessian nonzero count requires exact Hessian support.");
    }

    const NlpScaling scaling = problem.scaling();
    require_positive_finite(scaling.objective, "NLP objective scaling");
    require_optional_size(scaling.variables, static_cast<std::size_t>(variables), "NLP variable scaling");
    require_optional_size(scaling.constraints, static_cast<std::size_t>(constraints), "NLP constraint scaling");
    require_positive_finite(scaling.variables, "NLP variable scaling");
    require_positive_finite(scaling.constraints, "NLP constraint scaling");
}

}  // namespace epcsaft::native::equilibrium_nlp
