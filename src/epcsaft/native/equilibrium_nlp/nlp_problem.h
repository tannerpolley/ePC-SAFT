#pragma once

#include <map>
#include <string>
#include <vector>

namespace epcsaft::native::equilibrium_nlp {

struct NlpBounds {
    std::vector<double> variable_lower;
    std::vector<double> variable_upper;
    std::vector<double> constraint_lower;
    std::vector<double> constraint_upper;
};

struct NlpJacobianStructure {
    std::vector<int> rows;
    std::vector<int> cols;
};

struct NlpHessianStructure {
    std::vector<int> rows;
    std::vector<int> cols;
};

struct NlpScaling {
    double objective = 1.0;
    std::vector<double> variables;
    std::vector<double> constraints;
};

class NlpProblem {
public:
    virtual ~NlpProblem() = default;

    virtual std::string name() const = 0;
    virtual int variable_count() const = 0;
    virtual int constraint_count() const = 0;
    virtual int jacobian_nonzero_count() const = 0;

    virtual NlpBounds bounds() const = 0;
    virtual std::vector<double> initial_point() const = 0;
    virtual double objective(const std::vector<double>& variables) const = 0;
    virtual std::vector<double> objective_gradient(const std::vector<double>& variables) const = 0;
    virtual std::vector<double> constraints(const std::vector<double>& variables) const = 0;
    virtual NlpJacobianStructure jacobian_structure() const = 0;
    virtual std::vector<double> jacobian_values(const std::vector<double>& variables) const = 0;

    virtual bool has_exact_hessian() const {
        return false;
    }

    virtual int hessian_nonzero_count() const {
        return 0;
    }

    virtual NlpHessianStructure hessian_structure() const {
        return {};
    }

    virtual std::vector<double> hessian_values(
        const std::vector<double>& variables,
        double objective_factor,
        const std::vector<double>& constraint_multipliers
    ) const {
        (void)variables;
        (void)objective_factor;
        (void)constraint_multipliers;
        return {};
    }

    virtual std::string hessian_backend() const {
        return "";
    }

    virtual NlpScaling scaling() const {
        return {};
    }

    virtual std::map<std::string, std::string> diagnostics() const {
        return {};
    }
};

void validate_nlp_problem_shape(const NlpProblem& problem);

}  // namespace epcsaft::native::equilibrium_nlp
