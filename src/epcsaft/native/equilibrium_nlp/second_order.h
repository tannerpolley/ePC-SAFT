#pragma once

#include "nlp_problem.h"

#include <string>
#include <vector>

namespace epcsaft::native::equilibrium_nlp {

struct ObjectiveSecondOrderData {
    int variable_count = 0;
    double value = 0.0;
    std::vector<double> gradient;
    std::vector<double> hessian_row_major;
    std::string backend;
};

struct ConstraintSecondOrderData {
    int constraint_count = 0;
    int variable_count = 0;
    std::vector<double> values;
    std::vector<double> jacobian_row_major;
    std::vector<double> hessian_tensor_row_major;
    std::vector<bool> has_hessian;
    std::string backend;
};

struct ResidualSecondOrderData {
    int residual_count = 0;
    int variable_count = 0;
    std::vector<double> residuals;
    std::vector<double> jacobian_row_major;
    std::vector<double> residual_hessian_tensor_row_major;
    std::string backend;
};

struct VariableTransformSecondOrderData {
    int input_variable_count = 0;
    int output_variable_count = 0;
    std::vector<double> jacobian_row_major;
    std::vector<double> output_hessian_tensor_row_major;
    std::string backend;
};

int dense_lower_hessian_nonzero_count(int variable_count);

NlpHessianStructure dense_lower_hessian_structure(int variable_count);

std::vector<double> lower_triangle_values(
    const std::vector<double>& dense_row_major,
    int variable_count
);

ObjectiveSecondOrderData least_squares_objective_second_order(
    const ResidualSecondOrderData& residuals
);

ObjectiveSecondOrderData transformed_objective_second_order(
    const VariableTransformSecondOrderData& transform,
    const ObjectiveSecondOrderData& output_objective
);

ConstraintSecondOrderData transformed_constraint_second_order(
    const VariableTransformSecondOrderData& transform,
    const ConstraintSecondOrderData& output_constraints
);

class LagrangianHessianAssembler {
public:
    explicit LagrangianHessianAssembler(int variable_count);

    int variable_count() const {
        return variable_count_;
    }

    int nonzero_count() const;

    NlpHessianStructure structure() const;

    std::vector<double> values(
        double objective_factor,
        const ObjectiveSecondOrderData& objective,
        const ConstraintSecondOrderData& constraints,
        const std::vector<double>& constraint_multipliers
    ) const;

private:
    int variable_count_;
};

}  // namespace epcsaft::native::equilibrium_nlp
