#include "second_order.h"

#include "epcsaft_electrolyte.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <sstream>

namespace epcsaft::native::equilibrium_nlp {

namespace {

constexpr double kSymmetryTolerance = 1.0e-8;

void require_nonnegative_count(int value, const std::string& label) {
    if (value >= 0) {
        return;
    }
    throw ValueError(label + " cannot be negative.");
}

void require_positive_count(int value, const std::string& label) {
    if (value > 0) {
        return;
    }
    throw ValueError(label + " must be positive.");
}

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
    require_size(values, expected, label);
}

void require_bool_optional_size(const std::vector<bool>& values, std::size_t expected, const std::string& label) {
    if (values.empty() || values.size() == expected) {
        return;
    }
    std::ostringstream msg;
    msg << label << " has size " << values.size() << " but expected " << expected << ".";
    throw ValueError(msg.str());
}

void require_finite(double value, const std::string& label) {
    if (std::isfinite(value)) {
        return;
    }
    throw ValueError(label + " must be finite.");
}

void require_finite_values(const std::vector<double>& values, const std::string& label) {
    for (double value : values) {
        require_finite(value, label);
    }
}

void require_matrix_shape(const std::vector<double>& values, int rows, int cols, const std::string& label) {
    require_nonnegative_count(rows, label + " row count");
    require_nonnegative_count(cols, label + " column count");
    require_size(
        values,
        static_cast<std::size_t>(rows) * static_cast<std::size_t>(cols),
        label
    );
    require_finite_values(values, label);
}

void require_square_matrix_shape(const std::vector<double>& values, int variable_count, const std::string& label) {
    require_matrix_shape(values, variable_count, variable_count, label);
}

void require_tensor_shape(
    const std::vector<double>& values,
    int rows,
    int variable_count,
    const std::string& label
) {
    require_nonnegative_count(rows, label + " row count");
    require_nonnegative_count(variable_count, label + " variable count");
    require_size(
        values,
        static_cast<std::size_t>(rows)
            * static_cast<std::size_t>(variable_count)
            * static_cast<std::size_t>(variable_count),
        label
    );
    require_finite_values(values, label);
}

void validate_symmetric_matrix(const std::vector<double>& values, int variable_count, const std::string& label) {
    require_square_matrix_shape(values, variable_count, label);
    const std::size_t n = static_cast<std::size_t>(variable_count);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < row; ++col) {
            const double upper = values[col * n + row];
            const double lower = values[row * n + col];
            const double scale = std::max({1.0, std::abs(upper), std::abs(lower)});
            if (std::abs(upper - lower) > kSymmetryTolerance * scale) {
                throw ValueError(label + " must be symmetric.");
            }
        }
    }
}

void validate_symmetric_tensor(
    const std::vector<double>& values,
    int rows,
    int variable_count,
    const std::string& label
) {
    require_tensor_shape(values, rows, variable_count, label);
    const std::size_t n = static_cast<std::size_t>(variable_count);
    for (int tensor_row = 0; tensor_row < rows; ++tensor_row) {
        const std::size_t offset = static_cast<std::size_t>(tensor_row) * n * n;
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t col = 0; col < row; ++col) {
                const double upper = values[offset + col * n + row];
                const double lower = values[offset + row * n + col];
                const double scale = std::max({1.0, std::abs(upper), std::abs(lower)});
                if (std::abs(upper - lower) > kSymmetryTolerance * scale) {
                    throw ValueError(label + " Hessian blocks must be symmetric.");
                }
            }
        }
    }
}

std::vector<bool> active_hessian_rows(
    const std::vector<bool>& provided,
    int row_count,
    bool default_value
) {
    if (!provided.empty()) {
        require_bool_optional_size(provided, static_cast<std::size_t>(row_count), "second-order active row flags");
        return provided;
    }
    return std::vector<bool>(static_cast<std::size_t>(row_count), default_value);
}

void add_scaled_matrix(
    std::vector<double>& target,
    const std::vector<double>& source,
    int variable_count,
    double scale
) {
    if (scale == 0.0) {
        return;
    }
    require_square_matrix_shape(source, variable_count, "second-order source Hessian");
    if (target.size() != source.size()) {
        throw ValueError("second-order target Hessian shape does not match source Hessian.");
    }
    for (std::size_t index = 0; index < target.size(); ++index) {
        target[index] += scale * source[index];
    }
}

void symmetrize_square_matrix(std::vector<double>& values, int variable_count) {
    const std::size_t n = static_cast<std::size_t>(variable_count);
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col < row; ++col) {
            const double value = 0.5 * (values[row * n + col] + values[col * n + row]);
            values[row * n + col] = value;
            values[col * n + row] = value;
        }
    }
}

}  // namespace

int dense_lower_hessian_nonzero_count(int variable_count) {
    require_nonnegative_count(variable_count, "dense lower Hessian variable count");
    return variable_count * (variable_count + 1) / 2;
}

NlpHessianStructure dense_lower_hessian_structure(int variable_count) {
    require_nonnegative_count(variable_count, "dense lower Hessian variable count");
    NlpHessianStructure out;
    out.rows.reserve(static_cast<std::size_t>(dense_lower_hessian_nonzero_count(variable_count)));
    out.cols.reserve(static_cast<std::size_t>(dense_lower_hessian_nonzero_count(variable_count)));
    for (int row = 0; row < variable_count; ++row) {
        for (int col = 0; col <= row; ++col) {
            out.rows.push_back(row);
            out.cols.push_back(col);
        }
    }
    return out;
}

std::vector<double> lower_triangle_values(
    const std::vector<double>& dense_row_major,
    int variable_count
) {
    validate_symmetric_matrix(dense_row_major, variable_count, "dense Hessian");
    const std::size_t n = static_cast<std::size_t>(variable_count);
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(dense_lower_hessian_nonzero_count(variable_count)));
    for (std::size_t row = 0; row < n; ++row) {
        for (std::size_t col = 0; col <= row; ++col) {
            out.push_back(dense_row_major[row * n + col]);
        }
    }
    return out;
}

ObjectiveSecondOrderData residual_quadratic_objective_second_order(
    const ResidualSecondOrderData& residuals
) {
    require_positive_count(residuals.variable_count, "residual quadratic variable count");
    require_nonnegative_count(residuals.residual_count, "residual quadratic residual count");
    require_size(
        residuals.residuals,
        static_cast<std::size_t>(residuals.residual_count),
        "residual quadratic residual vector"
    );
    require_matrix_shape(
        residuals.jacobian_row_major,
        residuals.residual_count,
        residuals.variable_count,
        "residual quadratic residual Jacobian"
    );
    validate_symmetric_tensor(
        residuals.residual_hessian_tensor_row_major,
        residuals.residual_count,
        residuals.variable_count,
        "residual quadratic residual Hessian tensor"
    );
    require_finite_values(residuals.residuals, "residual quadratic residual vector");

    const std::size_t rows = static_cast<std::size_t>(residuals.residual_count);
    const std::size_t n = static_cast<std::size_t>(residuals.variable_count);
    ObjectiveSecondOrderData out;
    out.variable_count = residuals.variable_count;
    out.backend = residuals.backend;
    out.gradient.assign(n, 0.0);
    out.hessian_row_major.assign(n * n, 0.0);
    for (std::size_t residual_row = 0; residual_row < rows; ++residual_row) {
        const double residual = residuals.residuals[residual_row];
        out.value += 0.5 * residual * residual;
        for (std::size_t col = 0; col < n; ++col) {
            out.gradient[col] += residual * residuals.jacobian_row_major[residual_row * n + col];
        }
        for (std::size_t first = 0; first < n; ++first) {
            for (std::size_t second = 0; second < n; ++second) {
                out.hessian_row_major[first * n + second] +=
                    residuals.jacobian_row_major[residual_row * n + first]
                    * residuals.jacobian_row_major[residual_row * n + second];
                out.hessian_row_major[first * n + second] += residual
                    * residuals.residual_hessian_tensor_row_major[
                        residual_row * n * n + first * n + second
                    ];
            }
        }
    }
    validate_symmetric_matrix(out.hessian_row_major, residuals.variable_count, "residual quadratic objective Hessian");
    return out;
}

ObjectiveSecondOrderData transformed_objective_second_order(
    const VariableTransformSecondOrderData& transform,
    const ObjectiveSecondOrderData& output_objective
) {
    require_positive_count(transform.input_variable_count, "transform input variable count");
    require_positive_count(transform.output_variable_count, "transform output variable count");
    require_matrix_shape(
        transform.jacobian_row_major,
        transform.output_variable_count,
        transform.input_variable_count,
        "transform Jacobian"
    );
    validate_symmetric_tensor(
        transform.output_hessian_tensor_row_major,
        transform.output_variable_count,
        transform.input_variable_count,
        "transform output Hessian tensor"
    );
    if (output_objective.variable_count != transform.output_variable_count) {
        throw ValueError("transformed objective output variable count does not match transform output count.");
    }
    require_optional_size(
        output_objective.gradient,
        static_cast<std::size_t>(transform.output_variable_count),
        "transformed objective gradient"
    );
    validate_symmetric_matrix(
        output_objective.hessian_row_major,
        transform.output_variable_count,
        "transformed objective Hessian"
    );

    const std::size_t input = static_cast<std::size_t>(transform.input_variable_count);
    const std::size_t output = static_cast<std::size_t>(transform.output_variable_count);
    ObjectiveSecondOrderData out;
    out.variable_count = transform.input_variable_count;
    out.value = output_objective.value;
    out.backend = output_objective.backend.empty()
        ? transform.backend
        : output_objective.backend + "_through_" + transform.backend;
    out.gradient.assign(input, 0.0);
    out.hessian_row_major.assign(input * input, 0.0);

    if (!output_objective.gradient.empty()) {
        for (std::size_t out_index = 0; out_index < output; ++out_index) {
            for (std::size_t in_col = 0; in_col < input; ++in_col) {
                out.gradient[in_col] += output_objective.gradient[out_index]
                    * transform.jacobian_row_major[out_index * input + in_col];
            }
        }
    }

    for (std::size_t first = 0; first < input; ++first) {
        for (std::size_t second = 0; second < input; ++second) {
            double value = 0.0;
            for (std::size_t out_first = 0; out_first < output; ++out_first) {
                for (std::size_t out_second = 0; out_second < output; ++out_second) {
                    value += transform.jacobian_row_major[out_first * input + first]
                        * output_objective.hessian_row_major[out_first * output + out_second]
                        * transform.jacobian_row_major[out_second * input + second];
                }
                if (!output_objective.gradient.empty()) {
                    value += output_objective.gradient[out_first]
                        * transform.output_hessian_tensor_row_major[
                            out_first * input * input + first * input + second
                        ];
                }
            }
            out.hessian_row_major[first * input + second] = value;
        }
    }
    validate_symmetric_matrix(out.hessian_row_major, transform.input_variable_count, "transformed objective Hessian");
    return out;
}

ConstraintSecondOrderData transformed_constraint_second_order(
    const VariableTransformSecondOrderData& transform,
    const ConstraintSecondOrderData& output_constraints
) {
    require_positive_count(transform.input_variable_count, "transform input variable count");
    require_positive_count(transform.output_variable_count, "transform output variable count");
    require_matrix_shape(
        transform.jacobian_row_major,
        transform.output_variable_count,
        transform.input_variable_count,
        "transform Jacobian"
    );
    validate_symmetric_tensor(
        transform.output_hessian_tensor_row_major,
        transform.output_variable_count,
        transform.input_variable_count,
        "transform output Hessian tensor"
    );
    if (output_constraints.variable_count != transform.output_variable_count) {
        throw ValueError("transformed constraint output variable count does not match transform output count.");
    }
    require_nonnegative_count(output_constraints.constraint_count, "transformed constraint count");
    require_matrix_shape(
        output_constraints.jacobian_row_major,
        output_constraints.constraint_count,
        output_constraints.variable_count,
        "transformed constraint Jacobian"
    );
    validate_symmetric_tensor(
        output_constraints.hessian_tensor_row_major,
        output_constraints.constraint_count,
        output_constraints.variable_count,
        "transformed constraint Hessian tensor"
    );
    const std::vector<bool> active_rows = active_hessian_rows(
        output_constraints.has_hessian,
        output_constraints.constraint_count,
        true
    );

    const std::size_t input = static_cast<std::size_t>(transform.input_variable_count);
    const std::size_t output = static_cast<std::size_t>(transform.output_variable_count);
    const std::size_t constraints = static_cast<std::size_t>(output_constraints.constraint_count);
    ConstraintSecondOrderData out;
    out.constraint_count = output_constraints.constraint_count;
    out.variable_count = transform.input_variable_count;
    out.values = output_constraints.values;
    out.has_hessian = active_rows;
    out.backend = output_constraints.backend.empty()
        ? transform.backend
        : output_constraints.backend + "_through_" + transform.backend;
    out.jacobian_row_major.assign(constraints * input, 0.0);
    out.hessian_tensor_row_major.assign(constraints * input * input, 0.0);

    for (std::size_t constraint = 0; constraint < constraints; ++constraint) {
        for (std::size_t in_col = 0; in_col < input; ++in_col) {
            for (std::size_t out_col = 0; out_col < output; ++out_col) {
                out.jacobian_row_major[constraint * input + in_col] +=
                    output_constraints.jacobian_row_major[constraint * output + out_col]
                    * transform.jacobian_row_major[out_col * input + in_col];
            }
        }
        if (!active_rows[constraint]) {
            continue;
        }
        for (std::size_t first = 0; first < input; ++first) {
            for (std::size_t second = 0; second < input; ++second) {
                double value = 0.0;
                for (std::size_t out_first = 0; out_first < output; ++out_first) {
                    for (std::size_t out_second = 0; out_second < output; ++out_second) {
                        value += transform.jacobian_row_major[out_first * input + first]
                            * output_constraints.hessian_tensor_row_major[
                                constraint * output * output + out_first * output + out_second
                            ]
                            * transform.jacobian_row_major[out_second * input + second];
                    }
                    value += output_constraints.jacobian_row_major[constraint * output + out_first]
                        * transform.output_hessian_tensor_row_major[
                            out_first * input * input + first * input + second
                        ];
                }
                out.hessian_tensor_row_major[constraint * input * input + first * input + second] = value;
            }
        }
    }
    validate_symmetric_tensor(
        out.hessian_tensor_row_major,
        out.constraint_count,
        out.variable_count,
        "transformed constraint Hessian tensor"
    );
    return out;
}

LagrangianHessianAssembler::LagrangianHessianAssembler(int variable_count)
    : variable_count_(variable_count) {
    require_positive_count(variable_count_, "Lagrangian Hessian variable count");
}

int LagrangianHessianAssembler::nonzero_count() const {
    return dense_lower_hessian_nonzero_count(variable_count_);
}

NlpHessianStructure LagrangianHessianAssembler::structure() const {
    return dense_lower_hessian_structure(variable_count_);
}

std::vector<double> LagrangianHessianAssembler::values(
    double objective_factor,
    const ObjectiveSecondOrderData& objective,
    const ConstraintSecondOrderData& constraints,
    const std::vector<double>& constraint_multipliers
) const {
    require_finite(objective_factor, "Lagrangian objective factor");
    if (objective.variable_count != variable_count_) {
        throw ValueError("Lagrangian objective Hessian variable count does not match the NLP variable count.");
    }
    validate_symmetric_matrix(objective.hessian_row_major, variable_count_, "Lagrangian objective Hessian");
    require_nonnegative_count(constraints.constraint_count, "Lagrangian constraint count");
    if (constraints.variable_count != variable_count_) {
        throw ValueError("Lagrangian constraint Hessian variable count does not match the NLP variable count.");
    }
    require_size(
        constraint_multipliers,
        static_cast<std::size_t>(constraints.constraint_count),
        "Lagrangian constraint multiplier vector"
    );
    require_finite_values(constraint_multipliers, "Lagrangian constraint multiplier vector");
    validate_symmetric_tensor(
        constraints.hessian_tensor_row_major,
        constraints.constraint_count,
        variable_count_,
        "Lagrangian constraint Hessian tensor"
    );
    const std::vector<bool> active_rows = active_hessian_rows(
        constraints.has_hessian,
        constraints.constraint_count,
        true
    );

    const std::size_t n = static_cast<std::size_t>(variable_count_);
    std::vector<double> dense(n * n, 0.0);
    add_scaled_matrix(dense, objective.hessian_row_major, variable_count_, objective_factor);
    for (std::size_t constraint = 0; constraint < static_cast<std::size_t>(constraints.constraint_count); ++constraint) {
        if (!active_rows[constraint] || constraint_multipliers[constraint] == 0.0) {
            continue;
        }
        const std::size_t offset = constraint * n * n;
        for (std::size_t index = 0; index < n * n; ++index) {
            dense[index] += constraint_multipliers[constraint]
                * constraints.hessian_tensor_row_major[offset + index];
        }
    }
    symmetrize_square_matrix(dense, variable_count_);
    return lower_triangle_values(dense, variable_count_);
}

}  // namespace epcsaft::native::equilibrium_nlp
