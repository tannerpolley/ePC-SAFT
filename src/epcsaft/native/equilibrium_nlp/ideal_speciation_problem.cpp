#include "ideal_speciation_problem.h"

#include "epcsaft_chemical_equilibrium.h"
#include "epcsaft_electrolyte.h"
#include "gibbs_blocks.h"
#include "reaction_block.h"
#include "second_order.h"

#include <cppad/cppad.hpp>
#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>

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

double positive_scale_from_totals(const std::vector<double>& totals) {
    double scale = 0.0;
    for (double value : totals) {
        if (!std::isfinite(value)) {
            throw ValueError("Ideal speciation material totals must be finite.");
        }
        scale += std::abs(value);
    }
    return std::max(1.0, scale);
}

void validate_request(const IdealSpeciationRequest& request) {
    if (request.species_count <= 0) {
        throw ValueError("Ideal Ipopt speciation requires at least one species.");
    }
    if (request.balance_rows <= 0) {
        throw ValueError("Ideal Ipopt speciation requires at least one material balance.");
    }
    if (request.reaction_rows <= 0) {
        throw ValueError("Ideal Ipopt speciation requires at least one reaction.");
    }
    if (!(std::isfinite(request.min_mole_fraction) && request.min_mole_fraction > 0.0)) {
        throw ValueError("Ideal Ipopt speciation requires a positive min_mole_fraction.");
    }
    const auto species = static_cast<std::size_t>(request.species_count);
    require_size(
        request.balance_matrix_row_major,
        static_cast<std::size_t>(request.balance_rows) * species,
        "Ideal Ipopt speciation balance matrix"
    );
    require_size(
        request.total_vector,
        static_cast<std::size_t>(request.balance_rows),
        "Ideal Ipopt speciation total vector"
    );
    require_size(
        request.reaction_stoichiometry_row_major,
        static_cast<std::size_t>(request.reaction_rows) * species,
        "Ideal Ipopt speciation reaction matrix"
    );
    require_size(
        request.log_equilibrium_constants,
        static_cast<std::size_t>(request.reaction_rows),
        "Ideal Ipopt speciation log equilibrium constants"
    );
    require_size(request.initial_x, species, "Ideal Ipopt speciation initial composition");
    if (!request.charges.empty()) {
        require_size(request.charges, species, "Ideal Ipopt speciation charges");
    }
    for (double value : request.balance_matrix_row_major) {
        if (!std::isfinite(value)) {
            throw ValueError("Ideal Ipopt speciation balance matrix must contain finite values.");
        }
    }
    for (double value : request.reaction_stoichiometry_row_major) {
        if (!std::isfinite(value)) {
            throw ValueError("Ideal Ipopt speciation reaction matrix must contain finite values.");
        }
    }
    for (double value : request.log_equilibrium_constants) {
        if (!std::isfinite(value)) {
            throw ValueError("Ideal Ipopt speciation log equilibrium constants must be finite.");
        }
    }
    for (double value : request.charges) {
        if (!std::isfinite(value)) {
            throw ValueError("Ideal Ipopt speciation charges must be finite.");
        }
    }
}

std::vector<double> normalize_initial_x(const std::vector<double>& initial_x, double min_mole_fraction) {
    std::vector<double> out(initial_x.size(), min_mole_fraction);
    double total = 0.0;
    for (std::size_t index = 0; index < initial_x.size(); ++index) {
        const double value = initial_x[index];
        if (!(std::isfinite(value) && value >= -min_mole_fraction)) {
            throw ValueError("Ideal Ipopt speciation initial composition must be finite and non-negative.");
        }
        out[index] = std::max(0.0, value);
        total += out[index];
    }
    if (!(std::isfinite(total) && total > 0.0)) {
        throw ValueError("Ideal Ipopt speciation initial composition must have a positive sum.");
    }
    for (double& value : out) {
        value = std::max(value / total, min_mole_fraction);
    }
    const double clipped_total = std::accumulate(out.begin(), out.end(), 0.0);
    for (double& value : out) {
        value /= clipped_total;
    }
    return out;
}

std::vector<double> canonical_initial_amounts(const IdealSpeciationRequest& request) {
    const std::vector<double> x = normalize_initial_x(request.initial_x, request.min_mole_fraction);
    double numerator = 0.0;
    double denominator = 0.0;
    for (int row = 0; row < request.balance_rows; ++row) {
        double projected = 0.0;
        for (int col = 0; col < request.species_count; ++col) {
            projected += request.balance_matrix_row_major[
                static_cast<std::size_t>(row) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ] * x[static_cast<std::size_t>(col)];
        }
        numerator += projected * request.total_vector[static_cast<std::size_t>(row)];
        denominator += projected * projected;
    }
    const double total_scale = positive_scale_from_totals(request.total_vector);
    double scale = total_scale;
    if (std::isfinite(numerator) && std::isfinite(denominator) && denominator > 0.0) {
        scale = numerator / denominator;
    }
    if (!(std::isfinite(scale) && scale > 0.0)) {
        scale = total_scale;
    }
    const double floor = request.min_mole_fraction * std::max(1.0, scale);
    std::vector<double> amounts(x.size(), floor);
    for (std::size_t index = 0; index < x.size(); ++index) {
        amounts[index] = std::max(floor, x[index] * scale);
    }
    return amounts;
}

bool has_nonzero_charge(const IdealSpeciationRequest& request) {
    for (double charge : request.charges) {
        if (std::abs(charge) > 1.0e-12) {
            return true;
        }
    }
    return false;
}

bool charge_constraint_increases_rank(const IdealSpeciationRequest& request) {
    if (!has_nonzero_charge(request)) {
        return false;
    }
    Eigen::MatrixXd balances(request.balance_rows, request.species_count);
    for (int row = 0; row < request.balance_rows; ++row) {
        for (int col = 0; col < request.species_count; ++col) {
            balances(row, col) = request.balance_matrix_row_major[
                static_cast<std::size_t>(row) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ];
        }
    }
    Eigen::MatrixXd with_charge(request.balance_rows + 1, request.species_count);
    with_charge.topRows(request.balance_rows) = balances;
    for (int col = 0; col < request.species_count; ++col) {
        with_charge(request.balance_rows, col) = request.charges[static_cast<std::size_t>(col)];
    }
    const Eigen::FullPivLU<Eigen::MatrixXd> base_rank(balances);
    const Eigen::FullPivLU<Eigen::MatrixXd> charged_rank(with_charge);
    return charged_rank.rank() > base_rank.rank();
}

std::string ideal_derivative_backend_from_options(const std::string& requested) {
    if (requested == "cppad") {
        return "cppad";
    }
    if (requested == "auto" || requested == "analytic") {
        return "analytic";
    }
    throw ValueError("Ideal Ipopt speciation jacobian_backend must be 'auto', 'analytic', or 'cppad'.");
}

std::vector<double> cppad_ideal_reduced_gibbs_gradient(
    const std::vector<double>& amounts,
    const std::vector<double>& standard_mu_rt
) {
    using CppADScalar = CppAD::AD<double>;
    if (amounts.empty() || amounts.size() != standard_mu_rt.size()) {
        throw ValueError("CppAD ideal Gibbs gradient requires one standard chemical potential per species.");
    }
    const std::size_t species = amounts.size();
    std::vector<CppADScalar> variables(species);
    std::vector<double> point(species);
    for (std::size_t index = 0; index < species; ++index) {
        if (!(std::isfinite(amounts[index]) && amounts[index] > 0.0) || !std::isfinite(standard_mu_rt[index])) {
            throw ValueError("CppAD ideal Gibbs gradient requires positive amounts and finite standard potentials.");
        }
        variables[index] = amounts[index];
        point[index] = amounts[index];
    }
    CppAD::Independent(variables);

    CppADScalar total = CppADScalar(0.0);
    for (const CppADScalar& amount : variables) {
        total += amount;
    }
    CppADScalar objective = CppADScalar(0.0);
    for (std::size_t index = 0; index < species; ++index) {
        const CppADScalar mole_fraction = variables[index] / total;
        objective += variables[index] * (CppAD::log(mole_fraction) + standard_mu_rt[index]);
    }
    std::vector<CppADScalar> outputs = {objective};
    CppAD::ADFun<double> function(variables, outputs);
    std::vector<double> gradient = function.Jacobian(point);
    if (gradient.size() != species) {
        throw ValueError("CppAD ideal Gibbs gradient shape did not match the species count.");
    }
    for (double value : gradient) {
        if (!std::isfinite(value)) {
            throw ValueError("CppAD ideal Gibbs gradient produced a non-finite value.");
        }
    }
    return gradient;
}

class IdealSpeciationProblem final : public NlpProblem {
public:
    explicit IdealSpeciationProblem(IdealSpeciationRequest request, std::string derivative_backend)
        : request_(std::move(request)),
          derivative_backend_(std::move(derivative_backend)),
          standard_mu_rt_(standard_mu_rt_from_reactions(
              request_.species_count,
              request_.reaction_rows,
              request_.reaction_stoichiometry_row_major,
              request_.log_equilibrium_constants
          )),
          initial_amounts_(canonical_initial_amounts(request_)) {
        total_scale_ = positive_scale_from_totals(request_.total_vector);
        include_charge_constraint_ = charge_constraint_increases_rank(request_);
    }

    std::string name() const override {
        return "ideal_homogeneous_reactive_speciation";
    }

    int variable_count() const override {
        return request_.species_count;
    }

    int constraint_count() const override {
        return request_.balance_rows + (include_charge_constraint_ ? 1 : 0);
    }

    int jacobian_nonzero_count() const override {
        int balance_nonzeros = 0;
        for (double coefficient : request_.balance_matrix_row_major) {
            if (coefficient != 0.0) {
                ++balance_nonzeros;
            }
        }
        int charge_nonzeros = 0;
        if (include_charge_constraint_) {
            for (double charge : request_.charges) {
                if (charge != 0.0) {
                    ++charge_nonzeros;
                }
            }
        }
        return balance_nonzeros + charge_nonzeros;
    }

    NlpBounds bounds() const override {
        NlpBounds out;
        const double lower = request_.min_mole_fraction * total_scale_;
        const double upper = std::max(1.0, 10.0 * total_scale_);
        out.variable_lower.assign(static_cast<std::size_t>(request_.species_count), lower);
        out.variable_upper.assign(static_cast<std::size_t>(request_.species_count), upper);
        out.constraint_lower = request_.total_vector;
        out.constraint_upper = request_.total_vector;
        if (include_charge_constraint_) {
            out.constraint_lower.push_back(0.0);
            out.constraint_upper.push_back(0.0);
        }
        return out;
    }

    std::vector<double> initial_point() const override {
        return initial_amounts_;
    }

    double objective(const std::vector<double>& variables) const override {
        return evaluate_ideal_reduced_gibbs(variables, standard_mu_rt_, false).value;
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        if (derivative_backend_ == "cppad") {
            return cppad_ideal_reduced_gibbs_gradient(variables, standard_mu_rt_);
        }
        return evaluate_ideal_reduced_gibbs(variables, standard_mu_rt_, false).gradient;
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        require_size(variables, static_cast<std::size_t>(request_.species_count), "Ideal Ipopt variables");
        std::vector<double> out(static_cast<std::size_t>(constraint_count()), 0.0);
        for (int row = 0; row < request_.balance_rows; ++row) {
            for (int col = 0; col < request_.species_count; ++col) {
                out[static_cast<std::size_t>(row)] += request_.balance_matrix_row_major[
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(request_.species_count)
                    + static_cast<std::size_t>(col)
                ] * variables[static_cast<std::size_t>(col)];
            }
        }
        if (include_charge_constraint_) {
            for (int col = 0; col < request_.species_count; ++col) {
                out[static_cast<std::size_t>(request_.balance_rows)] +=
                    request_.charges[static_cast<std::size_t>(col)] * variables[static_cast<std::size_t>(col)];
            }
        }
        return out;
    }

    NlpJacobianStructure jacobian_structure() const override {
        NlpJacobianStructure out;
        out.rows.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        out.cols.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < request_.balance_rows; ++row) {
            for (int col = 0; col < request_.species_count; ++col) {
                const double coefficient = request_.balance_matrix_row_major[
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(request_.species_count)
                    + static_cast<std::size_t>(col)
                ];
                if (coefficient == 0.0) {
                    continue;
                }
                out.rows.push_back(row);
                out.cols.push_back(col);
            }
        }
        if (include_charge_constraint_) {
            for (int col = 0; col < request_.species_count; ++col) {
                if (request_.charges[static_cast<std::size_t>(col)] == 0.0) {
                    continue;
                }
                out.rows.push_back(request_.balance_rows);
                out.cols.push_back(col);
            }
        }
        return out;
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        (void)variables;
        std::vector<double> out;
        out.reserve(static_cast<std::size_t>(jacobian_nonzero_count()));
        for (int row = 0; row < request_.balance_rows; ++row) {
            for (int col = 0; col < request_.species_count; ++col) {
                const double coefficient = request_.balance_matrix_row_major[
                    static_cast<std::size_t>(row) * static_cast<std::size_t>(request_.species_count)
                    + static_cast<std::size_t>(col)
                ];
                if (coefficient == 0.0) {
                    continue;
                }
                out.push_back(coefficient);
            }
        }
        if (include_charge_constraint_) {
            for (double charge : request_.charges) {
                if (charge != 0.0) {
                    out.push_back(charge);
                }
            }
        }
        return out;
    }

    bool has_exact_hessian() const override {
        return true;
    }

    int hessian_nonzero_count() const override {
        return LagrangianHessianAssembler(request_.species_count).nonzero_count();
    }

    NlpHessianStructure hessian_structure() const override {
        return LagrangianHessianAssembler(request_.species_count).structure();
    }

    std::vector<double> hessian_values(
        const std::vector<double>& variables,
        double objective_factor,
        const std::vector<double>& constraint_multipliers
    ) const override {
        const IdealReducedGibbsResult gibbs = evaluate_ideal_reduced_gibbs(variables, standard_mu_rt_, true);
        const std::size_t n = static_cast<std::size_t>(request_.species_count);
        if (gibbs.hessian_row_major.size() != n * n) {
            throw ValueError("Ideal Ipopt speciation Hessian shape did not match the species count.");
        }
        ObjectiveSecondOrderData objective;
        objective.variable_count = request_.species_count;
        objective.hessian_row_major = gibbs.hessian_row_major;
        objective.backend = "analytic";

        ConstraintSecondOrderData constraints;
        constraints.constraint_count = constraint_count();
        constraints.variable_count = request_.species_count;
        constraints.hessian_tensor_row_major.assign(
            static_cast<std::size_t>(constraints.constraint_count)
                * static_cast<std::size_t>(request_.species_count)
                * static_cast<std::size_t>(request_.species_count),
            0.0
        );
        constraints.has_hessian.assign(static_cast<std::size_t>(constraints.constraint_count), false);
        constraints.backend = "linear";
        return LagrangianHessianAssembler(request_.species_count).values(
            objective_factor,
            objective,
            constraints,
            constraint_multipliers
        );
    }

    std::string hessian_backend() const override {
        return "analytic";
    }

    NlpScaling scaling() const override {
        NlpScaling out;
        out.objective = 1.0 / std::max(1.0, total_scale_);
        out.variables.assign(static_cast<std::size_t>(request_.species_count), 1.0 / total_scale_);
        out.constraints.reserve(static_cast<std::size_t>(request_.balance_rows));
        for (double total : request_.total_vector) {
            out.constraints.push_back(1.0 / std::max(1.0, std::abs(total)));
        }
        if (include_charge_constraint_) {
            out.constraints.push_back(1.0);
        }
        return out;
    }

    const std::vector<double>& standard_mu_rt() const {
        return standard_mu_rt_;
    }

private:
    IdealSpeciationRequest request_;
    std::string derivative_backend_;
    std::vector<double> standard_mu_rt_;
    std::vector<double> initial_amounts_;
    double total_scale_ = 1.0;
    bool include_charge_constraint_ = false;
};

std::vector<double> normalized_composition_from_amounts(const std::vector<double>& amounts) {
    const double total = std::accumulate(amounts.begin(), amounts.end(), 0.0);
    if (!(std::isfinite(total) && total > 0.0)) {
        throw SolutionError("Ideal Ipopt speciation produced non-positive total amount.");
    }
    std::vector<double> composition(amounts.size(), 0.0);
    for (std::size_t index = 0; index < amounts.size(); ++index) {
        composition[index] = amounts[index] / total;
    }
    return composition;
}

std::vector<double> mass_balance_residuals(
    const IdealSpeciationRequest& request,
    const std::vector<double>& amounts
) {
    std::vector<double> residuals(static_cast<std::size_t>(request.balance_rows), 0.0);
    for (int row = 0; row < request.balance_rows; ++row) {
        for (int col = 0; col < request.species_count; ++col) {
            residuals[static_cast<std::size_t>(row)] += request.balance_matrix_row_major[
                static_cast<std::size_t>(row) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ] * amounts[static_cast<std::size_t>(col)];
        }
        residuals[static_cast<std::size_t>(row)] -= request.total_vector[static_cast<std::size_t>(row)];
    }
    return residuals;
}

double charge_residual(const IdealSpeciationRequest& request, const std::vector<double>& amounts) {
    if (request.charges.empty()) {
        return 0.0;
    }
    double residual = 0.0;
    for (int index = 0; index < request.species_count; ++index) {
        residual += request.charges[static_cast<std::size_t>(index)] * amounts[static_cast<std::size_t>(index)];
    }
    return residual;
}

double max_abs(const std::vector<double>& values) {
    double out = 0.0;
    for (double value : values) {
        out = std::max(out, std::abs(value));
    }
    return out;
}

std::vector<double> dense_matrix_to_row_major(const Eigen::MatrixXd& matrix) {
    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(matrix.rows() * matrix.cols()));
    for (Eigen::Index row = 0; row < matrix.rows(); ++row) {
        for (Eigen::Index col = 0; col < matrix.cols(); ++col) {
            out.push_back(matrix(row, col));
        }
    }
    return out;
}

Eigen::MatrixXd ideal_log_amount_residual_jacobian(
    const IdealSpeciationRequest& request,
    const std::vector<double>& amounts,
    const std::vector<double>& composition
) {
    const int rows = request.balance_rows + 1 + request.reaction_rows;
    Eigen::MatrixXd jac = Eigen::MatrixXd::Zero(rows, request.species_count);
    for (int row = 0; row < request.balance_rows; ++row) {
        for (int col = 0; col < request.species_count; ++col) {
            jac(row, col) = request.balance_matrix_row_major[
                static_cast<std::size_t>(row) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ] * amounts[static_cast<std::size_t>(col)];
        }
    }
    const int charge_row = request.balance_rows;
    if (!request.charges.empty()) {
        for (int col = 0; col < request.species_count; ++col) {
            jac(charge_row, col) =
                request.charges[static_cast<std::size_t>(col)] * amounts[static_cast<std::size_t>(col)];
        }
    }
    const int reaction_offset = request.balance_rows + 1;
    for (int reaction = 0; reaction < request.reaction_rows; ++reaction) {
        double stoich_sum = 0.0;
        for (int col = 0; col < request.species_count; ++col) {
            stoich_sum += request.reaction_stoichiometry_row_major[
                static_cast<std::size_t>(reaction) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ];
        }
        for (int col = 0; col < request.species_count; ++col) {
            const double coefficient = request.reaction_stoichiometry_row_major[
                static_cast<std::size_t>(reaction) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ];
            jac(reaction_offset + reaction, col) =
                coefficient - stoich_sum * composition[static_cast<std::size_t>(col)];
        }
    }
    return jac;
}

Eigen::MatrixXd cppad_ideal_log_amount_residual_jacobian(
    const IdealSpeciationRequest& request,
    const std::vector<double>& amounts
) {
    using CppADScalar = CppAD::AD<double>;
    require_size(amounts, static_cast<std::size_t>(request.species_count), "Ideal CppAD residual amounts");
    const int rows = request.balance_rows + 1 + request.reaction_rows;
    std::vector<CppADScalar> variables(static_cast<std::size_t>(request.species_count));
    std::vector<double> point(static_cast<std::size_t>(request.species_count));
    for (int index = 0; index < request.species_count; ++index) {
        const double amount = amounts[static_cast<std::size_t>(index)];
        if (!(std::isfinite(amount) && amount > 0.0)) {
            throw ValueError("Ideal CppAD residual Jacobian requires positive species amounts.");
        }
        const double log_amount = std::log(amount);
        variables[static_cast<std::size_t>(index)] = log_amount;
        point[static_cast<std::size_t>(index)] = log_amount;
    }
    CppAD::Independent(variables);

    std::vector<CppADScalar> n(static_cast<std::size_t>(request.species_count));
    CppADScalar total = CppADScalar(0.0);
    for (int col = 0; col < request.species_count; ++col) {
        n[static_cast<std::size_t>(col)] = CppAD::exp(variables[static_cast<std::size_t>(col)]);
        total += n[static_cast<std::size_t>(col)];
    }
    std::vector<CppADScalar> x(static_cast<std::size_t>(request.species_count));
    for (int col = 0; col < request.species_count; ++col) {
        x[static_cast<std::size_t>(col)] = n[static_cast<std::size_t>(col)] / total;
    }

    std::vector<CppADScalar> outputs(static_cast<std::size_t>(rows), CppADScalar(0.0));
    for (int row = 0; row < request.balance_rows; ++row) {
        CppADScalar residual = CppADScalar(0.0);
        for (int col = 0; col < request.species_count; ++col) {
            residual += request.balance_matrix_row_major[
                static_cast<std::size_t>(row) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ] * n[static_cast<std::size_t>(col)];
        }
        outputs[static_cast<std::size_t>(row)] = residual - request.total_vector[static_cast<std::size_t>(row)];
    }

    CppADScalar charge = CppADScalar(0.0);
    if (request.charges.size() == static_cast<std::size_t>(request.species_count)) {
        for (int col = 0; col < request.species_count; ++col) {
            charge += request.charges[static_cast<std::size_t>(col)] * n[static_cast<std::size_t>(col)];
        }
    }
    outputs[static_cast<std::size_t>(request.balance_rows)] = charge;

    const int reaction_offset = request.balance_rows + 1;
    for (int reaction = 0; reaction < request.reaction_rows; ++reaction) {
        CppADScalar residual = -request.log_equilibrium_constants[static_cast<std::size_t>(reaction)];
        for (int col = 0; col < request.species_count; ++col) {
            const double coefficient = request.reaction_stoichiometry_row_major[
                static_cast<std::size_t>(reaction) * static_cast<std::size_t>(request.species_count)
                + static_cast<std::size_t>(col)
            ];
            residual += coefficient * CppAD::log(x[static_cast<std::size_t>(col)]);
        }
        outputs[static_cast<std::size_t>(reaction_offset + reaction)] = residual;
    }

    CppAD::ADFun<double> function(variables, outputs);
    std::vector<double> jacobian = function.Jacobian(point);
    if (jacobian.size() != static_cast<std::size_t>(rows * request.species_count)) {
        throw ValueError("Ideal CppAD residual Jacobian shape did not match the residual model.");
    }
    Eigen::MatrixXd out(rows, request.species_count);
    for (int row = 0; row < rows; ++row) {
        for (int col = 0; col < request.species_count; ++col) {
            const double value =
                jacobian[static_cast<std::size_t>(row * request.species_count + col)];
            if (!std::isfinite(value)) {
                throw ValueError("Ideal CppAD residual Jacobian produced a non-finite value.");
            }
            out(row, col) = value;
        }
    }
    return out;
}

void add_ideal_implicit_sensitivity_diagnostics(
    ChemicalEquilibriumResultNative& result,
    const IdealSpeciationRequest& request,
    const std::vector<double>& amounts,
    const std::vector<double>& composition,
    const std::vector<double>& residuals,
    const std::string& derivative_backend
) {
    Eigen::MatrixXd residual_state = derivative_backend == "cppad"
        ? cppad_ideal_log_amount_residual_jacobian(request, amounts)
        : ideal_log_amount_residual_jacobian(request, amounts, composition);
    Eigen::MatrixXd residual_parameter =
        Eigen::MatrixXd::Zero(residual_state.rows(), request.reaction_rows);
    const int reaction_offset = request.balance_rows + 1;
    for (int reaction = 0; reaction < request.reaction_rows; ++reaction) {
        residual_parameter(reaction_offset + reaction, reaction) = -1.0;
    }
    Eigen::MatrixXd sensitivity = residual_state.colPivHouseholderQr().solve(-residual_parameter);

    std::vector<double> log_amounts(amounts.size(), 0.0);
    for (std::size_t index = 0; index < amounts.size(); ++index) {
        log_amounts[index] = std::log(amounts[index]);
    }
    result.diagnostics_string["implicit_sensitivity_backend"] =
        derivative_backend == "cppad" ? "cppad_implicit" : "analytic_implicit";
    result.diagnostics_string["reactive_speciation_sensitivity_parameter"] = "log_equilibrium_constants";
    result.diagnostics_int["reactive_speciation_residual_rows"] = static_cast<int>(residual_state.rows());
    result.diagnostics_int["reactive_speciation_state_size"] = static_cast<int>(residual_state.cols());
    result.diagnostics_int["reactive_speciation_parameter_size"] = request.reaction_rows;
    result.diagnostics_vector["reactive_speciation_state"] = log_amounts;
    result.diagnostics_vector["reactive_speciation_residual"] = residuals;
    result.diagnostics_vector["reactive_speciation_residual_state_jacobian_row_major"] =
        dense_matrix_to_row_major(residual_state);
    result.diagnostics_vector["reactive_speciation_residual_parameter_jacobian_row_major"] =
        dense_matrix_to_row_major(residual_parameter);
    result.diagnostics_vector["reactive_speciation_log_amount_sensitivity_to_log_k_row_major"] =
        dense_matrix_to_row_major(sensitivity);
}

}  // namespace

IdealSpeciationIpoptResult solve_ideal_speciation_ipopt(
    const IdealSpeciationRequest& request,
    const IpoptSolveOptions& options,
    const std::string& derivative_backend
) {
    validate_request(request);
    IdealSpeciationProblem problem(request, derivative_backend);
    if (options.hessian_mode == "exact" && !problem.has_exact_hessian()) {
        throw ValueError("Ideal Ipopt speciation exact Hessian provider is unavailable.");
    }
    IpoptSolveResult ipopt = solve_ipopt_nlp(problem, options);
    if (!ipopt.accepted) {
        throw SolutionError("Ipopt did not accept the ideal reactive speciation NLP solution.");
    }
    if (ipopt.variables.size() != static_cast<std::size_t>(request.species_count)) {
        throw SolutionError("Ipopt returned an invalid ideal reactive speciation variable vector.");
    }

    IdealSpeciationIpoptResult out;
    out.ipopt = std::move(ipopt);
    out.amounts = out.ipopt.variables;
    out.composition = normalized_composition_from_amounts(out.amounts);
    out.activity_coefficients.assign(static_cast<std::size_t>(request.species_count), 1.0);
    out.mass_balance_residuals = mass_balance_residuals(request, out.amounts);
    out.charge_residual = charge_residual(request, out.amounts);
    out.reaction_residuals = evaluate_ideal_reaction_quotients(
        out.amounts,
        request.reaction_rows,
        request.reaction_stoichiometry_row_major,
        request.log_equilibrium_constants
    ).residuals;
    out.standard_mu_rt = problem.standard_mu_rt();
    return out;
}

ChemicalEquilibriumResultNative solve_ideal_speciation_chemical_equilibrium_ipopt(
    const IdealSpeciationRequest& request,
    const ChemicalEquilibriumOptionsNative& options
) {
    IpoptSolveOptions solve_options;
    solve_options.max_iterations = options.max_iterations;
    solve_options.tolerance = options.tolerance;
    solve_options.hessian_mode = options.hessian_mode;
    solve_options.iteration_history_limit = options.iteration_history_limit;
    solve_options.linear_solver = options.linear_solver;
    solve_options.acceptable_tolerance = options.acceptable_tolerance;
    solve_options.constraint_violation_tolerance = options.constraint_violation_tolerance;
    solve_options.dual_infeasibility_tolerance = options.dual_infeasibility_tolerance;
    solve_options.complementarity_tolerance = options.complementarity_tolerance;
    solve_options.initial_variables = options.initial_variables;
    solve_options.initial_bound_lower_multipliers = options.initial_bound_lower_multipliers;
    solve_options.initial_bound_upper_multipliers = options.initial_bound_upper_multipliers;
    solve_options.initial_constraint_multipliers = options.initial_constraint_multipliers;
    const std::string derivative_backend = ideal_derivative_backend_from_options(options.jacobian_backend);
    const IdealSpeciationIpoptResult ipopt_result =
        solve_ideal_speciation_ipopt(request, solve_options, derivative_backend);

    std::vector<double> residuals = ipopt_result.mass_balance_residuals;
    residuals.push_back(ipopt_result.charge_residual);
    residuals.insert(residuals.end(), ipopt_result.reaction_residuals.begin(), ipopt_result.reaction_residuals.end());
    const double residual_norm = max_abs(residuals);

    ChemicalEquilibriumResultNative result;
    result.success = ipopt_result.ipopt.accepted && residual_norm <= options.tolerance;
    if (result.success) {
        result.message = "converged";
    } else if (!ipopt_result.ipopt.accepted) {
        result.message = "Ipopt did not accept the ideal reactive speciation NLP solution.";
    } else {
        result.message = "Ipopt ideal reactive speciation residual acceptance gate failed";
    }
    result.composition = ipopt_result.composition;
    if (options.activity_output == "always") {
        result.activity_coefficients = ipopt_result.activity_coefficients;
    }
    result.mass_balance_residuals = ipopt_result.mass_balance_residuals;
    result.charge_residual = ipopt_result.charge_residual;
    result.reaction_residuals = ipopt_result.reaction_residuals;
    result.diagnostics_double = ipopt_result.ipopt.diagnostics_double;
    result.diagnostics_int = ipopt_result.ipopt.diagnostics_int;
    result.diagnostics_bool = ipopt_result.ipopt.diagnostics_bool;
    result.diagnostics_string = ipopt_result.ipopt.diagnostics_string;
    result.continuation_variables = ipopt_result.ipopt.variables;
    result.continuation_bound_lower_multipliers = ipopt_result.ipopt.bound_lower_multipliers;
    result.continuation_bound_upper_multipliers = ipopt_result.ipopt.bound_upper_multipliers;
    result.continuation_constraint_multipliers = ipopt_result.ipopt.constraint_multipliers;
    result.iteration_history = ipopt_result.ipopt.iteration_history;
    result.diagnostics_string["solver_language"] = "c++";
    result.diagnostics_string["native_entrypoint"] = "_solve_chemical_equilibrium_native";
    result.diagnostics_string["problem_class"] = "homogeneous_ideal_gibbs_speciation";
    result.diagnostics_string["activity_model"] = "ideal";
    result.diagnostics_string["activity_output"] = options.activity_output;
    result.diagnostics_string["activity_basis"] = "ideal_mole_fraction";
    result.diagnostics_string["phase"] = options.phase;
    result.diagnostics_string["requested_jacobian_backend"] = options.jacobian_backend;
    result.diagnostics_string["jacobian_backend"] = derivative_backend;
    result.diagnostics_string["derivative_backend"] = derivative_backend;
    result.diagnostics_string["derivative_capability_path"] = derivative_backend == "cppad"
        ? "chemical_equilibrium:ideal_mole_fraction:ipopt_amount_gibbs_cppad"
        : "chemical_equilibrium:ideal_mole_fraction:ipopt_amount_gibbs";
    result.diagnostics_string["selected_solver_backend"] = "native_ipopt";
    result.diagnostics_string["solver_selection_reason"] = "explicit_request";
    result.diagnostics_string["ipopt_solver_status"] = ipopt_result.ipopt.solver_status;
    result.diagnostics_string["ipopt_application_status"] = ipopt_result.ipopt.application_status;
    result.diagnostics_bool["derivative_available"] = true;
    result.diagnostics_bool["jacobian_available"] = true;
    result.diagnostics_bool["activity_coefficients_evaluated"] = !result.activity_coefficients.empty();
    result.diagnostics_bool["charge_constraint_in_nlp"] = charge_constraint_increases_rank(request);
    result.diagnostics_bool["ipopt_solver_ran"] = ipopt_result.ipopt.solver_ran;
    result.diagnostics_bool["ipopt_accepted"] = ipopt_result.ipopt.accepted;
    result.diagnostics_double["residual_norm"] = residual_norm;
    result.diagnostics_double["tolerance"] = options.tolerance;
    result.diagnostics_double["objective"] = ipopt_result.ipopt.objective;
    result.diagnostics_vector["history"] = {};
    result.diagnostics_vector["phase_handoff_composition"] = ipopt_result.composition;
    result.diagnostics_vector["ideal_gibbs_standard_mu_rt"] = ipopt_result.standard_mu_rt;
    add_ideal_implicit_sensitivity_diagnostics(
        result,
        request,
        ipopt_result.amounts,
        ipopt_result.composition,
        residuals,
        derivative_backend
    );
    return result;
}

}  // namespace epcsaft::native::equilibrium_nlp
