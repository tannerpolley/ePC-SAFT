#include "ipopt_adapter.h"

#include "epcsaft_electrolyte.h"

#include <algorithm>
#include <cmath>
#include <sstream>

#ifdef EPCSAFT_HAS_IPOPT
#if __has_include(<coin-or/IpIpoptApplication.hpp>)
#include <coin-or/IpIpoptApplication.hpp>
#include <coin-or/IpTNLP.hpp>
#elif __has_include(<IpIpoptApplication.hpp>)
#include <IpIpoptApplication.hpp>
#include <IpTNLP.hpp>
#else
#error "EPCSAFT_HAS_IPOPT is enabled, but Ipopt C++ headers were not found."
#endif
#endif

namespace epcsaft::native::equilibrium_nlp {

namespace {

class QuadraticSmokeProblem final : public NlpProblem {
public:
    std::string name() const override {
        return "quadratic_linear_constraint_smoke";
    }

    int variable_count() const override {
        return 2;
    }

    int constraint_count() const override {
        return 1;
    }

    int jacobian_nonzero_count() const override {
        return 2;
    }

    NlpBounds bounds() const override {
        return {
            {-10.0, -10.0},
            {10.0, 10.0},
            {3.0},
            {3.0},
        };
    }

    std::vector<double> initial_point() const override {
        return {0.5, 2.5};
    }

    double objective(const std::vector<double>& variables) const override {
        const double dx = variables[0] - 1.0;
        const double dy = variables[1] - 2.0;
        return dx * dx + dy * dy;
    }

    std::vector<double> objective_gradient(const std::vector<double>& variables) const override {
        return {2.0 * (variables[0] - 1.0), 2.0 * (variables[1] - 2.0)};
    }

    std::vector<double> constraints(const std::vector<double>& variables) const override {
        return {variables[0] + variables[1]};
    }

    NlpJacobianStructure jacobian_structure() const override {
        return {{0, 0}, {0, 1}};
    }

    std::vector<double> jacobian_values(const std::vector<double>& variables) const override {
        (void)variables;
        return {1.0, 1.0};
    }

    NlpScaling scaling() const override {
        return {1.0, {1.0, 1.0}, {1.0}};
    }

    std::map<std::string, std::string> diagnostics() const override {
        return {
            {"smoke_problem", "quadratic_linear_constraint"},
            {"gradient_backend", "analytic"},
            {"jacobian_backend", "analytic"},
        };
    }
};

#ifdef EPCSAFT_HAS_IPOPT
bool all_finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

std::vector<double> vector_from_raw(const Ipopt::Number* values, Ipopt::Index count) {
    return std::vector<double>(values, values + count);
}

std::string solver_status_name(Ipopt::SolverReturn status) {
    switch (status) {
        case Ipopt::SUCCESS:
            return "success";
        case Ipopt::MAXITER_EXCEEDED:
            return "max_iterations_exceeded";
        case Ipopt::STOP_AT_ACCEPTABLE_POINT:
            return "acceptable_point";
        case Ipopt::LOCAL_INFEASIBILITY:
            return "local_infeasibility";
        case Ipopt::USER_REQUESTED_STOP:
            return "user_requested_stop";
        case Ipopt::DIVERGING_ITERATES:
            return "diverging_iterates";
        case Ipopt::RESTORATION_FAILURE:
            return "restoration_failure";
        case Ipopt::ERROR_IN_STEP_COMPUTATION:
            return "error_in_step_computation";
        case Ipopt::INVALID_NUMBER_DETECTED:
            return "invalid_number_detected";
        default:
            return "ipopt_status_" + std::to_string(static_cast<int>(status));
    }
}

std::string application_status_name(Ipopt::ApplicationReturnStatus status) {
    switch (status) {
        case Ipopt::Solve_Succeeded:
            return "solve_succeeded";
        case Ipopt::Solved_To_Acceptable_Level:
            return "solved_to_acceptable_level";
        case Ipopt::Maximum_Iterations_Exceeded:
            return "maximum_iterations_exceeded";
        case Ipopt::Infeasible_Problem_Detected:
            return "infeasible_problem_detected";
        case Ipopt::Invalid_Problem_Definition:
            return "invalid_problem_definition";
        case Ipopt::Invalid_Option:
            return "invalid_option";
        case Ipopt::Invalid_Number_Detected:
            return "invalid_number_detected";
        case Ipopt::Unrecoverable_Exception:
            return "unrecoverable_exception";
        default:
            return "ipopt_application_status_" + std::to_string(static_cast<int>(status));
    }
}

class IpoptTnlpAdapter final : public Ipopt::TNLP {
public:
    explicit IpoptTnlpAdapter(const NlpProblem& problem, const IpoptSolveOptions& options)
        : problem_(problem),
          options_(options),
          bounds_(problem_.bounds()),
          initial_(problem_.initial_point()),
          jacobian_structure_(problem_.jacobian_structure()),
          scaling_(problem_.scaling()) {
        result_.variables = initial_;
        result_.hessian_strategy = options_.limited_memory_hessian ? "limited_memory" : "exact";
        result_.diagnostics_string["problem_name"] = problem_.name();
        for (const auto& item : problem_.diagnostics()) {
            result_.diagnostics_string[item.first] = item.second;
        }
    }

    bool get_nlp_info(
        Ipopt::Index& n,
        Ipopt::Index& m,
        Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag,
        IndexStyleEnum& index_style
    ) override {
        n = problem_.variable_count();
        m = problem_.constraint_count();
        nnz_jac_g = problem_.jacobian_nonzero_count();
        nnz_h_lag = 0;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(
        Ipopt::Index n,
        Ipopt::Number* x_l,
        Ipopt::Number* x_u,
        Ipopt::Index m,
        Ipopt::Number* g_l,
        Ipopt::Number* g_u
    ) override {
        if (n != problem_.variable_count() || m != problem_.constraint_count()) {
            return false;
        }
        std::copy(bounds_.variable_lower.begin(), bounds_.variable_lower.end(), x_l);
        std::copy(bounds_.variable_upper.begin(), bounds_.variable_upper.end(), x_u);
        std::copy(bounds_.constraint_lower.begin(), bounds_.constraint_lower.end(), g_l);
        std::copy(bounds_.constraint_upper.begin(), bounds_.constraint_upper.end(), g_u);
        return all_finite(bounds_.variable_lower) && all_finite(bounds_.variable_upper)
            && all_finite(bounds_.constraint_lower) && all_finite(bounds_.constraint_upper);
    }

    bool get_starting_point(
        Ipopt::Index n,
        bool init_x,
        Ipopt::Number* x,
        bool init_z,
        Ipopt::Number* z_L,
        Ipopt::Number* z_U,
        Ipopt::Index m,
        bool init_lambda,
        Ipopt::Number* lambda
    ) override {
        (void)z_L;
        (void)z_U;
        (void)m;
        (void)lambda;
        if (!init_x || init_z || init_lambda || n != problem_.variable_count()) {
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        return all_finite(initial_);
    }

    bool get_scaling_parameters(
        Ipopt::Number& obj_scaling,
        bool& use_x_scaling,
        Ipopt::Index n,
        Ipopt::Number* x_scaling,
        bool& use_g_scaling,
        Ipopt::Index m,
        Ipopt::Number* g_scaling
    ) override {
        obj_scaling = scaling_.objective;
        use_x_scaling = !scaling_.variables.empty();
        use_g_scaling = !scaling_.constraints.empty();
        if (use_x_scaling) {
            if (n != problem_.variable_count()) {
                return false;
            }
            std::copy(scaling_.variables.begin(), scaling_.variables.end(), x_scaling);
        }
        if (use_g_scaling) {
            if (m != problem_.constraint_count()) {
                return false;
            }
            std::copy(scaling_.constraints.begin(), scaling_.constraints.end(), g_scaling);
        }
        return std::isfinite(obj_scaling) && obj_scaling > 0.0 && all_finite(scaling_.variables)
            && all_finite(scaling_.constraints);
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number& obj_value
    ) override {
        (void)new_x;
        if (n != problem_.variable_count()) {
            return false;
        }
        obj_value = problem_.objective(vector_from_raw(x, n));
        return std::isfinite(obj_value);
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number* grad_f
    ) override {
        (void)new_x;
        if (n != problem_.variable_count()) {
            return false;
        }
        const std::vector<double> gradient = problem_.objective_gradient(vector_from_raw(x, n));
        if (gradient.size() != static_cast<std::size_t>(n) || !all_finite(gradient)) {
            return false;
        }
        std::copy(gradient.begin(), gradient.end(), grad_f);
        return true;
    }

    bool eval_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Index m,
        Ipopt::Number* g
    ) override {
        (void)new_x;
        if (n != problem_.variable_count() || m != problem_.constraint_count()) {
            return false;
        }
        const std::vector<double> values = problem_.constraints(vector_from_raw(x, n));
        if (values.size() != static_cast<std::size_t>(m) || !all_finite(values)) {
            return false;
        }
        std::copy(values.begin(), values.end(), g);
        return true;
    }

    bool eval_jac_g(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Index m,
        Ipopt::Index nele_jac,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values
    ) override {
        (void)new_x;
        if (n != problem_.variable_count() || m != problem_.constraint_count()
            || nele_jac != problem_.jacobian_nonzero_count()) {
            return false;
        }
        if (values == nullptr) {
            for (std::size_t index = 0; index < jacobian_structure_.rows.size(); ++index) {
                iRow[index] = jacobian_structure_.rows[index];
                jCol[index] = jacobian_structure_.cols[index];
            }
            return true;
        }
        const std::vector<double> jacobian = problem_.jacobian_values(vector_from_raw(x, n));
        if (jacobian.size() != static_cast<std::size_t>(nele_jac) || !all_finite(jacobian)) {
            return false;
        }
        std::copy(jacobian.begin(), jacobian.end(), values);
        return true;
    }

    bool eval_h(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number obj_factor,
        Ipopt::Index m,
        const Ipopt::Number* lambda,
        bool new_lambda,
        Ipopt::Index nele_hess,
        Ipopt::Index* iRow,
        Ipopt::Index* jCol,
        Ipopt::Number* values
    ) override {
        (void)n;
        (void)x;
        (void)new_x;
        (void)obj_factor;
        (void)m;
        (void)lambda;
        (void)new_lambda;
        (void)nele_hess;
        (void)iRow;
        (void)jCol;
        (void)values;
        return false;
    }

    void finalize_solution(
        Ipopt::SolverReturn status,
        Ipopt::Index n,
        const Ipopt::Number* x,
        const Ipopt::Number* z_L,
        const Ipopt::Number* z_U,
        Ipopt::Index m,
        const Ipopt::Number* g,
        const Ipopt::Number* lambda,
        Ipopt::Number obj_value,
        const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq
    ) override {
        (void)z_L;
        (void)z_U;
        (void)lambda;
        (void)ip_data;
        (void)ip_cq;
        result_.solver_ran = true;
        result_.solver_status = solver_status_name(status);
        result_.solved = status == Ipopt::SUCCESS;
        result_.acceptable = status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        result_.accepted = result_.solved || result_.acceptable;
        result_.objective = obj_value;
        result_.variables = vector_from_raw(x, n);
        result_.constraints = vector_from_raw(g, m);
        result_.diagnostics_int["solver_status_code"] = static_cast<int>(status);
        result_.diagnostics_int["variables"] = static_cast<int>(n);
        result_.diagnostics_int["constraints"] = static_cast<int>(m);
        result_.diagnostics_bool["exact_gradient_required"] = true;
        result_.diagnostics_bool["exact_jacobian_required"] = true;
        result_.diagnostics_bool["exact_hessian_required"] = false;
    }

    const IpoptSolveResult& result() const {
        return result_;
    }

private:
    const NlpProblem& problem_;
    IpoptSolveOptions options_;
    NlpBounds bounds_;
    std::vector<double> initial_;
    NlpJacobianStructure jacobian_structure_;
    NlpScaling scaling_;
    IpoptSolveResult result_;
};
#endif

}  // namespace

IpoptAdapterInfo native_ipopt_adapter_info() {
    IpoptAdapterInfo info;
#ifdef EPCSAFT_HAS_IPOPT
    info.compiled = true;
    info.adapter_available = true;
    info.status = "enabled_available";
#else
    info.compiled = false;
    info.adapter_available = false;
    info.status = "disabled";
#endif
    return info;
}

IpoptSolveResult solve_ipopt_nlp(
    const NlpProblem& problem,
    const IpoptSolveOptions& options
) {
    validate_nlp_problem_shape(problem);
    if (!options.limited_memory_hessian) {
        throw ValueError("Native Ipopt adapter currently supports only Ipopt limited-memory Hessian mode.");
    }
#ifndef EPCSAFT_HAS_IPOPT
    (void)problem;
    (void)options;
    throw SolutionError("Native Ipopt adapter requires a build configured with EPCSAFT_ENABLE_IPOPT=ON.");
#else
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetIntegerValue("print_level", options.print_level);
    app->Options()->SetIntegerValue("max_iter", options.max_iterations);
    app->Options()->SetNumericValue("tol", options.tolerance);
    app->Options()->SetNumericValue("acceptable_tol", options.acceptable_tolerance);
    app->Options()->SetStringValue("jacobian_approximation", "exact");
    app->Options()->SetStringValue("gradient_approximation", "exact");
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("nlp_scaling_method", "user-scaling");

    const Ipopt::ApplicationReturnStatus init_status = app->Initialize();
    if (init_status != Ipopt::Solve_Succeeded) {
        std::ostringstream msg;
        msg << "Ipopt initialization failed: " << application_status_name(init_status);
        throw SolutionError(msg.str());
    }

    auto* adapter = new IpoptTnlpAdapter(problem, options);
    Ipopt::SmartPtr<Ipopt::TNLP> tnlp = adapter;
    const Ipopt::ApplicationReturnStatus solve_status = app->OptimizeTNLP(tnlp);
    IpoptSolveResult result = adapter->result();
    result.application_status = application_status_name(solve_status);
    result.diagnostics_int["application_status_code"] = static_cast<int>(solve_status);
    return result;
#endif
}

IpoptSolveResult solve_ipopt_quadratic_smoke() {
    QuadraticSmokeProblem problem;
    IpoptSolveOptions options;
    options.max_iterations = 50;
    options.tolerance = 1.0e-10;
    options.acceptable_tolerance = 1.0e-8;
    return solve_ipopt_nlp(problem, options);
}

}  // namespace epcsaft::native::equilibrium_nlp
