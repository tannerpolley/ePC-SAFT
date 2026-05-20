#include "ipopt_adapter.h"

#include "epcsaft_electrolyte.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <exception>
#include <limits>
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

    bool has_exact_hessian() const override {
        return true;
    }

    int hessian_nonzero_count() const override {
        return 3;
    }

    NlpHessianStructure hessian_structure() const override {
        return {{0, 1, 1}, {0, 0, 1}};
    }

    std::vector<double> hessian_values(
        const std::vector<double>& variables,
        double objective_factor,
        const std::vector<double>& constraint_multipliers
    ) const override {
        (void)variables;
        (void)constraint_multipliers;
        return {2.0 * objective_factor, 0.0, 2.0 * objective_factor};
    }

    std::string hessian_backend() const override {
        return "analytic";
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

std::string normalize_hessian_mode(const IpoptSolveOptions& options) {
    std::string mode = options.hessian_mode;
    std::replace(mode.begin(), mode.end(), '_', '-');
    if (mode == "limited-memory" || mode == "exact" || mode == "auto") {
        return mode;
    }
    throw ValueError("Native Ipopt adapter hessian_mode must be 'auto', 'exact', or 'limited-memory'.");
}

double normalize_positive_tolerance(
    double value,
    double default_value,
    const std::string& label
) {
    if (!std::isfinite(default_value) || default_value <= 0.0) {
        throw ValueError("Native Ipopt adapter requires a positive default tolerance for " + label + ".");
    }
    if (value <= 0.0) {
        return default_value;
    }
    if (!std::isfinite(value)) {
        throw ValueError("Native Ipopt adapter " + label + " must be finite.");
    }
    return value;
}

std::string normalize_linear_solver(std::string value) {
    if (value.empty()) {
        return "auto";
    }
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    if (value == "auto") {
        return value;
    }
    if (value.find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_-") != std::string::npos) {
        throw ValueError(
            "Native Ipopt adapter linear_solver must be 'auto' or a non-empty Ipopt solver token."
        );
    }
    return value;
}

IpoptSolveOptions normalize_ipopt_solve_options(const IpoptSolveOptions& options) {
    IpoptSolveOptions normalized = options;
    if (!std::isfinite(normalized.tolerance) || normalized.tolerance <= 0.0) {
        throw ValueError("Native Ipopt adapter tolerance must be positive and finite.");
    }
    normalized.acceptable_tolerance = normalize_positive_tolerance(
        normalized.acceptable_tolerance,
        std::max(100.0 * normalized.tolerance, 1.0e-10),
        "acceptable_tolerance"
    );
    normalized.constraint_violation_tolerance = normalize_positive_tolerance(
        normalized.constraint_violation_tolerance,
        normalized.tolerance,
        "constraint_violation_tolerance"
    );
    normalized.dual_infeasibility_tolerance = normalize_positive_tolerance(
        normalized.dual_infeasibility_tolerance,
        normalized.tolerance,
        "dual_infeasibility_tolerance"
    );
    normalized.complementarity_tolerance = normalize_positive_tolerance(
        normalized.complementarity_tolerance,
        normalized.tolerance,
        "complementarity_tolerance"
    );
    if (!std::isfinite(normalized.bound_push) || normalized.bound_push < 0.0) {
        throw ValueError("Native Ipopt adapter bound_push must be non-negative and finite.");
    }
    if (!std::isfinite(normalized.bound_frac) || normalized.bound_frac < 0.0) {
        throw ValueError("Native Ipopt adapter bound_frac must be non-negative and finite.");
    }
    normalized.linear_solver = normalize_linear_solver(normalized.linear_solver);
    return normalized;
}

#ifdef EPCSAFT_HAS_IPOPT
bool all_finite(const std::vector<double>& values) {
    return std::all_of(values.begin(), values.end(), [](double value) {
        return std::isfinite(value);
    });
}

std::vector<double> vector_from_raw(const Ipopt::Number* values, Ipopt::Index count) {
    return std::vector<double>(values, values + count);
}

void require_vector_size(
    const std::vector<double>& values,
    std::size_t expected,
    const std::string& label
) {
    if (values.empty() || values.size() == expected) {
        return;
    }
    std::ostringstream msg;
    msg << label << " has size " << values.size() << " but expected " << expected << ".";
    throw ValueError(msg.str());
}

double vector_min_or_zero(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    return *std::min_element(values.begin(), values.end());
}

double vector_max_or_zero(const std::vector<double>& values) {
    if (values.empty()) {
        return 0.0;
    }
    return *std::max_element(values.begin(), values.end());
}

std::string solver_status_name(Ipopt::SolverReturn status) {
    switch (status) {
        case Ipopt::SUCCESS:
            return "success";
        case Ipopt::MAXITER_EXCEEDED:
            return "max_iterations_exceeded";
        case Ipopt::CPUTIME_EXCEEDED:
            return "cpu_time_exceeded";
        case Ipopt::WALLTIME_EXCEEDED:
            return "wall_time_exceeded";
        case Ipopt::STOP_AT_TINY_STEP:
            return "tiny_step_detected";
        case Ipopt::STOP_AT_ACCEPTABLE_POINT:
            return "acceptable_point";
        case Ipopt::LOCAL_INFEASIBILITY:
            return "local_infeasibility";
        case Ipopt::USER_REQUESTED_STOP:
            return "user_requested_stop";
        case Ipopt::FEASIBLE_POINT_FOUND:
            return "feasible_point_found";
        case Ipopt::DIVERGING_ITERATES:
            return "diverging_iterates";
        case Ipopt::RESTORATION_FAILURE:
            return "restoration_failure";
        case Ipopt::ERROR_IN_STEP_COMPUTATION:
            return "error_in_step_computation";
        case Ipopt::INVALID_NUMBER_DETECTED:
            return "invalid_number_detected";
        case Ipopt::TOO_FEW_DEGREES_OF_FREEDOM:
            return "too_few_degrees_of_freedom";
        case Ipopt::INVALID_OPTION:
            return "invalid_option";
        case Ipopt::OUT_OF_MEMORY:
            return "out_of_memory";
        case Ipopt::INTERNAL_ERROR:
            return "internal_error";
        case Ipopt::UNASSIGNED:
            return "unassigned";
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
        case Ipopt::Search_Direction_Becomes_Too_Small:
            return "search_direction_too_small";
        case Ipopt::Diverging_Iterates:
            return "diverging_iterates";
        case Ipopt::User_Requested_Stop:
            return "user_requested_stop";
        case Ipopt::Feasible_Point_Found:
            return "feasible_point_found";
        case Ipopt::Maximum_CpuTime_Exceeded:
            return "maximum_cpu_time_exceeded";
        case Ipopt::Maximum_WallTime_Exceeded:
            return "maximum_wall_time_exceeded";
        case Ipopt::Not_Enough_Degrees_Of_Freedom:
            return "not_enough_degrees_of_freedom";
        case Ipopt::Invalid_Problem_Definition:
            return "invalid_problem_definition";
        case Ipopt::Invalid_Option:
            return "invalid_option";
        case Ipopt::Invalid_Number_Detected:
            return "invalid_number_detected";
        case Ipopt::Unrecoverable_Exception:
            return "unrecoverable_exception";
        case Ipopt::NonIpopt_Exception_Thrown:
            return "non_ipopt_exception_thrown";
        case Ipopt::Insufficient_Memory:
            return "insufficient_memory";
        case Ipopt::Internal_Error:
            return "internal_error";
        default:
            return "ipopt_application_status_" + std::to_string(static_cast<int>(status));
    }
}

class IpoptTnlpAdapter final : public Ipopt::TNLP {
public:
    explicit IpoptTnlpAdapter(const NlpProblem& problem, const IpoptSolveOptions& options)
        : problem_(problem),
          options_(options),
          selected_hessian_mode_(select_hessian_mode(problem_, options_)),
          bounds_(problem_.bounds()),
          initial_(options_.initial_variables.empty() ? problem_.initial_point() : options_.initial_variables),
          jacobian_structure_(problem_.jacobian_structure()),
          hessian_structure_(selected_hessian_mode_ == "exact" ? problem_.hessian_structure() : NlpHessianStructure{}),
          scaling_(problem_.scaling()) {
        const auto variable_count = static_cast<std::size_t>(problem_.variable_count());
        const auto constraint_count = static_cast<std::size_t>(problem_.constraint_count());
        require_vector_size(options_.initial_variables, variable_count, "continuation_state.variables");
        require_vector_size(
            options_.initial_bound_lower_multipliers,
            variable_count,
            "continuation_state.bound_lower_multipliers"
        );
        require_vector_size(
            options_.initial_bound_upper_multipliers,
            variable_count,
            "continuation_state.bound_upper_multipliers"
        );
        require_vector_size(
            options_.initial_constraint_multipliers,
            constraint_count,
            "continuation_state.constraint_multipliers"
        );
        result_.variables = initial_;
        result_.diagnostics_string["problem_name"] = problem_.name();
        result_.diagnostics_string["gradient_approximation"] = "exact";
        result_.diagnostics_string["jacobian_approximation"] = "exact";
        result_.diagnostics_string["hessian_approximation"] = selected_hessian_mode_;
        result_.diagnostics_string["hessian_backend"] =
            selected_hessian_mode_ == "exact" ? problem_.hessian_backend() : "limited-memory";
        result_.diagnostics_string["scaling_method"] = "user-scaling";
        result_.diagnostics_string["linear_solver_requested"] = options_.linear_solver;
        result_.diagnostics_string["linear_solver_selected"] =
            options_.linear_solver == "auto" ? "default" : options_.linear_solver;
        result_.diagnostics_int["iteration_history_limit"] = std::max(0, options_.iteration_history_limit);
        result_.diagnostics_int["variable_scaling_count"] = static_cast<int>(scaling_.variables.size());
        result_.diagnostics_int["constraint_scaling_count"] = static_cast<int>(scaling_.constraints.size());
        result_.diagnostics_double["objective_scaling"] = scaling_.objective;
        result_.diagnostics_double["acceptable_tolerance"] = options_.acceptable_tolerance;
        result_.diagnostics_double["constraint_violation_tolerance"] = options_.constraint_violation_tolerance;
        result_.diagnostics_double["dual_infeasibility_tolerance"] = options_.dual_infeasibility_tolerance;
        result_.diagnostics_double["complementarity_tolerance"] = options_.complementarity_tolerance;
        result_.diagnostics_double["bound_push"] = options_.bound_push;
        result_.diagnostics_double["bound_frac"] = options_.bound_frac;
        result_.diagnostics_double["variable_scaling_min"] = vector_min_or_zero(scaling_.variables);
        result_.diagnostics_double["variable_scaling_max"] = vector_max_or_zero(scaling_.variables);
        result_.diagnostics_double["constraint_scaling_min"] = vector_min_or_zero(scaling_.constraints);
        result_.diagnostics_double["constraint_scaling_max"] = vector_max_or_zero(scaling_.constraints);
        result_.diagnostics_bool["exact_hessian_available"] = problem_.has_exact_hessian();
        result_.diagnostics_bool["warm_start_requested"] = has_warm_start();
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
        nnz_h_lag = selected_hessian_mode_ == "exact" ? problem_.hessian_nonzero_count() : 0;
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
            record_callback_failure("get_bounds_info", "dimension mismatch");
            return false;
        }
        std::copy(bounds_.variable_lower.begin(), bounds_.variable_lower.end(), x_l);
        std::copy(bounds_.variable_upper.begin(), bounds_.variable_upper.end(), x_u);
        std::copy(bounds_.constraint_lower.begin(), bounds_.constraint_lower.end(), g_l);
        std::copy(bounds_.constraint_upper.begin(), bounds_.constraint_upper.end(), g_u);
        const bool valid = all_finite(bounds_.variable_lower) && all_finite(bounds_.variable_upper)
            && all_finite(bounds_.constraint_lower) && all_finite(bounds_.constraint_upper);
        if (!valid) {
            record_callback_failure("get_bounds_info", "nonfinite bounds");
        }
        return valid;
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
        if (!init_x || n != problem_.variable_count() || m != problem_.constraint_count()) {
            record_callback_failure("get_starting_point", "dimension mismatch");
            return false;
        }
        std::copy(initial_.begin(), initial_.end(), x);
        if (init_z) {
            if (
                options_.initial_bound_lower_multipliers.size() != initial_.size()
                || options_.initial_bound_upper_multipliers.size() != initial_.size()
            ) {
                record_callback_failure("get_starting_point", "bound multiplier size mismatch");
                return false;
            }
            std::copy(options_.initial_bound_lower_multipliers.begin(), options_.initial_bound_lower_multipliers.end(), z_L);
            std::copy(options_.initial_bound_upper_multipliers.begin(), options_.initial_bound_upper_multipliers.end(), z_U);
        }
        if (init_lambda) {
            if (
                options_.initial_constraint_multipliers.size()
                != static_cast<std::size_t>(problem_.constraint_count())
            ) {
                record_callback_failure("get_starting_point", "constraint multiplier size mismatch");
                return false;
            }
            std::copy(
                options_.initial_constraint_multipliers.begin(),
                options_.initial_constraint_multipliers.end(),
                lambda
            );
        }
        result_.diagnostics_bool["warm_start_used"] = init_z || init_lambda;
        const bool valid = all_finite(initial_);
        if (!valid) {
            record_callback_failure("get_starting_point", "nonfinite initial point");
        }
        return valid;
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
                record_callback_failure("get_scaling_parameters", "variable dimension mismatch");
                return false;
            }
            std::copy(scaling_.variables.begin(), scaling_.variables.end(), x_scaling);
        }
        if (use_g_scaling) {
            if (m != problem_.constraint_count()) {
                record_callback_failure("get_scaling_parameters", "constraint dimension mismatch");
                return false;
            }
            std::copy(scaling_.constraints.begin(), scaling_.constraints.end(), g_scaling);
        }
        const bool valid = std::isfinite(obj_scaling) && obj_scaling > 0.0 && all_finite(scaling_.variables)
            && all_finite(scaling_.constraints);
        if (!valid) {
            record_callback_failure("get_scaling_parameters", "nonpositive or nonfinite scaling");
        }
        return valid;
    }

    bool eval_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number& obj_value
    ) override {
        (void)new_x;
        if (n != problem_.variable_count()) {
            record_callback_failure("eval_f", "dimension mismatch");
            return false;
        }
        try {
            obj_value = problem_.objective(vector_from_raw(x, n));
            const bool valid = std::isfinite(obj_value);
            if (!valid) {
                record_callback_failure("eval_f", "nonfinite objective");
            }
            return valid;
        } catch (const std::exception& exc) {
            record_callback_exception("eval_f", exc.what());
            return false;
        } catch (...) {
            record_callback_exception("eval_f", "unknown exception");
            return false;
        }
    }

    bool eval_grad_f(
        Ipopt::Index n,
        const Ipopt::Number* x,
        bool new_x,
        Ipopt::Number* grad_f
    ) override {
        (void)new_x;
        if (n != problem_.variable_count()) {
            record_callback_failure("eval_grad_f", "dimension mismatch");
            return false;
        }
        try {
            const std::vector<double> gradient = problem_.objective_gradient(vector_from_raw(x, n));
            if (gradient.size() != static_cast<std::size_t>(n) || !all_finite(gradient)) {
                record_callback_failure("eval_grad_f", "invalid gradient shape or nonfinite value");
                return false;
            }
            std::copy(gradient.begin(), gradient.end(), grad_f);
            return true;
        } catch (const std::exception& exc) {
            record_callback_exception("eval_grad_f", exc.what());
            return false;
        } catch (...) {
            record_callback_exception("eval_grad_f", "unknown exception");
            return false;
        }
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
            record_callback_failure("eval_g", "dimension mismatch");
            return false;
        }
        try {
            const std::vector<double> values = problem_.constraints(vector_from_raw(x, n));
            if (values.size() != static_cast<std::size_t>(m) || !all_finite(values)) {
                record_callback_failure("eval_g", "invalid constraint shape or nonfinite value");
                return false;
            }
            std::copy(values.begin(), values.end(), g);
            return true;
        } catch (const std::exception& exc) {
            record_callback_exception("eval_g", exc.what());
            return false;
        } catch (...) {
            record_callback_exception("eval_g", "unknown exception");
            return false;
        }
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
            record_callback_failure("eval_jac_g", "dimension mismatch");
            return false;
        }
        if (values == nullptr) {
            if (jacobian_structure_.rows.size() != static_cast<std::size_t>(nele_jac)
                || jacobian_structure_.cols.size() != static_cast<std::size_t>(nele_jac)) {
                record_callback_failure("eval_jac_g", "invalid Jacobian structure shape");
                return false;
            }
            for (std::size_t index = 0; index < jacobian_structure_.rows.size(); ++index) {
                iRow[index] = jacobian_structure_.rows[index];
                jCol[index] = jacobian_structure_.cols[index];
            }
            return true;
        }
        try {
            const std::vector<double> jacobian = problem_.jacobian_values(vector_from_raw(x, n));
            if (jacobian.size() != static_cast<std::size_t>(nele_jac) || !all_finite(jacobian)) {
                record_callback_failure("eval_jac_g", "invalid Jacobian value shape or nonfinite value");
                return false;
            }
            std::copy(jacobian.begin(), jacobian.end(), values);
            return true;
        } catch (const std::exception& exc) {
            record_callback_exception("eval_jac_g", exc.what());
            return false;
        } catch (...) {
            record_callback_exception("eval_jac_g", "unknown exception");
            return false;
        }
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
        (void)new_x;
        (void)new_lambda;
        if (selected_hessian_mode_ != "exact" || n != problem_.variable_count() || m != problem_.constraint_count()
            || nele_hess != problem_.hessian_nonzero_count()) {
            record_callback_failure("eval_h", "dimension mismatch");
            return false;
        }
        if (values == nullptr) {
            if (hessian_structure_.rows.size() != static_cast<std::size_t>(nele_hess)
                || hessian_structure_.cols.size() != static_cast<std::size_t>(nele_hess)) {
                record_callback_failure("eval_h", "invalid Hessian structure shape");
                return false;
            }
            for (std::size_t index = 0; index < hessian_structure_.rows.size(); ++index) {
                iRow[index] = hessian_structure_.rows[index];
                jCol[index] = hessian_structure_.cols[index];
            }
            return true;
        }
        try {
            const std::vector<double> hessian = problem_.hessian_values(
                vector_from_raw(x, n),
                obj_factor,
                vector_from_raw(lambda, m)
            );
            if (hessian.size() != static_cast<std::size_t>(nele_hess) || !all_finite(hessian)) {
                record_callback_failure("eval_h", "invalid Hessian value shape or nonfinite value");
                return false;
            }
            std::copy(hessian.begin(), hessian.end(), values);
            result_.diagnostics_int["eval_h_calls"] += 1;
            return true;
        } catch (const std::exception& exc) {
            record_callback_exception("eval_h", exc.what());
            return false;
        } catch (...) {
            record_callback_exception("eval_h", "unknown exception");
            return false;
        }
    }

    bool intermediate_callback(
        Ipopt::AlgorithmMode mode,
        Ipopt::Index iter,
        Ipopt::Number obj_value,
        Ipopt::Number inf_pr,
        Ipopt::Number inf_du,
        Ipopt::Number mu,
        Ipopt::Number d_norm,
        Ipopt::Number regularization_size,
        Ipopt::Number alpha_du,
        Ipopt::Number alpha_pr,
        Ipopt::Index ls_trials,
        const Ipopt::IpoptData* ip_data,
        Ipopt::IpoptCalculatedQuantities* ip_cq
    ) override {
        (void)d_norm;
        (void)ip_data;
        (void)ip_cq;
        result_.diagnostics_int["iteration_count"] = static_cast<int>(iter);
        const int limit = std::max(0, options_.iteration_history_limit);
        if (limit > 0) {
            IpoptIterationRecord record;
            record.iteration = static_cast<int>(iter);
            record.objective = obj_value;
            record.primal_infeasibility = inf_pr;
            record.dual_infeasibility = inf_du;
            record.barrier_parameter = mu;
            record.step_size_dual = alpha_du;
            record.step_size_primal = alpha_pr;
            record.regularization_size = regularization_size;
            record.step_trial_count = static_cast<int>(ls_trials);
            record.restoration_phase = mode == Ipopt::RestorationPhaseMode;
            result_.iteration_history.push_back(record);
            while (result_.iteration_history.size() > static_cast<std::size_t>(limit)) {
                result_.iteration_history.erase(result_.iteration_history.begin());
            }
        }
        return true;
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
        (void)ip_data;
        (void)ip_cq;
        result_.solver_ran = true;
        result_.solver_status = solver_status_name(status);
        result_.solved = status == Ipopt::SUCCESS;
        result_.acceptable = status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
        result_.feasible_point = status == Ipopt::FEASIBLE_POINT_FOUND;
        result_.accepted = result_.solved || result_.acceptable;
        result_.objective = obj_value;
        result_.variables = vector_from_raw(x, n);
        result_.constraints = vector_from_raw(g, m);
        result_.bound_lower_multipliers = vector_from_raw(z_L, n);
        result_.bound_upper_multipliers = vector_from_raw(z_U, n);
        result_.constraint_multipliers = vector_from_raw(lambda, m);
        result_.diagnostics_int["solver_status_code"] = static_cast<int>(status);
        result_.diagnostics_int["variables"] = static_cast<int>(n);
        result_.diagnostics_int["constraints"] = static_cast<int>(m);
        result_.diagnostics_int["iteration_history_size"] = static_cast<int>(result_.iteration_history.size());
        result_.diagnostics_bool["exact_gradient_required"] = true;
        result_.diagnostics_bool["exact_jacobian_required"] = true;
        result_.diagnostics_bool["feasible_point_found"] = result_.feasible_point;
    }

    const IpoptSolveResult& result() const {
        return result_;
    }

private:
    void record_callback_failure(const std::string& callback, const std::string& message) {
        result_.diagnostics_string["last_callback_failure"] = callback + ": " + message;
    }

    void record_callback_exception(const std::string& callback, const std::string& message) {
        result_.diagnostics_string["last_callback_exception"] = callback + ": " + message;
    }

    bool has_warm_start() const {
        return !options_.initial_variables.empty()
            || !options_.initial_bound_lower_multipliers.empty()
            || !options_.initial_bound_upper_multipliers.empty()
            || !options_.initial_constraint_multipliers.empty();
    }

    static std::string select_hessian_mode(const NlpProblem& problem, const IpoptSolveOptions& options) {
        const std::string mode = normalize_hessian_mode(options);
        if (mode == "auto") {
            if (problem.has_exact_hessian()) {
                return "exact";
            }
            throw ValueError("Native Ipopt auto Hessian mode requires an NLP Hessian provider; request a supported non-exact mode explicitly to opt out.");
        }
        if (mode == "exact" && !problem.has_exact_hessian()) {
            throw ValueError("Native Ipopt exact Hessian mode requires an NLP Hessian provider.");
        }
        return mode;
    }

    const NlpProblem& problem_;
    IpoptSolveOptions options_;
    std::string selected_hessian_mode_;
    NlpBounds bounds_;
    std::vector<double> initial_;
    NlpJacobianStructure jacobian_structure_;
    NlpHessianStructure hessian_structure_;
    NlpScaling scaling_;
    IpoptSolveResult result_;
};
#endif

}  // namespace

std::string solve_diagnostic_string(
    const IpoptSolveResult& solve,
    const std::string& key,
    const std::string& default_value
) {
    const auto item = solve.diagnostics_string.find(key);
    return item == solve.diagnostics_string.end() ? default_value : item->second;
}

int solve_diagnostic_int(
    const IpoptSolveResult& solve,
    const std::string& key,
    int default_value
) {
    const auto item = solve.diagnostics_int.find(key);
    return item == solve.diagnostics_int.end() ? default_value : item->second;
}

double solve_diagnostic_double(
    const IpoptSolveResult& solve,
    const std::string& key,
    double default_value
) {
    const auto item = solve.diagnostics_double.find(key);
    return item == solve.diagnostics_double.end() ? default_value : item->second;
}

bool solve_diagnostic_bool(
    const IpoptSolveResult& solve,
    const std::string& key,
    bool default_value
) {
    const auto item = solve.diagnostics_bool.find(key);
    return item == solve.diagnostics_bool.end() ? default_value : item->second;
}

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
    const IpoptSolveOptions normalized_options = normalize_ipopt_solve_options(options);
    const std::string normalized_hessian_mode = normalize_hessian_mode(normalized_options);
    std::string selected_hessian_mode = normalized_hessian_mode;
    if (normalized_hessian_mode == "auto") {
        if (!problem.has_exact_hessian()) {
            throw ValueError("Native Ipopt auto Hessian mode requires an NLP Hessian provider; request a supported non-exact mode explicitly to opt out.");
        }
        selected_hessian_mode = "exact";
    }
    if (selected_hessian_mode == "exact" && !problem.has_exact_hessian()) {
        throw ValueError("Native Ipopt exact Hessian mode requires an NLP Hessian provider.");
    }
    if (!std::isfinite(normalized_options.max_wall_time_seconds) || normalized_options.max_wall_time_seconds < 0.0) {
        throw ValueError("Native Ipopt adapter requires a non-negative wall-clock timeout.");
    }
#ifndef EPCSAFT_HAS_IPOPT
    (void)problem;
    (void)options;
    throw SolutionError("Native Ipopt adapter requires a build configured with EPCSAFT_ENABLE_IPOPT=ON.");
#else
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    app->Options()->SetIntegerValue("print_level", normalized_options.print_level);
    app->Options()->SetIntegerValue("max_iter", normalized_options.max_iterations);
    app->Options()->SetNumericValue("tol", normalized_options.tolerance);
    app->Options()->SetNumericValue("acceptable_tol", normalized_options.acceptable_tolerance);
    app->Options()->SetNumericValue("constr_viol_tol", normalized_options.constraint_violation_tolerance);
    app->Options()->SetNumericValue("dual_inf_tol", normalized_options.dual_infeasibility_tolerance);
    app->Options()->SetNumericValue("compl_inf_tol", normalized_options.complementarity_tolerance);
    app->Options()->SetStringValue("jacobian_approximation", "exact");
    app->Options()->SetStringValue("gradient_approximation", "exact");
    if (normalized_options.bound_push > 0.0) {
        app->Options()->SetNumericValue("bound_push", normalized_options.bound_push);
    }
    if (normalized_options.bound_frac > 0.0) {
        app->Options()->SetNumericValue("bound_frac", normalized_options.bound_frac);
    }
    if (selected_hessian_mode == "limited-memory") {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }
    if (normalized_options.linear_solver != "auto") {
        app->Options()->SetStringValue("linear_solver", normalized_options.linear_solver);
    }
    app->Options()->SetStringValue("nlp_scaling_method", "user-scaling");
    if (
        !normalized_options.initial_bound_lower_multipliers.empty()
        || !normalized_options.initial_bound_upper_multipliers.empty()
        || !normalized_options.initial_constraint_multipliers.empty()
    ) {
        app->Options()->SetStringValue("warm_start_init_point", "yes");
    }
    if (normalized_options.max_wall_time_seconds > 0.0) {
        app->Options()->SetNumericValue("max_wall_time", normalized_options.max_wall_time_seconds);
    }

    const Ipopt::ApplicationReturnStatus init_status = app->Initialize();
    if (init_status != Ipopt::Solve_Succeeded) {
        std::ostringstream msg;
        if (normalized_options.linear_solver != "auto") {
            msg << "Ipopt initialization failed for linear_solver='" << normalized_options.linear_solver
                << "': " << application_status_name(init_status);
        } else {
            msg << "Ipopt initialization failed: " << application_status_name(init_status);
        }
        throw SolutionError(msg.str());
    }

    auto* adapter = new IpoptTnlpAdapter(problem, normalized_options);
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

IpoptSolveResult solve_ipopt_quadratic_smoke(const IpoptSolveOptions& options) {
    QuadraticSmokeProblem problem;
    return solve_ipopt_nlp(problem, options);
}

}  // namespace epcsaft::native::equilibrium_nlp
