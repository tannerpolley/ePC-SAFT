#pragma once

#include "nlp_problem.h"

#include <map>
#include <string>
#include <vector>

namespace epcsaft::native::equilibrium_nlp {

struct IpoptAdapterInfo {
    bool compiled = false;
    bool adapter_available = false;
    std::string backend = "ipopt";
    std::string status = "disabled";
    std::string adapter_kind = "native_tnlp_adapter";
    std::string hessian_strategy = "limited_memory";
    bool exact_gradient_required = true;
    bool exact_jacobian_required = true;
    bool exact_hessian_required = false;
};

struct IpoptSolveOptions {
    int max_iterations = 100;
    int print_level = 0;
    double tolerance = 1.0e-8;
    double acceptable_tolerance = 1.0e-6;
    bool limited_memory_hessian = true;
};

struct IpoptSolveResult {
    bool solver_ran = false;
    bool solved = false;
    bool acceptable = false;
    bool accepted = false;
    std::string backend = "ipopt";
    std::string adapter_kind = "native_tnlp_adapter";
    std::string solver_status;
    std::string application_status;
    std::string hessian_strategy = "limited_memory";
    double objective = 0.0;
    std::vector<double> variables;
    std::vector<double> constraints;
    std::map<std::string, int> diagnostics_int;
    std::map<std::string, double> diagnostics_double;
    std::map<std::string, bool> diagnostics_bool;
    std::map<std::string, std::string> diagnostics_string;
};

IpoptAdapterInfo native_ipopt_adapter_info();

IpoptSolveResult solve_ipopt_nlp(
    const NlpProblem& problem,
    const IpoptSolveOptions& options
);

IpoptSolveResult solve_ipopt_quadratic_smoke();

}  // namespace epcsaft::native::equilibrium_nlp
