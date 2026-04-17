#include "epcsaft_core_internal.h"

#include <algorithm>
#include <cmath>
#include <random>
#include <sstream>
#include <string>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>

using Ipopt::AlgorithmMode;
using Ipopt::ApplicationReturnStatus;
using Ipopt::Index;
using Ipopt::IpoptApplication;
using Ipopt::Number;
using Ipopt::SmartPtr;
using Ipopt::SolverReturn;
using Ipopt::TNLP;
using namespace thermo_detail;

namespace {

constexpr double kRegressionGradientFloor = 1.0e-300;

template <typename Scalar>
struct PureNeutralStateScalar {
    Scalar den = scalar_constant<Scalar>(0.0);
    Scalar d = scalar_constant<Scalar>(0.0);
    Scalar zeta0 = scalar_constant<Scalar>(0.0);
    Scalar zeta1 = scalar_constant<Scalar>(0.0);
    Scalar zeta2 = scalar_constant<Scalar>(0.0);
    Scalar zeta3 = scalar_constant<Scalar>(0.0);
    Scalar eta = scalar_constant<Scalar>(0.0);
    Scalar pair_diameter = scalar_constant<Scalar>(0.0);
    Scalar ghs = scalar_constant<Scalar>(0.0);
    Scalar ares_hs = scalar_constant<Scalar>(0.0);
    Scalar ares_hc = scalar_constant<Scalar>(0.0);
    Scalar I1 = scalar_constant<Scalar>(0.0);
    Scalar I2 = scalar_constant<Scalar>(0.0);
    Scalar dEtaI1_deta = scalar_constant<Scalar>(0.0);
    Scalar dEtaI2_deta = scalar_constant<Scalar>(0.0);
    Scalar C1 = scalar_constant<Scalar>(0.0);
    Scalar C2 = scalar_constant<Scalar>(0.0);
    Scalar m2es3 = scalar_constant<Scalar>(0.0);
    Scalar m2e2s3 = scalar_constant<Scalar>(0.0);
    Scalar ares_disp = scalar_constant<Scalar>(0.0);
    Scalar ares_total = scalar_constant<Scalar>(0.0);
    Scalar zraw_hc = scalar_constant<Scalar>(0.0);
    Scalar zraw_disp = scalar_constant<Scalar>(0.0);
    Scalar zraw_total = scalar_constant<Scalar>(0.0);
    Scalar Z = scalar_constant<Scalar>(0.0);
    Scalar pressure = scalar_constant<Scalar>(0.0);
    Scalar lnfug = scalar_constant<Scalar>(0.0);
};

struct PureNeutralObjectiveEvaluation {
    double objective = 0.0;
    std::array<double, 3> gradient{0.0, 0.0, 0.0};
    vector<double> density_raw_residuals;
    vector<double> pure_vle_raw_residuals;
};

void validate_pure_neutral_base_args_cpp(const add_args &base_args) {
    if (base_args.m.size() != 1 || base_args.s.size() != 1 || base_args.e.size() != 1) {
        throw ValueError("Native IPOPT pure-neutral regression requires exactly one component.");
    }
    if (!base_args.z.empty() && (base_args.z.size() != 1 || std::abs(base_args.z[0]) > 1.0e-12)) {
        throw ValueError("Native IPOPT pure-neutral regression currently supports only neutral single-component models.");
    }
    if (!base_args.assoc_num.empty() || !base_args.assoc_matrix.empty() || !base_args.k_hb.empty() || !base_args.e_assoc.empty() || !base_args.vol_a.empty()) {
        throw ValueError("Native IPOPT pure-neutral regression currently supports only nonassociating single-component models.");
    }
    if (base_args.mw.size() != 1) {
        throw ValueError("Native IPOPT pure-neutral regression requires a single MW value in the fixed parameter payload.");
    }
}

template <typename Scalar>
PureNeutralStateScalar<Scalar> pure_neutral_state_scalar_cpp(
    double t,
    const Scalar &rho,
    const Scalar &m,
    const Scalar &s,
    const Scalar &e
) {
    PureNeutralStateScalar<Scalar> state;
    state.den = rho * (N_AV / 1.0e30);
    state.d = s * (1.0 - 0.12 * scalar_exp(-3.0 * e / t));

    const Scalar prefactor = PI / 6.0 * state.den * m;
    state.zeta0 = prefactor;
    state.zeta1 = prefactor * state.d;
    state.zeta2 = prefactor * scalar_pow(state.d, 2.0);
    state.zeta3 = prefactor * scalar_pow(state.d, 3.0);
    state.eta = state.zeta3;
    state.pair_diameter = state.d / 2.0;
    state.ghs = 1.0 / (1.0 - state.zeta3)
        + state.pair_diameter * 3.0 * state.zeta2 / scalar_pow(1.0 - state.zeta3, 2.0)
        + scalar_pow(state.pair_diameter, 2.0) * 2.0 * state.zeta2 * state.zeta2 / scalar_pow(1.0 - state.zeta3, 3.0);

    state.ares_hs = 1.0 / state.zeta0 * (
        3.0 * state.zeta1 * state.zeta2 / (1.0 - state.zeta3)
        + scalar_pow(state.zeta2, 3.0) / (state.zeta3 * scalar_pow(1.0 - state.zeta3, 2.0))
        + (scalar_pow(state.zeta2, 3.0) / scalar_pow(state.zeta3, 2.0) - state.zeta0) * scalar_log(1.0 - state.zeta3)
    );
    state.ares_hc = m * state.ares_hs - (m - 1.0) * scalar_log(state.ghs);

    const Scalar c1 = (m - 1.0) / m;
    const Scalar c2 = (m - 2.0) / m;
    for (size_t i = 0; i < kDispersionA0.size(); ++i) {
        Scalar a_i = kDispersionA0[i] + c1 * kDispersionA1[i] + c1 * c2 * kDispersionA2[i];
        Scalar b_i = kDispersionB0[i] + c1 * kDispersionB1[i] + c1 * c2 * kDispersionB2[i];
        state.I1 += a_i * scalar_pow(state.eta, static_cast<int>(i));
        state.I2 += b_i * scalar_pow(state.eta, static_cast<int>(i));
        state.dEtaI1_deta += a_i * static_cast<double>(i + 1) * scalar_pow(state.eta, static_cast<int>(i));
        state.dEtaI2_deta += b_i * static_cast<double>(i + 1) * scalar_pow(state.eta, static_cast<int>(i));
    }

    state.C1 = 1.0 / (
        1.0
        + m * (8.0 * state.eta - 2.0 * state.eta * state.eta) / scalar_pow(1.0 - state.eta, 4.0)
        + (1.0 - m) * (
            20.0 * state.eta
            - 27.0 * state.eta * state.eta
            + 12.0 * scalar_pow(state.eta, 3.0)
            - 2.0 * scalar_pow(state.eta, 4.0)
        ) / scalar_pow((1.0 - state.eta) * (2.0 - state.eta), 2.0)
    );
    state.C2 = -state.C1 * state.C1 * (
        m * (-4.0 * state.eta * state.eta + 20.0 * state.eta + 8.0) / scalar_pow(1.0 - state.eta, 5.0)
        + (1.0 - m) * (
            2.0 * scalar_pow(state.eta, 3.0)
            + 12.0 * state.eta * state.eta
            - 48.0 * state.eta
            + 40.0
        ) / scalar_pow((1.0 - state.eta) * (2.0 - state.eta), 3.0)
    );

    state.m2es3 = m * m * (e / t) * scalar_pow(s, 3.0);
    state.m2e2s3 = m * m * scalar_pow(e / t, 2.0) * scalar_pow(s, 3.0);
    state.ares_disp = -2.0 * PI * state.den * state.I1 * state.m2es3
        - PI * state.den * m * state.C1 * state.I2 * state.m2e2s3;
    state.ares_total = state.ares_hc + state.ares_disp;

    const Scalar dghs_drho = state.zeta3 / scalar_pow(1.0 - state.zeta3, 2.0)
        + state.pair_diameter * (
            3.0 * state.zeta2 / scalar_pow(1.0 - state.zeta3, 2.0)
            + 6.0 * state.zeta2 * state.zeta3 / scalar_pow(1.0 - state.zeta3, 3.0)
        )
        + scalar_pow(state.pair_diameter, 2.0) * (
            4.0 * state.zeta2 * state.zeta2 / scalar_pow(1.0 - state.zeta3, 3.0)
            + 6.0 * state.zeta2 * state.zeta2 * state.zeta3 / scalar_pow(1.0 - state.zeta3, 4.0)
        );
    const Scalar dadrho_hs = state.zeta3 / (1.0 - state.zeta3)
        + 3.0 * state.zeta1 * state.zeta2 / state.zeta0 / scalar_pow(1.0 - state.zeta3, 2.0)
        + (3.0 * scalar_pow(state.zeta2, 3.0) - state.zeta3 * scalar_pow(state.zeta2, 3.0))
            / state.zeta0 / scalar_pow(1.0 - state.zeta3, 3.0);
    state.zraw_hc = m * dadrho_hs - (m - 1.0) * dghs_drho / state.ghs;
    state.zraw_disp = -2.0 * PI * state.den * state.dEtaI1_deta * state.m2es3
        - PI * state.den * m * (state.C1 * state.dEtaI2_deta + state.C2 * state.eta * state.I2) * state.m2e2s3;
    state.zraw_total = state.zraw_hc + state.zraw_disp;
    state.Z = 1.0 + state.zraw_total;
    if (!(scalar_value(state.Z) > 0.0)) {
        throw ValueError("Encountered non-positive compressibility factor during native IPOPT regression evaluation.");
    }
    state.pressure = state.Z * kb * t * state.den * 1.0e30;
    state.lnfug = state.ares_total + state.zraw_total - scalar_log(state.Z);
    return state;
}

add_args pure_neutral_args_with_theta_cpp(const add_args &base_args, const vector<double> &x) {
    if (x.size() != 3) {
        throw ValueError("Native IPOPT pure-neutral regression expects exactly three optimization variables.");
    }
    add_args args = base_args;
    args.m[0] = x[0];
    args.s[0] = x[1];
    args.e[0] = x[2];
    return args;
}

double pure_neutral_pressure_cpp(double t, double rho, const vector<double> &x) {
    auto state = pure_neutral_state_scalar_cpp<double>(t, rho, x[0], x[1], x[2]);
    return state.pressure;
}

double pure_neutral_lnfug_cpp(double t, double rho, const vector<double> &x) {
    auto state = pure_neutral_state_scalar_cpp<double>(t, rho, x[0], x[1], x[2]);
    return state.lnfug;
}

double pure_neutral_pressure_rho_derivative_cpp(double t, double rho, const vector<double> &x) {
    AutoDual rho_dual = make_autodiff_scalar(rho, 1.0);
    auto state = pure_neutral_state_scalar_cpp<AutoDual>(t, rho_dual, x[0], x[1], x[2]);
    return scalar_derivative(state.pressure);
}

double pure_neutral_lnfug_rho_derivative_cpp(double t, double rho, const vector<double> &x) {
    AutoDual rho_dual = make_autodiff_scalar(rho, 1.0);
    auto state = pure_neutral_state_scalar_cpp<AutoDual>(t, rho_dual, x[0], x[1], x[2]);
    return scalar_derivative(state.lnfug);
}

double pure_neutral_pressure_parameter_derivative_cpp(double t, double rho, const vector<double> &x, int parameter_index) {
    AutoDual m = make_autodiff_scalar(x[0], parameter_index == 0 ? 1.0 : 0.0);
    AutoDual s = make_autodiff_scalar(x[1], parameter_index == 1 ? 1.0 : 0.0);
    AutoDual e = make_autodiff_scalar(x[2], parameter_index == 2 ? 1.0 : 0.0);
    auto state = pure_neutral_state_scalar_cpp<AutoDual>(t, rho, m, s, e);
    return scalar_derivative(state.pressure);
}

double pure_neutral_lnfug_parameter_derivative_cpp(double t, double rho, const vector<double> &x, int parameter_index) {
    AutoDual m = make_autodiff_scalar(x[0], parameter_index == 0 ? 1.0 : 0.0);
    AutoDual s = make_autodiff_scalar(x[1], parameter_index == 1 ? 1.0 : 0.0);
    AutoDual e = make_autodiff_scalar(x[2], parameter_index == 2 ? 1.0 : 0.0);
    auto state = pure_neutral_state_scalar_cpp<AutoDual>(t, rho, m, s, e);
    return scalar_derivative(state.lnfug);
}

double clip_start_value_cpp(double value, double lower, double upper) {
    double clipped = value;
    double margin = 1.0e-8 * std::max(1.0, std::abs(value));
    if (clipped <= lower) {
        clipped = lower + margin;
    }
    if (clipped >= upper) {
        clipped = upper - margin;
    }
    return clipped;
}

vector<vector<double>> candidate_starts_cpp(
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart
) {
    if (x0.size() != lower.size() || x0.size() != upper.size()) {
        throw ValueError("Initial guess and bounds must have matching lengths for native IPOPT regression.");
    }
    vector<vector<double>> starts;
    vector<double> first = x0;
    for (size_t i = 0; i < first.size(); ++i) {
        first[i] = clip_start_value_cpp(first[i], lower[i], upper[i]);
    }
    starts.push_back(first);

    std::mt19937 rng(12345);
    for (int k = 0; k < multistart; ++k) {
        vector<double> point(x0.size(), 0.0);
        for (size_t i = 0; i < x0.size(); ++i) {
            std::uniform_real_distribution<double> uniform(lower[i], upper[i]);
            point[i] = clip_start_value_cpp(uniform(rng), lower[i], upper[i]);
        }
        starts.push_back(point);
    }
    return starts;
}

double rms_metric_cpp(const vector<double> &values) {
    if (values.empty()) {
        return 0.0;
    }
    double accum = 0.0;
    for (double value : values) {
        accum += value * value;
    }
    return std::sqrt(accum / static_cast<double>(values.size()));
}

bool ipopt_status_success_cpp(SolverReturn status) {
    return status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT;
}

std::string ipopt_status_message_cpp(SolverReturn status) {
    switch (status) {
        case Ipopt::SUCCESS:
            return "Ipopt terminated successfully.";
        case Ipopt::STOP_AT_ACCEPTABLE_POINT:
            return "Ipopt terminated at an acceptable point.";
        case Ipopt::MAXITER_EXCEEDED:
            return "Ipopt exceeded the maximum iteration count.";
        case Ipopt::LOCAL_INFEASIBILITY:
            return "Ipopt converged to a locally infeasible point.";
        case Ipopt::INVALID_NUMBER_DETECTED:
            return "Ipopt detected an invalid number during evaluation.";
        case Ipopt::DIVERGING_ITERATES:
            return "Ipopt detected diverging iterates.";
        case Ipopt::RESTORATION_FAILURE:
            return "Ipopt restoration phase failed.";
        case Ipopt::ERROR_IN_STEP_COMPUTATION:
            return "Ipopt encountered an unrecoverable step-computation error.";
        default:
            return "Ipopt terminated without a successful solution.";
    }
    return "Ipopt terminated without a successful solution.";
}

PureNeutralObjectiveEvaluation evaluate_pure_neutral_objective_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x
) {
    PureNeutralObjectiveEvaluation eval;
    const vector<double> one_x = {1.0};
    add_args args = pure_neutral_args_with_theta_cpp(base_args, x);

    eval.density_raw_residuals.reserve(density_records.size());
    eval.pure_vle_raw_residuals.reserve(pure_vle_records.size());

    for (const auto &record : density_records) {
        double rho_calc = den_cpp(record.t, record.p, one_x, record.phase, args);
        double denom = std::max(std::abs(record.rho_exp), kRegressionGradientFloor);
        double raw = (rho_calc - record.rho_exp) / denom;
        eval.density_raw_residuals.push_back(raw);
        double scaled = density_scale * raw;
        eval.objective += 0.5 * scaled * scaled;

        double dpdrho = pure_neutral_pressure_rho_derivative_cpp(record.t, rho_calc, x);
        if (!(std::isfinite(dpdrho) && std::abs(dpdrho) > 0.0)) {
            throw ValueError("Encountered invalid exact dp/drho during native IPOPT density residual evaluation.");
        }
        for (int j = 0; j < 3; ++j) {
            double dpdtheta = pure_neutral_pressure_parameter_derivative_cpp(record.t, rho_calc, x, j);
            double drho_dtheta = -dpdtheta / dpdrho;
            double raw_grad = drho_dtheta / denom;
            eval.gradient[j] += scaled * (density_scale * raw_grad);
        }
    }

    for (const auto &record : pure_vle_records) {
        double rho_liq = den_cpp(record.t, record.p, one_x, 0, args);
        double rho_vap = den_cpp(record.t, record.p, one_x, 1, args);

        double lnphi_liq = pure_neutral_lnfug_cpp(record.t, rho_liq, x);
        double lnphi_vap = pure_neutral_lnfug_cpp(record.t, rho_vap, x);
        double raw = lnphi_liq - lnphi_vap;
        eval.pure_vle_raw_residuals.push_back(raw);
        double scaled = pure_vle_scale * raw;
        eval.objective += 0.5 * scaled * scaled;

        double dlnphi_liq_drho = pure_neutral_lnfug_rho_derivative_cpp(record.t, rho_liq, x);
        double dlnphi_vap_drho = pure_neutral_lnfug_rho_derivative_cpp(record.t, rho_vap, x);
        double dpdrho_liq = pure_neutral_pressure_rho_derivative_cpp(record.t, rho_liq, x);
        double dpdrho_vap = pure_neutral_pressure_rho_derivative_cpp(record.t, rho_vap, x);
        if (!(std::isfinite(dpdrho_liq) && std::abs(dpdrho_liq) > 0.0 && std::isfinite(dpdrho_vap) && std::abs(dpdrho_vap) > 0.0)) {
            throw ValueError("Encountered invalid exact dp/drho during native IPOPT fugacity-balance evaluation.");
        }

        for (int j = 0; j < 3; ++j) {
            double dlnphi_liq_dtheta = pure_neutral_lnfug_parameter_derivative_cpp(record.t, rho_liq, x, j);
            double dlnphi_vap_dtheta = pure_neutral_lnfug_parameter_derivative_cpp(record.t, rho_vap, x, j);
            double dpdtheta_liq = pure_neutral_pressure_parameter_derivative_cpp(record.t, rho_liq, x, j);
            double dpdtheta_vap = pure_neutral_pressure_parameter_derivative_cpp(record.t, rho_vap, x, j);
            double drho_liq_dtheta = -dpdtheta_liq / dpdrho_liq;
            double drho_vap_dtheta = -dpdtheta_vap / dpdrho_vap;
            double total_grad = (
                dlnphi_liq_dtheta + dlnphi_liq_drho * drho_liq_dtheta
                - dlnphi_vap_dtheta - dlnphi_vap_drho * drho_vap_dtheta
            );
            eval.gradient[j] += scaled * (pure_vle_scale * total_grad);
        }
    }

    return eval;
}

class PureNeutralRegressionTNLP : public TNLP {
public:
    PureNeutralRegressionTNLP(
        add_args base_args,
        vector<PureNeutralRegressionDensityRecord> density_records,
        double density_scale,
        vector<PureNeutralRegressionVLERecord> pure_vle_records,
        double pure_vle_scale,
        vector<double> start,
        vector<double> lower,
        vector<double> upper
    )
        : base_args_(std::move(base_args)),
          density_records_(std::move(density_records)),
          density_scale_(density_scale),
          pure_vle_records_(std::move(pure_vle_records)),
          pure_vle_scale_(pure_vle_scale),
          start_(std::move(start)),
          lower_(std::move(lower)),
          upper_(std::move(upper)) {}

    bool get_nlp_info(Index &n, Index &m, Index &nnz_jac_g, Index &nnz_h_lag, IndexStyleEnum &index_style) override {
        n = 3;
        m = 0;
        nnz_jac_g = 0;
        nnz_h_lag = 0;
        index_style = TNLP::C_STYLE;
        return true;
    }

    bool get_bounds_info(Index n, Number *x_l, Number *x_u, Index m, Number *g_l, Number *g_u) override {
        (void)m;
        (void)g_l;
        (void)g_u;
        for (Index i = 0; i < n; ++i) {
            x_l[i] = lower_[static_cast<size_t>(i)];
            x_u[i] = upper_[static_cast<size_t>(i)];
        }
        return true;
    }

    bool get_starting_point(Index n, bool init_x, Number *x, bool init_z, Number *z_L, Number *z_U, Index m, bool init_lambda, Number *lambda) override {
        (void)init_z;
        (void)z_L;
        (void)z_U;
        (void)m;
        (void)init_lambda;
        (void)lambda;
        if (!init_x) {
            return false;
        }
        for (Index i = 0; i < n; ++i) {
            x[i] = start_[static_cast<size_t>(i)];
        }
        return true;
    }

    bool eval_f(Index n, const Number *x, bool new_x, Number &obj_value) override {
        if (!ensure_cache(n, x, new_x)) {
            return false;
        }
        obj_value = cache_.objective;
        return true;
    }

    bool eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) override {
        if (!ensure_cache(n, x, new_x)) {
            return false;
        }
        for (Index i = 0; i < n; ++i) {
            grad_f[i] = cache_.gradient[static_cast<size_t>(i)];
        }
        return true;
    }

    bool eval_g(Index n, const Number *x, bool new_x, Index m, Number *g) override {
        (void)n;
        (void)x;
        (void)new_x;
        (void)m;
        (void)g;
        return true;
    }

    bool eval_jac_g(Index n, const Number *x, bool new_x, Index m, Index nele_jac, Index *iRow, Index *jCol, Number *values) override {
        (void)n;
        (void)x;
        (void)new_x;
        (void)m;
        (void)nele_jac;
        (void)iRow;
        (void)jCol;
        (void)values;
        return true;
    }

    bool get_list_of_nonlinear_variables(Index num_nonlin_vars, Index *pos_nonlin_vars) override {
        if (num_nonlin_vars != 3) {
            return false;
        }
        pos_nonlin_vars[0] = 0;
        pos_nonlin_vars[1] = 1;
        pos_nonlin_vars[2] = 2;
        return true;
    }

    Index get_number_of_nonlinear_variables() override {
        return 3;
    }

    bool intermediate_callback(
        AlgorithmMode mode,
        Index iter,
        Number obj_value,
        Number inf_pr,
        Number inf_du,
        Number mu,
        Number d_norm,
        Number regularization_size,
        Number alpha_du,
        Number alpha_pr,
        Index ls_trials,
        const Ipopt::IpoptData *ip_data,
        Ipopt::IpoptCalculatedQuantities *ip_cq
    ) override {
        (void)mode;
        (void)obj_value;
        (void)inf_pr;
        (void)inf_du;
        (void)mu;
        (void)d_norm;
        (void)regularization_size;
        (void)alpha_du;
        (void)alpha_pr;
        (void)ls_trials;
        (void)ip_data;
        (void)ip_cq;
        iterations_ = static_cast<int>(iter);
        return true;
    }

    void finalize_solution(
        SolverReturn status,
        Index n,
        const Number *x,
        const Number *z_L,
        const Number *z_U,
        Index m,
        const Number *g,
        const Number *lambda,
        Number obj_value,
        const Ipopt::IpoptData *ip_data,
        Ipopt::IpoptCalculatedQuantities *ip_cq
    ) override {
        (void)z_L;
        (void)z_U;
        (void)m;
        (void)g;
        (void)lambda;
        (void)ip_data;
        (void)ip_cq;
        status_ = status;
        obj_value_ = obj_value;
        solution_.assign(x, x + n);
        if (!cache_valid_ || last_x_ != solution_) {
            try {
                cache_ = evaluate_pure_neutral_objective_cpp(
                    base_args_,
                    density_records_,
                    density_scale_,
                    pure_vle_records_,
                    pure_vle_scale_,
                    solution_
                );
                cache_valid_ = true;
                last_x_ = solution_;
            } catch (...) {
                cache_valid_ = false;
            }
        }
    }

    PureNeutralRegressionResult result() const {
        PureNeutralRegressionResult out;
        out.x = solution_.empty() ? start_ : solution_;
        out.success = ipopt_status_success_cpp(status_);
        out.status = static_cast<int>(status_);
        out.message = ipopt_status_message_cpp(status_);
        out.nfev = evaluations_;
        out.iterations = iterations_;
        if (cache_valid_) {
            out.cost = cache_.objective;
            out.residual_norm = std::sqrt(std::max(0.0, 2.0 * cache_.objective));
            out.density_metric = rms_metric_cpp(cache_.density_raw_residuals);
            out.pure_vle_metric = rms_metric_cpp(cache_.pure_vle_raw_residuals);
        } else {
            out.cost = std::isfinite(obj_value_) ? obj_value_ : HUGE_DBL;
            out.residual_norm = HUGE_DBL;
            out.density_metric = HUGE_DBL;
            out.pure_vle_metric = HUGE_DBL;
        }
        return out;
    }

private:
    bool ensure_cache(Index n, const Number *x, bool new_x) {
        vector<double> current(static_cast<size_t>(n), 0.0);
        for (Index i = 0; i < n; ++i) {
            current[static_cast<size_t>(i)] = x[i];
        }
        if (!cache_valid_ || new_x || current != last_x_) {
            try {
                cache_ = evaluate_pure_neutral_objective_cpp(
                    base_args_,
                    density_records_,
                    density_scale_,
                    pure_vle_records_,
                    pure_vle_scale_,
                    current
                );
            } catch (...) {
                return false;
            }
            cache_valid_ = true;
            last_x_ = std::move(current);
            ++evaluations_;
        }
        return true;
    }

    add_args base_args_;
    vector<PureNeutralRegressionDensityRecord> density_records_;
    double density_scale_ = 1.0;
    vector<PureNeutralRegressionVLERecord> pure_vle_records_;
    double pure_vle_scale_ = 1.0;
    vector<double> start_;
    vector<double> lower_;
    vector<double> upper_;
    PureNeutralObjectiveEvaluation cache_;
    vector<double> last_x_;
    vector<double> solution_;
    bool cache_valid_ = false;
    int evaluations_ = 0;
    int iterations_ = 0;
    SolverReturn status_ = Ipopt::INTERNAL_ERROR;
    double obj_value_ = HUGE_DBL;
};

PureNeutralRegressionResult solve_one_start_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &start,
    const vector<double> &lower,
    const vector<double> &upper,
    bool derivative_test
) {
    SmartPtr<PureNeutralRegressionTNLP> nlp = new PureNeutralRegressionTNLP(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        start,
        lower,
        upper
    );
    SmartPtr<IpoptApplication> app = ::IpoptApplicationFactory();
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("sb", "yes");
    app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetIntegerValue("max_iter", 500);
    app->Options()->SetNumericValue("tol", 1.0e-8);
    app->Options()->SetNumericValue("acceptable_tol", 1.0e-6);
    if (derivative_test) {
        app->Options()->SetStringValue("derivative_test", "first-order");
        app->Options()->SetNumericValue("derivative_test_perturbation", 1.0e-7);
    }
    ApplicationReturnStatus init_status = app->Initialize();
    if (init_status != Ipopt::Solve_Succeeded) {
        throw ValueError("Failed to initialize IPOPT application for native pure-neutral regression.");
    }
    app->OptimizeTNLP(GetRawPtr(nlp));
    return nlp->result();
}

}  // namespace

PureNeutralRegressionDebugResult evaluate_pure_neutral_objective_debug_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x
) {
    validate_pure_neutral_base_args_cpp(base_args);
    PureNeutralObjectiveEvaluation eval = evaluate_pure_neutral_objective_cpp(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        x
    );
    PureNeutralRegressionDebugResult out;
    out.objective = eval.objective;
    out.gradient.assign(eval.gradient.begin(), eval.gradient.end());
    out.density_raw_residuals = std::move(eval.density_raw_residuals);
    out.pure_vle_raw_residuals = std::move(eval.pure_vle_raw_residuals);
    return out;
}

PureNeutralRegressionResult fit_pure_neutral_ipopt_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart,
    bool derivative_test
) {
    validate_pure_neutral_base_args_cpp(base_args);
    if (density_records.empty() || pure_vle_records.empty()) {
        throw ValueError("Native IPOPT pure-neutral regression requires both density and pure-VLE record families.");
    }
    if (x0.size() != 3 || lower.size() != 3 || upper.size() != 3) {
        throw ValueError("Native IPOPT pure-neutral regression requires 3-variable starts and bounds for m, s, and e.");
    }

    vector<vector<double>> starts = candidate_starts_cpp(x0, lower, upper, multistart);
    bool have_result = false;
    PureNeutralRegressionResult best;
    for (const auto &start : starts) {
        PureNeutralRegressionResult candidate = solve_one_start_cpp(
            base_args,
            density_records,
            density_scale,
            pure_vle_records,
            pure_vle_scale,
            start,
            lower,
            upper,
            derivative_test
        );
        if (!have_result) {
            best = candidate;
            have_result = true;
            continue;
        }
        if (candidate.success && !best.success) {
            best = candidate;
            continue;
        }
        if (candidate.success == best.success && candidate.cost < best.cost) {
            best = candidate;
        }
    }
    if (!have_result) {
        throw ValueError("Native IPOPT pure-neutral regression did not generate any candidate starts.");
    }
    return best;
}
