#include "epcsaft_core_internal.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <sstream>
#include <string>

#include <coin/IpIpoptApplication.hpp>
#include <coin/IpSolveStatistics.hpp>
#include <coin/IpTNLP.hpp>
#include <unsupported/Eigen/LevenbergMarquardt>

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
constexpr double kPresolveAcceptDensityRms = 1.0e-2;
constexpr double kPresolveAcceptPureVleRms = 1.0e-2;
constexpr int kThetaSize = 3;
constexpr int kSeedRho = 0;
constexpr int kSeedM = 1;
constexpr int kSeedS = 2;
constexpr int kSeedE = 3;
constexpr int kSeedCount = 4;

using ParamDerivative = Eigen::VectorXd;
using ParamDual = Eigen::AutoDiffScalar<ParamDerivative>;
using ResidualVector = Eigen::VectorXd;
using ResidualJacobian = Eigen::MatrixXd;
using LMInputVector = Eigen::VectorXd;

ParamDual make_param_dual(double value, int seed_index = -1);

template <typename Scalar>
Scalar regression_scalar_constant(double value) {
    return scalar_constant<Scalar>(value);
}

template <>
AutoDual regression_scalar_constant<AutoDual>(double value) {
    return make_autodiff_scalar(value, 0.0);
}

template <>
ParamDual regression_scalar_constant<ParamDual>(double value) {
    return make_param_dual(value);
}

template <typename Scalar>
struct PureNeutralStateScalar {
    Scalar den = regression_scalar_constant<Scalar>(0.0);
    Scalar d = regression_scalar_constant<Scalar>(0.0);
    Scalar zeta0 = regression_scalar_constant<Scalar>(0.0);
    Scalar zeta1 = regression_scalar_constant<Scalar>(0.0);
    Scalar zeta2 = regression_scalar_constant<Scalar>(0.0);
    Scalar zeta3 = regression_scalar_constant<Scalar>(0.0);
    Scalar eta = regression_scalar_constant<Scalar>(0.0);
    Scalar pair_diameter = regression_scalar_constant<Scalar>(0.0);
    Scalar ghs = regression_scalar_constant<Scalar>(0.0);
    Scalar ares_hs = regression_scalar_constant<Scalar>(0.0);
    Scalar ares_hc = regression_scalar_constant<Scalar>(0.0);
    Scalar I1 = regression_scalar_constant<Scalar>(0.0);
    Scalar I2 = regression_scalar_constant<Scalar>(0.0);
    Scalar dEtaI1_deta = regression_scalar_constant<Scalar>(0.0);
    Scalar dEtaI2_deta = regression_scalar_constant<Scalar>(0.0);
    Scalar C1 = regression_scalar_constant<Scalar>(0.0);
    Scalar C2 = regression_scalar_constant<Scalar>(0.0);
    Scalar m2es3 = regression_scalar_constant<Scalar>(0.0);
    Scalar m2e2s3 = regression_scalar_constant<Scalar>(0.0);
    Scalar ares_disp = regression_scalar_constant<Scalar>(0.0);
    Scalar ares_total = regression_scalar_constant<Scalar>(0.0);
    Scalar zraw_hc = regression_scalar_constant<Scalar>(0.0);
    Scalar zraw_disp = regression_scalar_constant<Scalar>(0.0);
    Scalar zraw_total = regression_scalar_constant<Scalar>(0.0);
    Scalar Z = regression_scalar_constant<Scalar>(0.0);
    Scalar pressure = regression_scalar_constant<Scalar>(0.0);
    Scalar lnfug = regression_scalar_constant<Scalar>(0.0);
};

struct RegressionProfilingStats {
    int objective_evaluations = 0;
    int gradient_evaluations = 0;
    int residual_evaluations = 0;
    int density_solves = 0;
    int fused_state_evaluations = 0;
    double callback_wall_time_s = 0.0;
};

struct PureNeutralFusedState {
    double pressure = 0.0;
    double lnfug = 0.0;
    double Z = 0.0;
    double dpdrho = 0.0;
    double dlnfug_drho = 0.0;
    std::array<double, kThetaSize> dpdtheta{0.0, 0.0, 0.0};
    std::array<double, kThetaSize> dlnfugdtheta{0.0, 0.0, 0.0};
};

struct PureNeutralResidualEvaluation {
    ResidualVector residuals;
    ResidualJacobian jacobian;
    vector<double> density_raw_residuals;
    vector<double> pure_vle_raw_residuals;
};

struct PureNeutralObjectiveEvaluation {
    double objective = 0.0;
    std::array<double, kThetaSize> gradient{0.0, 0.0, 0.0};
    vector<double> density_raw_residuals;
    vector<double> pure_vle_raw_residuals;
};

ParamDual make_param_dual(double value, int seed_index) {
    ParamDual x;
    x.value() = value;
    x.derivatives() = ParamDerivative::Zero(kThetaSize);
    if (seed_index >= 0) {
        x.derivatives()[seed_index] = 1.0;
    }
    return x;
}

double scalar_derivative_at(double, int) {
    return 0.0;
}

double scalar_derivative_at(const AutoDual &x, int idx) {
    return idx == 0 ? scalar_derivative(x) : 0.0;
}

double scalar_derivative_at(const ParamDual &x, int idx) {
    return (idx < x.derivatives().size()) ? x.derivatives()[idx] : 0.0;
}

void validate_pure_neutral_base_args_cpp(const add_args &base_args) {
    if (base_args.m.size() != 1 || base_args.s.size() != 1 || base_args.e.size() != 1) {
        throw ValueError("Native pure-neutral regression requires exactly one component.");
    }
    if (!base_args.z.empty() && (base_args.z.size() != 1 || std::abs(base_args.z[0]) > 1.0e-12)) {
        throw ValueError("Native pure-neutral regression currently supports only neutral single-component models.");
    }
    if (!base_args.assoc_num.empty() || !base_args.assoc_matrix.empty() || !base_args.k_hb.empty() || !base_args.e_assoc.empty() || !base_args.vol_a.empty()) {
        throw ValueError("Native pure-neutral regression currently supports only nonassociating single-component models.");
    }
    if (base_args.mw.size() != 1) {
        throw ValueError("Native pure-neutral regression requires a single MW value in the fixed parameter payload.");
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
    const Scalar s2 = s * s;
    const Scalar s3 = s2 * s;
    state.den = rho * (N_AV / 1.0e30);
    state.d = s * (1.0 - 0.12 * scalar_exp(-3.0 * e / t));
    const Scalar d2 = state.d * state.d;
    const Scalar d3 = d2 * state.d;

    const Scalar prefactor = PI / 6.0 * state.den * m;
    state.zeta0 = prefactor;
    state.zeta1 = prefactor * state.d;
    state.zeta2 = prefactor * d2;
    state.zeta3 = prefactor * d3;
    state.eta = state.zeta3;
    state.pair_diameter = state.d / 2.0;
    const Scalar pair_diameter2 = state.pair_diameter * state.pair_diameter;
    const Scalar zeta2_sq = state.zeta2 * state.zeta2;
    const Scalar zeta2_cu = zeta2_sq * state.zeta2;
    const Scalar zeta3_sq = state.zeta3 * state.zeta3;
    state.ghs = 1.0 / (1.0 - state.zeta3)
        + state.pair_diameter * 3.0 * state.zeta2 / scalar_pow(1.0 - state.zeta3, 2.0)
        + pair_diameter2 * 2.0 * zeta2_sq / scalar_pow(1.0 - state.zeta3, 3.0);

    state.ares_hs = 1.0 / state.zeta0 * (
        3.0 * state.zeta1 * state.zeta2 / (1.0 - state.zeta3)
        + zeta2_cu / (state.zeta3 * scalar_pow(1.0 - state.zeta3, 2.0))
        + (zeta2_cu / zeta3_sq - state.zeta0) * scalar_log(1.0 - state.zeta3)
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

    const Scalar e_over_t = e / t;
    const Scalar e_over_t2 = e_over_t * e_over_t;
    state.m2es3 = m * m * e_over_t * s3;
    state.m2e2s3 = m * m * e_over_t2 * s3;
    state.ares_disp = -2.0 * PI * state.den * state.I1 * state.m2es3
        - PI * state.den * m * state.C1 * state.I2 * state.m2e2s3;
    state.ares_total = state.ares_hc + state.ares_disp;

    const Scalar dghs_drho = state.zeta3 / scalar_pow(1.0 - state.zeta3, 2.0)
        + state.pair_diameter * (
            3.0 * state.zeta2 / scalar_pow(1.0 - state.zeta3, 2.0)
            + 6.0 * state.zeta2 * state.zeta3 / scalar_pow(1.0 - state.zeta3, 3.0)
        )
        + pair_diameter2 * (
            4.0 * zeta2_sq / scalar_pow(1.0 - state.zeta3, 3.0)
            + 6.0 * zeta2_sq * state.zeta3 / scalar_pow(1.0 - state.zeta3, 4.0)
        );
    const Scalar dadrho_hs = state.zeta3 / (1.0 - state.zeta3)
        + 3.0 * state.zeta1 * state.zeta2 / state.zeta0 / scalar_pow(1.0 - state.zeta3, 2.0)
        + (3.0 * zeta2_cu - state.zeta3 * zeta2_cu)
            / state.zeta0 / scalar_pow(1.0 - state.zeta3, 3.0);
    state.zraw_hc = m * dadrho_hs - (m - 1.0) * dghs_drho / state.ghs;
    state.zraw_disp = -2.0 * PI * state.den * state.dEtaI1_deta * state.m2es3
        - PI * state.den * m * (state.C1 * state.dEtaI2_deta + state.C2 * state.eta * state.I2) * state.m2e2s3;
    state.zraw_total = state.zraw_hc + state.zraw_disp;
    state.Z = 1.0 + state.zraw_total;
    if (!(scalar_value(state.Z) > 0.0)) {
        throw ValueError("Encountered non-positive compressibility factor during native regression evaluation.");
    }
    state.pressure = state.Z * kb * t * state.den * 1.0e30;
    state.lnfug = state.ares_total + state.zraw_total - scalar_log(state.Z);
    return state;
}

PureNeutralFusedState evaluate_fused_state_cpp(double t, double rho, const vector<double> &x) {
    AutoDual rho_dual = make_autodiff_scalar(rho, 1.0);
    auto rho_state = pure_neutral_state_scalar_cpp<AutoDual>(t, rho_dual, x[0], x[1], x[2]);
    constexpr std::array<double, kThetaSize> kParameterSteps = {1.0e-6, 1.0e-6, 1.0e-5};

    PureNeutralFusedState out;
    out.pressure = scalar_value(rho_state.pressure);
    out.lnfug = scalar_value(rho_state.lnfug);
    out.Z = scalar_value(rho_state.Z);
    out.dpdrho = scalar_derivative_at(rho_state.pressure, 0);
    out.dlnfug_drho = scalar_derivative_at(rho_state.lnfug, 0);
    for (int j = 0; j < kThetaSize; ++j) {
        vector<double> x_forward = x;
        vector<double> x_backward = x;
        x_forward[static_cast<size_t>(j)] += kParameterSteps[static_cast<size_t>(j)];
        x_backward[static_cast<size_t>(j)] -= kParameterSteps[static_cast<size_t>(j)];
        auto forward_state = pure_neutral_state_scalar_cpp<double>(t, rho, x_forward[0], x_forward[1], x_forward[2]);
        auto backward_state = pure_neutral_state_scalar_cpp<double>(t, rho, x_backward[0], x_backward[1], x_backward[2]);
        out.dpdtheta[static_cast<size_t>(j)] =
            (forward_state.pressure - backward_state.pressure) / (2.0 * kParameterSteps[static_cast<size_t>(j)]);
        out.dlnfugdtheta[static_cast<size_t>(j)] =
            (forward_state.lnfug - backward_state.lnfug) / (2.0 * kParameterSteps[static_cast<size_t>(j)]);
    }
    return out;
}

add_args pure_neutral_args_with_theta_cpp(const add_args &base_args, const vector<double> &x) {
    if (x.size() != kThetaSize) {
        throw ValueError("Native pure-neutral regression expects exactly three optimization variables.");
    }
    add_args args = base_args;
    args.m[0] = x[0];
    args.s[0] = x[1];
    args.e[0] = x[2];
    return args;
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

bool same_start_cpp(const vector<double> &a, const vector<double> &b) {
    if (a.size() != b.size()) {
        return false;
    }
    for (size_t i = 0; i < a.size(); ++i) {
        double scale = std::max({1.0, std::abs(a[i]), std::abs(b[i])});
        if (std::abs(a[i] - b[i]) > 1.0e-12 * scale) {
            return false;
        }
    }
    return true;
}

void append_start_if_distinct_cpp(vector<vector<double>> &starts, vector<double> point) {
    for (const auto &existing : starts) {
        if (same_start_cpp(existing, point)) {
            return;
        }
    }
    starts.push_back(std::move(point));
}

vector<vector<double>> candidate_starts_cpp(
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart
) {
    if (x0.size() != lower.size() || x0.size() != upper.size()) {
        throw ValueError("Initial guess and bounds must have matching lengths for native regression.");
    }
    vector<vector<double>> starts;
    vector<double> first = x0;
    for (size_t i = 0; i < first.size(); ++i) {
        first[i] = clip_start_value_cpp(first[i], lower[i], upper[i]);
    }
    append_start_if_distinct_cpp(starts, first);

    constexpr std::array<double, kThetaSize> kCounterBiasFactors = {0.93, 1.04, 0.95};
    vector<double> counter_bias = first;
    for (size_t i = 0; i < counter_bias.size(); ++i) {
        counter_bias[i] = clip_start_value_cpp(
            first[i] * kCounterBiasFactors[i],
            lower[i],
            upper[i]
        );
    }
    append_start_if_distinct_cpp(starts, std::move(counter_bias));

    std::mt19937 rng(12345);
    for (int k = 0; k < multistart; ++k) {
        vector<double> point(x0.size(), 0.0);
        for (size_t i = 0; i < x0.size(); ++i) {
            std::uniform_real_distribution<double> uniform(lower[i], upper[i]);
            point[i] = clip_start_value_cpp(uniform(rng), lower[i], upper[i]);
        }
        append_start_if_distinct_cpp(starts, std::move(point));
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
}

std::string least_squares_status_message_cpp(Eigen::LevenbergMarquardtSpace::Status status) {
    using Eigen::LevenbergMarquardtSpace::Status;
    switch (status) {
        case Status::RelativeReductionTooSmall:
            return "Levenberg-Marquardt terminated because the relative reduction became small.";
        case Status::RelativeErrorTooSmall:
            return "Levenberg-Marquardt terminated because the relative step became small.";
        case Status::RelativeErrorAndReductionTooSmall:
            return "Levenberg-Marquardt terminated because the relative step and reduction became small.";
        case Status::CosinusTooSmall:
            return "Levenberg-Marquardt terminated because the gradient cosine became small.";
        case Status::TooManyFunctionEvaluation:
            return "Levenberg-Marquardt exceeded the maximum function evaluation count.";
        case Status::FtolTooSmall:
            return "Levenberg-Marquardt ftol is too small for further progress.";
        case Status::XtolTooSmall:
            return "Levenberg-Marquardt xtol is too small for further progress.";
        case Status::GtolTooSmall:
            return "Levenberg-Marquardt gtol is too small for further progress.";
        case Status::ImproperInputParameters:
            return "Levenberg-Marquardt reported improper input parameters.";
        case Status::UserAsked:
            return "Levenberg-Marquardt terminated early by user request.";
        case Status::Running:
        case Status::NotStarted:
        default:
            return "Levenberg-Marquardt terminated without a clear success code.";
    }
}

bool least_squares_status_success_cpp(Eigen::LevenbergMarquardtSpace::Status status) {
    using Eigen::LevenbergMarquardtSpace::Status;
    return status == Status::RelativeReductionTooSmall
        || status == Status::RelativeErrorTooSmall
        || status == Status::RelativeErrorAndReductionTooSmall
        || status == Status::CosinusTooSmall;
}

void validate_fused_state_cpp(const PureNeutralFusedState &state, const char *label) {
    if (!(std::isfinite(state.dpdrho) && std::abs(state.dpdrho) > 0.0)) {
        throw ValueError(std::string("Encountered invalid exact dp/drho during native regression ") + label + " evaluation.");
    }
    if (!(std::isfinite(state.dlnfug_drho) && std::isfinite(state.pressure) && std::isfinite(state.lnfug) && std::isfinite(state.Z) && state.Z > 0.0)) {
        throw ValueError(std::string("Encountered invalid fused state during native regression ") + label + " evaluation.");
    }
}

PureNeutralResidualEvaluation evaluate_residual_jacobian_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x,
    RegressionProfilingStats *stats
) {
    auto callback_start = std::chrono::steady_clock::now();
    if (stats != nullptr) {
        ++stats->residual_evaluations;
    }

    PureNeutralResidualEvaluation eval;
    const vector<double> one_x = {1.0};
    add_args args = pure_neutral_args_with_theta_cpp(base_args, x);

    const Index density_count = static_cast<Index>(density_records.size());
    const Index pure_vle_count = static_cast<Index>(pure_vle_records.size());
    const Index total_count = density_count + pure_vle_count;
    eval.residuals = ResidualVector::Zero(total_count);
    eval.jacobian = ResidualJacobian::Zero(total_count, kThetaSize);
    eval.density_raw_residuals.reserve(density_records.size());
    eval.pure_vle_raw_residuals.reserve(pure_vle_records.size());

    Index row = 0;
    for (const auto &record : density_records) {
        double rho_calc = den_cpp(record.t, record.p, one_x, record.phase, args);
        if (stats != nullptr) {
            ++stats->density_solves;
        }
        double denom = std::max(std::abs(record.rho_exp), kRegressionGradientFloor);
        double raw = (rho_calc - record.rho_exp) / denom;
        eval.density_raw_residuals.push_back(raw);
        eval.residuals[row] = density_scale * raw;

        PureNeutralFusedState state = evaluate_fused_state_cpp(record.t, rho_calc, x);
        if (stats != nullptr) {
            ++stats->fused_state_evaluations;
        }
        validate_fused_state_cpp(state, "density residual");
        for (int j = 0; j < kThetaSize; ++j) {
            double drho_dtheta = -state.dpdtheta[static_cast<size_t>(j)] / state.dpdrho;
            eval.jacobian(row, j) = density_scale * (drho_dtheta / denom);
        }
        ++row;
    }

    for (const auto &record : pure_vle_records) {
        double rho_liq = den_cpp(record.t, record.p, one_x, 0, args);
        double rho_vap = den_cpp(record.t, record.p, one_x, 1, args);
        if (stats != nullptr) {
            stats->density_solves += 2;
        }

        PureNeutralFusedState liq = evaluate_fused_state_cpp(record.t, rho_liq, x);
        PureNeutralFusedState vap = evaluate_fused_state_cpp(record.t, rho_vap, x);
        if (stats != nullptr) {
            stats->fused_state_evaluations += 2;
        }
        validate_fused_state_cpp(liq, "liquid fugacity-balance");
        validate_fused_state_cpp(vap, "vapor fugacity-balance");

        double raw = liq.lnfug - vap.lnfug;
        eval.pure_vle_raw_residuals.push_back(raw);
        eval.residuals[row] = pure_vle_scale * raw;

        for (int j = 0; j < kThetaSize; ++j) {
            double drho_liq_dtheta = -liq.dpdtheta[static_cast<size_t>(j)] / liq.dpdrho;
            double drho_vap_dtheta = -vap.dpdtheta[static_cast<size_t>(j)] / vap.dpdrho;
            double total_grad =
                liq.dlnfugdtheta[static_cast<size_t>(j)] + liq.dlnfug_drho * drho_liq_dtheta
                - vap.dlnfugdtheta[static_cast<size_t>(j)] - vap.dlnfug_drho * drho_vap_dtheta;
            eval.jacobian(row, j) = pure_vle_scale * total_grad;
        }
        ++row;
    }

    if (stats != nullptr) {
        stats->callback_wall_time_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - callback_start).count();
    }
    return eval;
}

PureNeutralObjectiveEvaluation objective_from_residual_eval_cpp(const PureNeutralResidualEvaluation &residual_eval) {
    PureNeutralObjectiveEvaluation out;
    out.objective = 0.5 * residual_eval.residuals.squaredNorm();
    Eigen::Matrix<double, kThetaSize, 1> gradient = residual_eval.jacobian.transpose() * residual_eval.residuals;
    for (int j = 0; j < kThetaSize; ++j) {
        out.gradient[static_cast<size_t>(j)] = gradient[j];
    }
    out.density_raw_residuals = residual_eval.density_raw_residuals;
    out.pure_vle_raw_residuals = residual_eval.pure_vle_raw_residuals;
    return out;
}

PureNeutralObjectiveEvaluation evaluate_pure_neutral_objective_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x,
    RegressionProfilingStats *stats
) {
    PureNeutralResidualEvaluation residual_eval = evaluate_residual_jacobian_cpp(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        x,
        stats
    );
    return objective_from_residual_eval_cpp(residual_eval);
}

struct BoundedTransformResult {
    vector<double> x;
    std::array<double, kThetaSize> dxdy{1.0, 1.0, 1.0};
};

double logistic_cpp(double y) {
    if (y >= 0.0) {
        double exp_neg = std::exp(-y);
        return 1.0 / (1.0 + exp_neg);
    }
    double exp_pos = std::exp(y);
    return exp_pos / (1.0 + exp_pos);
}

double logit_cpp(double p) {
    return std::log(p / (1.0 - p));
}

BoundedTransformResult unconstrained_to_bounded_cpp(
    const LMInputVector &y,
    const vector<double> &lower,
    const vector<double> &upper
) {
    BoundedTransformResult out;
    out.x.resize(kThetaSize, 0.0);
    for (int i = 0; i < kThetaSize; ++i) {
        double lo = lower[static_cast<size_t>(i)];
        double hi = upper[static_cast<size_t>(i)];
        double sigma = logistic_cpp(y[i]);
        double width = hi - lo;
        out.x[static_cast<size_t>(i)] = lo + width * sigma;
        out.dxdy[static_cast<size_t>(i)] = width * sigma * (1.0 - sigma);
    }
    return out;
}

LMInputVector bounded_to_unconstrained_cpp(
    const vector<double> &x,
    const vector<double> &lower,
    const vector<double> &upper
) {
    LMInputVector y(kThetaSize);
    for (int i = 0; i < kThetaSize; ++i) {
        double lo = lower[static_cast<size_t>(i)];
        double hi = upper[static_cast<size_t>(i)];
        double clipped = clip_start_value_cpp(x[static_cast<size_t>(i)], lo, hi);
        double width = hi - lo;
        double p = (clipped - lo) / width;
        if (p < 1.0e-12) {
            p = 1.0e-12;
        } else if (p > 1.0 - 1.0e-12) {
            p = 1.0 - 1.0e-12;
        }
        y[i] = logit_cpp(p);
    }
    return y;
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
        n = kThetaSize;
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
        ++profiling_.objective_evaluations;
        if (!ensure_cache(n, x, new_x)) {
            return false;
        }
        obj_value = cache_.objective;
        return true;
    }

    bool eval_grad_f(Index n, const Number *x, bool new_x, Number *grad_f) override {
        ++profiling_.gradient_evaluations;
        if (!ensure_cache(n, x, new_x)) {
            return false;
        }
        for (Index i = 0; i < n; ++i) {
            grad_f[i] = cache_.gradient[static_cast<size_t>(i)];
        }
        return true;
    }

    bool eval_g(Index, const Number *, bool, Index, Number *) override {
        return true;
    }

    bool eval_jac_g(Index, const Number *, bool, Index, Index, Index *, Index *, Number *) override {
        return true;
    }

    bool get_list_of_nonlinear_variables(Index num_nonlin_vars, Index *pos_nonlin_vars) override {
        if (num_nonlin_vars != kThetaSize) {
            return false;
        }
        for (int i = 0; i < kThetaSize; ++i) {
            pos_nonlin_vars[i] = i;
        }
        return true;
    }

    Index get_number_of_nonlinear_variables() override {
        return kThetaSize;
    }

    bool intermediate_callback(
        AlgorithmMode,
        Index iter,
        Number,
        Number,
        Number,
        Number,
        Number,
        Number,
        Number,
        Number,
        Index,
        const Ipopt::IpoptData *,
        Ipopt::IpoptCalculatedQuantities *
    ) override {
        iterations_ = static_cast<int>(iter);
        return true;
    }

    void finalize_solution(
        SolverReturn status,
        Index n,
        const Number *x,
        const Number *,
        const Number *,
        Index,
        const Number *,
        const Number *,
        Number obj_value,
        const Ipopt::IpoptData *,
        Ipopt::IpoptCalculatedQuantities *
    ) override {
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
                    solution_,
                    &profiling_
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
        out.nfev = profiling_.objective_evaluations + profiling_.gradient_evaluations;
        out.iterations = iterations_;
        out.objective_evaluations = profiling_.objective_evaluations;
        out.gradient_evaluations = profiling_.gradient_evaluations;
        out.residual_evaluations = profiling_.residual_evaluations;
        out.density_solves = profiling_.density_solves;
        out.fused_state_evaluations = profiling_.fused_state_evaluations;
        out.callback_wall_time_s = profiling_.callback_wall_time_s;
        out.backend = "ipopt_native";
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
                    current,
                    &profiling_
                );
            } catch (...) {
                return false;
            }
            cache_valid_ = true;
            last_x_ = std::move(current);
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
    int iterations_ = 0;
    SolverReturn status_ = Ipopt::INTERNAL_ERROR;
    double obj_value_ = HUGE_DBL;
    RegressionProfilingStats profiling_;
};

PureNeutralRegressionResult solve_one_start_ipopt_cpp(
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
    auto solve_start = std::chrono::steady_clock::now();
    PureNeutralRegressionResult warm_start;
    bool have_warm_start = false;
    vector<double> seeded_start = start;
    try {
        warm_start = fit_pure_neutral_least_squares_cpp(
            base_args,
            density_records,
            density_scale,
            pure_vle_records,
            pure_vle_scale,
            start,
            lower,
            upper,
            0
        );
        if (warm_start.success && warm_start.x.size() == kThetaSize) {
            if (
                warm_start.density_metric < kPresolveAcceptDensityRms
                && warm_start.pure_vle_metric < kPresolveAcceptPureVleRms
            ) {
                warm_start.backend = "ipopt_native";
                warm_start.message = "Accepted native least-squares presolve solution.";
                return warm_start;
            }
            seeded_start = warm_start.x;
            have_warm_start = true;
        }
    } catch (...) {
    }

    SmartPtr<PureNeutralRegressionTNLP> nlp = new PureNeutralRegressionTNLP(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        seeded_start,
        lower,
        upper
    );
    SmartPtr<IpoptApplication> app = ::IpoptApplicationFactory();
    app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    app->Options()->SetStringValue("sb", "yes");
    app->Options()->SetIntegerValue("print_level", 0);
    app->Options()->SetIntegerValue("max_iter", 500);
    app->Options()->SetIntegerValue("acceptable_iter", 3);
    app->Options()->SetNumericValue("tol", 1.0e-6);
    app->Options()->SetNumericValue("acceptable_tol", 1.0e-5);
    if (derivative_test) {
        app->Options()->SetStringValue("derivative_test", "first-order");
        app->Options()->SetNumericValue("derivative_test_perturbation", 1.0e-7);
    }
    ApplicationReturnStatus init_status = app->Initialize();
    if (init_status != Ipopt::Solve_Succeeded) {
        throw ValueError("Failed to initialize IPOPT application for native pure-neutral regression.");
    }
    app->OptimizeTNLP(GetRawPtr(nlp));
    PureNeutralRegressionResult out = nlp->result();
    out.solve_wall_time_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start).count();
    if (have_warm_start) {
        out.nfev += warm_start.nfev;
        out.objective_evaluations += warm_start.objective_evaluations;
        out.gradient_evaluations += warm_start.gradient_evaluations;
        out.residual_evaluations += warm_start.residual_evaluations;
        out.density_solves += warm_start.density_solves;
        out.fused_state_evaluations += warm_start.fused_state_evaluations;
        out.callback_wall_time_s += warm_start.callback_wall_time_s;
    }
    return out;
}

struct PureNeutralLeastSquaresFunctor : Eigen::DenseFunctor<double> {
    PureNeutralLeastSquaresFunctor(
        add_args base_args,
        vector<PureNeutralRegressionDensityRecord> density_records,
        double density_scale,
        vector<PureNeutralRegressionVLERecord> pure_vle_records,
        double pure_vle_scale,
        vector<double> lower,
        vector<double> upper
    )
        : Eigen::DenseFunctor<double>(
              kThetaSize,
              static_cast<int>(density_records.size() + pure_vle_records.size())
          ),
          base_args(std::move(base_args)),
          density_records(std::move(density_records)),
          density_scale(density_scale),
          pure_vle_records(std::move(pure_vle_records)),
          pure_vle_scale(pure_vle_scale),
          lower(std::move(lower)),
          upper(std::move(upper)) {}

    int operator()(const LMInputVector &y, ResidualVector &fvec) {
        if (!ensure_cache(y)) {
            return -1;
        }
        fvec = cache_eval.residuals;
        return 0;
    }

    int df(const LMInputVector &y, ResidualJacobian &fjac) {
        if (!ensure_cache(y)) {
            return -1;
        }
        fjac = cache_eval.jacobian;
        for (int j = 0; j < kThetaSize; ++j) {
            fjac.col(j) *= cache_transform.dxdy[static_cast<size_t>(j)];
        }
        return 0;
    }

    RegressionProfilingStats profiling;
    BoundedTransformResult cache_transform;
    PureNeutralResidualEvaluation cache_eval;
    bool cache_valid = false;
    LMInputVector last_y = LMInputVector::Zero(kThetaSize);

private:
    bool ensure_cache(const LMInputVector &y) {
        if (!cache_valid || !y.isApprox(last_y, 0.0)) {
            try {
                cache_transform = unconstrained_to_bounded_cpp(y, lower, upper);
                cache_eval = evaluate_residual_jacobian_cpp(
                    base_args,
                    density_records,
                    density_scale,
                    pure_vle_records,
                    pure_vle_scale,
                    cache_transform.x,
                    &profiling
                );
            } catch (...) {
                return false;
            }
            cache_valid = true;
            last_y = y;
        }
        return true;
    }

    add_args base_args;
    vector<PureNeutralRegressionDensityRecord> density_records;
    double density_scale = 1.0;
    vector<PureNeutralRegressionVLERecord> pure_vle_records;
    double pure_vle_scale = 1.0;
    vector<double> lower;
    vector<double> upper;
};

PureNeutralRegressionResult solve_one_start_least_squares_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &start,
    const vector<double> &lower,
    const vector<double> &upper
) {
    auto solve_start = std::chrono::steady_clock::now();
    PureNeutralLeastSquaresFunctor functor(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        lower,
        upper
    );
    Eigen::LevenbergMarquardt<PureNeutralLeastSquaresFunctor> lm(functor);
    lm.setFtol(1.0e-6);
    lm.setXtol(1.0e-6);
    lm.setGtol(0.0);
    lm.setFactor(10.0);
    lm.setMaxfev(200);

    LMInputVector y = bounded_to_unconstrained_cpp(start, lower, upper);
    Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(y);
    BoundedTransformResult final_transform = unconstrained_to_bounded_cpp(y, lower, upper);

    PureNeutralRegressionResult out;
    out.x = final_transform.x;
    out.success = least_squares_status_success_cpp(status);
    out.status = static_cast<int>(status);
    out.message = least_squares_status_message_cpp(status);
    out.nfev = static_cast<int>(lm.nfev() + lm.njev());
    out.iterations = static_cast<int>(lm.iterations());
    out.objective_evaluations = static_cast<int>(lm.nfev());
    out.gradient_evaluations = static_cast<int>(lm.njev());
    out.residual_evaluations = functor.profiling.residual_evaluations;
    out.density_solves = functor.profiling.density_solves;
    out.fused_state_evaluations = functor.profiling.fused_state_evaluations;
    out.callback_wall_time_s = functor.profiling.callback_wall_time_s;
    out.backend = "least_squares_native";
    out.solve_wall_time_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start).count();

    PureNeutralResidualEvaluation final_eval = evaluate_residual_jacobian_cpp(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        out.x,
        &functor.profiling
    );
    PureNeutralObjectiveEvaluation final_objective = objective_from_residual_eval_cpp(final_eval);
    out.cost = final_objective.objective;
    out.residual_norm = std::sqrt(std::max(0.0, 2.0 * final_objective.objective));
    out.density_metric = rms_metric_cpp(final_objective.density_raw_residuals);
    out.pure_vle_metric = rms_metric_cpp(final_objective.pure_vle_raw_residuals);
    out.residual_evaluations = functor.profiling.residual_evaluations;
    out.density_solves = functor.profiling.density_solves;
    out.fused_state_evaluations = functor.profiling.fused_state_evaluations;
    out.callback_wall_time_s = functor.profiling.callback_wall_time_s;
    return out;
}

PureNeutralRegressionResult choose_better_result_cpp(
    bool have_result,
    const PureNeutralRegressionResult &best,
    const PureNeutralRegressionResult &candidate
) {
    if (!have_result) {
        return candidate;
    }
    if (candidate.success && !best.success) {
        return candidate;
    }
    if (candidate.success == best.success && candidate.cost < best.cost) {
        return candidate;
    }
    return best;
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
    constexpr std::array<double, kThetaSize> kObjectiveSteps = {1.0e-6, 1.0e-6, 1.0e-5};
    RegressionProfilingStats profiling;
    PureNeutralResidualEvaluation residual_eval = evaluate_residual_jacobian_cpp(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        x,
        &profiling
    );
    PureNeutralObjectiveEvaluation eval = objective_from_residual_eval_cpp(residual_eval);
    PureNeutralRegressionDebugResult out;
    out.objective = eval.objective;
    out.gradient.assign(static_cast<size_t>(kThetaSize), 0.0);
    for (int j = 0; j < kThetaSize; ++j) {
        vector<double> x_forward = x;
        vector<double> x_backward = x;
        x_forward[static_cast<size_t>(j)] += kObjectiveSteps[static_cast<size_t>(j)];
        x_backward[static_cast<size_t>(j)] -= kObjectiveSteps[static_cast<size_t>(j)];
        PureNeutralObjectiveEvaluation forward_eval = evaluate_pure_neutral_objective_cpp(
            base_args,
            density_records,
            density_scale,
            pure_vle_records,
            pure_vle_scale,
            x_forward,
            nullptr
        );
        PureNeutralObjectiveEvaluation backward_eval = evaluate_pure_neutral_objective_cpp(
            base_args,
            density_records,
            density_scale,
            pure_vle_records,
            pure_vle_scale,
            x_backward,
            nullptr
        );
        out.gradient[static_cast<size_t>(j)] =
            (forward_eval.objective - backward_eval.objective) / (2.0 * kObjectiveSteps[static_cast<size_t>(j)]);
    }
    out.residuals.assign(
        residual_eval.residuals.data(),
        residual_eval.residuals.data() + residual_eval.residuals.size()
    );
    out.jacobian_rows = static_cast<int>(residual_eval.jacobian.rows());
    out.jacobian_cols = static_cast<int>(residual_eval.jacobian.cols());
    out.jacobian_row_major.resize(static_cast<size_t>(residual_eval.jacobian.size()));
    for (Eigen::Index i = 0; i < residual_eval.jacobian.rows(); ++i) {
        for (Eigen::Index j = 0; j < residual_eval.jacobian.cols(); ++j) {
            out.jacobian_row_major[static_cast<size_t>(i * residual_eval.jacobian.cols() + j)] = residual_eval.jacobian(i, j);
        }
    }
    out.density_raw_residuals = std::move(eval.density_raw_residuals);
    out.pure_vle_raw_residuals = std::move(eval.pure_vle_raw_residuals);
    out.residual_evaluations = profiling.residual_evaluations;
    out.density_solves = profiling.density_solves;
    out.fused_state_evaluations = profiling.fused_state_evaluations;
    out.callback_wall_time_s = profiling.callback_wall_time_s;
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
        throw ValueError("Native pure-neutral regression requires both density and pure-VLE record families.");
    }
    if (x0.size() != kThetaSize || lower.size() != kThetaSize || upper.size() != kThetaSize) {
        throw ValueError("Native pure-neutral regression requires 3-variable starts and bounds for m, s, and e.");
    }

    vector<vector<double>> starts = candidate_starts_cpp(x0, lower, upper, multistart);
    bool have_result = false;
    PureNeutralRegressionResult best;
    for (const auto &start : starts) {
        try {
            PureNeutralRegressionResult candidate = solve_one_start_ipopt_cpp(
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
            best = choose_better_result_cpp(have_result, best, candidate);
            have_result = true;
        } catch (...) {
        }
    }
    if (!have_result) {
        throw ValueError("Native pure-neutral IPOPT regression did not generate any candidate starts.");
    }
    return best;
}

PureNeutralRegressionResult fit_pure_neutral_least_squares_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart
) {
    validate_pure_neutral_base_args_cpp(base_args);
    if (density_records.empty() || pure_vle_records.empty()) {
        throw ValueError("Native pure-neutral regression requires both density and pure-VLE record families.");
    }
    if (x0.size() != kThetaSize || lower.size() != kThetaSize || upper.size() != kThetaSize) {
        throw ValueError("Native pure-neutral regression requires 3-variable starts and bounds for m, s, and e.");
    }

    vector<vector<double>> starts = candidate_starts_cpp(x0, lower, upper, multistart);
    bool have_result = false;
    PureNeutralRegressionResult best;
    for (const auto &start : starts) {
        try {
            PureNeutralRegressionResult candidate = solve_one_start_least_squares_cpp(
                base_args,
                density_records,
                density_scale,
                pure_vle_records,
                pure_vle_scale,
                start,
                lower,
                upper
            );
            best = choose_better_result_cpp(have_result, best, candidate);
            have_result = true;
        } catch (...) {
        }
    }
    if (!have_result) {
        throw ValueError("Native pure-neutral least-squares regression did not generate any candidate starts.");
    }
    return best;
}
