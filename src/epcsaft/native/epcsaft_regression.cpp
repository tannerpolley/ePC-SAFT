#include "epcsaft_core_internal.h"
#include "autodiff/ad_scalar.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <sstream>
#include <string>

#include <unsupported/Eigen/LevenbergMarquardt>

using Index = Eigen::Index;
using thermo_detail::kDispersionA0;
using thermo_detail::kDispersionA1;
using thermo_detail::kDispersionA2;
using thermo_detail::kDispersionB0;
using thermo_detail::kDispersionB1;
using thermo_detail::kDispersionB2;

namespace regression_detail {

constexpr double kRegressionGradientFloor = 1.0e-300;
constexpr double kTransformFiniteBound = 25.0;
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

double scalar_value(double x) {
    return x;
}

template <typename DerType>
double scalar_value(const Eigen::AutoDiffScalar<DerType> &x) {
    return x.value();
}

#ifdef EPCSAFT_HAS_CPPAD
template <typename Base>
double scalar_value(const CppAD::AD<Base> &x) {
    return CppAD::Value(x);
}
#endif

double scalar_log(double x) {
    return std::log(x);
}

template <typename DerType>
auto scalar_log(const Eigen::AutoDiffScalar<DerType> &x) -> decltype(log(x)) {
    using std::log;
    return log(x);
}

#ifdef EPCSAFT_HAS_CPPAD
template <typename Base>
CppAD::AD<Base> scalar_log(const CppAD::AD<Base> &x) {
    return CppAD::log(x);
}
#endif

double scalar_exp(double x) {
    return std::exp(x);
}

template <typename DerType>
auto scalar_exp(const Eigen::AutoDiffScalar<DerType> &x) -> decltype(exp(x)) {
    using std::exp;
    return exp(x);
}

#ifdef EPCSAFT_HAS_CPPAD
template <typename Base>
CppAD::AD<Base> scalar_exp(const CppAD::AD<Base> &x) {
    return CppAD::exp(x);
}
#endif

double scalar_pow(double x, int exponent) {
    return std::pow(x, exponent);
}

template <typename DerType>
auto scalar_pow(const Eigen::AutoDiffScalar<DerType> &x, int exponent) -> decltype(pow(x, static_cast<double>(exponent))) {
    using std::pow;
    return pow(x, static_cast<double>(exponent));
}

#ifdef EPCSAFT_HAS_CPPAD
template <typename Base>
CppAD::AD<Base> scalar_pow(const CppAD::AD<Base> &x, int exponent) {
    return CppAD::pow(x, static_cast<double>(exponent));
}
#endif

double scalar_pow(double x, double exponent) {
    return std::pow(x, exponent);
}

template <typename DerType>
auto scalar_pow(const Eigen::AutoDiffScalar<DerType> &x, double exponent) -> decltype(pow(x, exponent)) {
    using std::pow;
    return pow(x, exponent);
}

#ifdef EPCSAFT_HAS_CPPAD
template <typename Base>
CppAD::AD<Base> scalar_pow(const CppAD::AD<Base> &x, double exponent) {
    return CppAD::pow(x, exponent);
}
#endif

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
    PureNeutralFusedState out;
#ifdef EPCSAFT_HAS_CPPAD
    using epcsaft::autodiff::CppADScalar;
    std::vector<CppADScalar> independent(1);
    independent[0] = rho;
    CppAD::Independent(independent);
    auto rho_state = pure_neutral_state_scalar_cpp<CppADScalar>(t, independent[0], x[0], x[1], x[2]);
    std::vector<CppADScalar> dependent(2);
    dependent[0] = rho_state.pressure;
    dependent[1] = rho_state.lnfug;
    CppAD::ADFun<double> tape(independent, dependent);
    const std::vector<double> jacobian = tape.Jacobian(std::vector<double>{rho});
    out.pressure = scalar_value(rho_state.pressure);
    out.lnfug = scalar_value(rho_state.lnfug);
    out.Z = scalar_value(rho_state.Z);
    out.dpdrho = jacobian[0];
    out.dlnfug_drho = jacobian[1];
#else
    AutoDual rho_dual = make_autodiff_scalar(rho, 1.0);
    auto rho_state = pure_neutral_state_scalar_cpp<AutoDual>(t, rho_dual, x[0], x[1], x[2]);
    out.pressure = scalar_value(rho_state.pressure);
    out.lnfug = scalar_value(rho_state.lnfug);
    out.Z = scalar_value(rho_state.Z);
    out.dpdrho = scalar_derivative_at(rho_state.pressure, 0);
    out.dlnfug_drho = scalar_derivative_at(rho_state.lnfug, 0);
#endif
    auto theta_state = pure_neutral_state_scalar_cpp<ParamDual>(
        t,
        make_param_dual(rho),
        make_param_dual(x[0], 0),
        make_param_dual(x[1], 1),
        make_param_dual(x[2], 2)
    );
    for (int j = 0; j < kThetaSize; ++j) {
        out.dpdtheta[static_cast<size_t>(j)] = scalar_derivative_at(theta_state.pressure, j);
        out.dlnfugdtheta[static_cast<size_t>(j)] = scalar_derivative_at(theta_state.lnfug, j);
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

constexpr std::array<std::array<double, kThetaSize>, 5> kDeterministicSeedFactors{{
    {{1.00, 1.00, 1.00}},
    {{0.93, 1.04, 0.95}},
    {{1.07, 0.96, 1.05}},
    {{0.97, 1.02, 0.98}},
    {{1.03, 0.98, 1.02}},
}};

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

double validate_positive_bound_cpp(double value, const char *label) {
    if (!(value > 0.0) || !std::isfinite(value)) {
        throw ValueError(std::string("Native pure-neutral regression requires positive finite ") + label + " bounds.");
    }
    return value;
}

BoundedTransformResult unconstrained_to_bounded_cpp(
    const LMInputVector &y,
    const vector<double> &lower,
    const vector<double> &upper
) {
    BoundedTransformResult out;
    out.x.resize(kThetaSize, 0.0);
    for (int i = 0; i < kThetaSize; ++i) {
        double lo = validate_positive_bound_cpp(lower[static_cast<size_t>(i)], "lower");
        double hi = validate_positive_bound_cpp(upper[static_cast<size_t>(i)], "upper");
        if (!(hi > lo)) {
            throw ValueError("Native pure-neutral regression requires strictly increasing bounds for transformed variables.");
        }
        double log_lo = std::log(lo);
        double log_hi = std::log(hi);
        double sigma = logistic_cpp(y[i]);
        double log_x = log_lo + (log_hi - log_lo) * sigma;
        out.x[static_cast<size_t>(i)] = std::exp(log_x);
        out.dxdy[static_cast<size_t>(i)] = out.x[static_cast<size_t>(i)] * (log_hi - log_lo) * sigma * (1.0 - sigma);
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
        double lo = validate_positive_bound_cpp(lower[static_cast<size_t>(i)], "lower");
        double hi = validate_positive_bound_cpp(upper[static_cast<size_t>(i)], "upper");
        if (!(hi > lo)) {
            throw ValueError("Native pure-neutral regression requires strictly increasing bounds for transformed variables.");
        }
        double clipped = clip_start_value_cpp(x[static_cast<size_t>(i)], lo, hi);
        double log_lo = std::log(lo);
        double log_hi = std::log(hi);
        double p = (std::log(clipped) - log_lo) / (log_hi - log_lo);
        if (p < 1.0e-12) {
            p = 1.0e-12;
        } else if (p > 1.0 - 1.0e-12) {
            p = 1.0 - 1.0e-12;
        }
        y[i] = logit_cpp(p);
    }
    return y;
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
    for (const auto &factors : kDeterministicSeedFactors) {
        vector<double> point = first;
        for (int i = 0; i < kThetaSize; ++i) {
            point[static_cast<size_t>(i)] = clip_start_value_cpp(
                first[static_cast<size_t>(i)] * factors[static_cast<size_t>(i)],
                lower[static_cast<size_t>(i)],
                upper[static_cast<size_t>(i)]
            );
        }
        append_start_if_distinct_cpp(starts, std::move(point));
    }

    std::mt19937 rng(12345);
    for (int k = 0; k < multistart; ++k) {
        vector<double> point(kThetaSize, 0.0);
        for (int i = 0; i < kThetaSize; ++i) {
            double lo = validate_positive_bound_cpp(lower[static_cast<size_t>(i)], "lower");
            double hi = validate_positive_bound_cpp(upper[static_cast<size_t>(i)], "upper");
            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            double log_x = std::log(lo) + uniform(rng) * (std::log(hi) - std::log(lo));
            point[static_cast<size_t>(i)] = std::exp(log_x);
        }
        append_start_if_distinct_cpp(starts, std::move(point));
    }
    return starts;
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
    PureNeutralObjectiveEvaluation initial_eval = evaluate_pure_neutral_objective_cpp(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        start,
        nullptr
    );
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
    out.initial_cost = initial_eval.objective;
    out.initial_density_metric = rms_metric_cpp(initial_eval.density_raw_residuals);
    out.initial_pure_vle_metric = rms_metric_cpp(initial_eval.pure_vle_raw_residuals);
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

constexpr int kGenericTargetM = 0;
constexpr int kGenericTargetS = 1;
constexpr int kGenericTargetE = 2;
constexpr int kGenericTargetEAssoc = 3;
constexpr int kGenericTargetVolA = 4;
constexpr int kGenericTargetDBorn = 5;
constexpr int kGenericTargetKIJ = 6;
constexpr int kGenericTargetLIJ = 7;
constexpr int kGenericTargetKHB = 8;

constexpr int kGenericTermDensity = 1;
constexpr int kGenericTermPureVLE = 2;
constexpr int kGenericTermOsmotic = 3;
constexpr int kGenericTermMIAC = 4;
constexpr int kGenericTermBinaryVLE = 5;
constexpr int kGenericTermComponentLnFugacity = 6;

double relative_residual_cpp(double calc, double exp) {
    double denom = std::max(std::abs(exp), 1.0e-8);
    return (calc - exp) / denom;
}

double relative_or_absolute_residual_cpp(double calc, double exp) {
    if (std::abs(exp) <= 1.0e-8) {
        return calc - exp;
    }
    return (calc - exp) / std::abs(exp);
}

std::string generic_term_name_cpp(const GenericRegressionRecord &record) {
    if (!record.term_name.empty()) {
        return record.term_name;
    }
    switch (record.term) {
        case kGenericTermDensity:
            return "density";
        case kGenericTermPureVLE:
            return "pure_vle_fugacity_balance";
        case kGenericTermOsmotic:
            return "osmotic_coefficient";
        case kGenericTermMIAC:
            return "mean_ionic_activity";
        case kGenericTermBinaryVLE:
            return "binary_vle_fugacity_balance";
        case kGenericTermComponentLnFugacity:
            return "component_lnfugacity";
        default:
            return "unknown";
    }
}

void set_vector_value_cpp(vector<double> &values, int index, double value, const char *label) {
    if (index < 0 || static_cast<size_t>(index) >= values.size()) {
        throw ValueError(std::string("Native generic regression target index is out of range for ") + label + ".");
    }
    values[static_cast<size_t>(index)] = value;
}

void apply_generic_targets_cpp(
    add_args &args,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &theta
) {
    if (target_kinds.size() != theta.size() || target_indices.size() != theta.size() || target_indices_2.size() != theta.size()) {
        throw ValueError("Native generic regression target metadata must match theta length.");
    }
    const size_t n = args.m.size();
    for (size_t j = 0; j < theta.size(); ++j) {
        const double value = theta[j];
        const int index = target_indices[j];
        switch (target_kinds[j]) {
            case kGenericTargetM:
                set_vector_value_cpp(args.m, index, value, "m");
                break;
            case kGenericTargetS:
                set_vector_value_cpp(args.s, index, value, "s");
                break;
            case kGenericTargetE:
                set_vector_value_cpp(args.e, index, value, "e");
                break;
            case kGenericTargetEAssoc:
                set_vector_value_cpp(args.e_assoc, index, value, "e_assoc");
                break;
            case kGenericTargetVolA:
                set_vector_value_cpp(args.vol_a, index, value, "vol_a");
                break;
            case kGenericTargetDBorn:
                set_vector_value_cpp(args.d_born, index, value, "d_born");
                break;
            case kGenericTargetKIJ: {
                const int other = target_indices_2[j];
                if (index < 0 || other < 0 || static_cast<size_t>(index) >= n || static_cast<size_t>(other) >= n) {
                    throw ValueError("Native generic regression k_ij target index is out of range.");
                }
                if (args.k_ij.size() != n * n) {
                    throw ValueError("Native generic regression requires a dense k_ij matrix for binary targets.");
                }
                args.k_ij[static_cast<size_t>(index) * n + static_cast<size_t>(other)] = value;
                args.k_ij[static_cast<size_t>(other) * n + static_cast<size_t>(index)] = value;
                break;
            }
            case kGenericTargetLIJ: {
                const int other = target_indices_2[j];
                if (index < 0 || other < 0 || static_cast<size_t>(index) >= n || static_cast<size_t>(other) >= n) {
                    throw ValueError("Native generic regression l_ij target index is out of range.");
                }
                if (args.l_ij.size() != n * n) {
                    throw ValueError("Native generic regression requires a dense l_ij matrix for binary targets.");
                }
                args.l_ij[static_cast<size_t>(index) * n + static_cast<size_t>(other)] = value;
                args.l_ij[static_cast<size_t>(other) * n + static_cast<size_t>(index)] = value;
                break;
            }
            case kGenericTargetKHB: {
                const int other = target_indices_2[j];
                if (index < 0 || other < 0 || static_cast<size_t>(index) >= n || static_cast<size_t>(other) >= n) {
                    throw ValueError("Native generic regression k_hb_ij target index is out of range.");
                }
                if (args.k_hb.size() != n * n) {
                    throw ValueError("Native generic regression requires a dense k_hb matrix for binary targets.");
                }
                args.k_hb[static_cast<size_t>(index) * n + static_cast<size_t>(other)] = value;
                args.k_hb[static_cast<size_t>(other) * n + static_cast<size_t>(index)] = value;
                break;
            }
            default:
                throw ValueError("Unknown native generic regression target kind.");
        }
    }
}

double generic_mass_density_cpp(ePCSAFTStateNative &state, const add_args &args, const vector<double> &x) {
    if (args.mw.size() != x.size()) {
        throw ValueError("Native generic regression mass density requires one MW value per component.");
    }
    double mw_mix = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        mw_mix += x[i] * args.mw[i];
    }
    return state.density() * mw_mix;
}

void append_generic_residuals_cpp(
    const add_args &base_args,
    const GenericRegressionRecord &record,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &theta,
    vector<double> &scaled,
    std::map<std::string, vector<double>> &raw_by_term
) {
    add_args args = base_args;
    apply_generic_targets_cpp(args, target_kinds, target_indices, target_indices_2, theta);
    auto mixture = std::make_shared<ePCSAFTMixtureNative>(args);
    const std::string term_name = generic_term_name_cpp(record);
    try {
        if (record.term == kGenericTermDensity) {
            auto state = mixture->state(record.t, record.x, record.phase, true, record.p, false, 0.0);
            double calc = record.density_kind == 1 ? generic_mass_density_cpp(*state, args, record.x) : state->density();
            double residual = relative_residual_cpp(calc, record.target);
            raw_by_term[term_name].push_back(residual);
            scaled.push_back(record.scale * residual);
            return;
        }
        if (record.term == kGenericTermPureVLE) {
            auto liquid = mixture->state(record.t, record.x, 0, true, record.p, false, 0.0);
            auto vapor = mixture->state(record.t, record.x, 1, true, record.p, false, 0.0);
            vector<double> lnphi_liq = liquid->ln_fugacity_coefficient();
            vector<double> lnphi_vap = vapor->ln_fugacity_coefficient();
            double residual = lnphi_liq.at(0) - lnphi_vap.at(0);
            raw_by_term[term_name].push_back(residual);
            scaled.push_back(record.scale * residual);
            return;
        }
        if (record.term == kGenericTermOsmotic) {
            auto state = mixture->state(record.t, record.x, record.phase, true, record.p, false, 0.0);
            double residual = relative_or_absolute_residual_cpp(state->osmotic_coefficient(), record.target);
            raw_by_term[term_name].push_back(residual);
            scaled.push_back(record.scale * residual);
            return;
        }
        if (record.term == kGenericTermMIAC) {
            auto state = mixture->state(record.t, record.x, record.phase, true, record.p, false, 0.0);
            ActivityCoefficientNative activity = state->activity_coefficient_native(
                false,
                record.solvent_index >= 0,
                record.solvent_index
            );
            int pair_index = -1;
            for (size_t i = 0; i < activity.pair_cation_indices.size(); ++i) {
                bool cation_matches = record.target_index < 0 || activity.pair_cation_indices[i] == record.target_index;
                bool anion_matches = record.target_index_2 < 0 || activity.pair_anion_indices[i] == record.target_index_2;
                if (cation_matches && anion_matches) {
                    if (pair_index >= 0 && record.target_index < 0 && record.target_index_2 < 0) {
                        throw ValueError("Native generic regression MIAC records require a pair when multiple ion pairs are present.");
                    }
                    pair_index = static_cast<int>(i);
                }
            }
            if (pair_index < 0) {
                throw ValueError("Native generic regression MIAC pair was not present in the activity coefficient state.");
            }
            const vector<double> &values = record.activity_basis == 1
                ? activity.mean_ionic_activity_coefficients_molality
                : activity.mean_ionic_activity_coefficients_mole_fraction;
            double residual = relative_or_absolute_residual_cpp(values.at(static_cast<size_t>(pair_index)), record.target);
            raw_by_term[term_name].push_back(residual);
            scaled.push_back(record.scale * residual);
            return;
        }
        if (record.term == kGenericTermBinaryVLE) {
            auto liquid = mixture->state(record.t, record.x, 0, true, record.p, false, 0.0);
            auto vapor = mixture->state(record.t, record.y, 1, true, record.p, false, 0.0);
            vector<double> lnphi_liq = liquid->ln_fugacity_coefficient();
            vector<double> lnphi_vap = vapor->ln_fugacity_coefficient();
            for (int index : {record.target_index, record.target_index_2}) {
                if (index < 0) {
                    continue;
                }
                const double xi = record.x.at(static_cast<size_t>(index));
                const double yi = record.y.at(static_cast<size_t>(index));
                if (!(xi > 0.0) || !(yi > 0.0)) {
                    throw ValueError("Native generic regression binary VLE records require positive x/y values.");
                }
                double residual = std::log(xi) + lnphi_liq.at(static_cast<size_t>(index))
                    - std::log(yi) - lnphi_vap.at(static_cast<size_t>(index));
                raw_by_term[term_name].push_back(residual);
                scaled.push_back(record.scale * residual);
            }
            return;
        }
        if (record.term == kGenericTermComponentLnFugacity) {
            auto state = mixture->state(record.t, record.x, record.phase, true, record.p, false, 0.0);
            vector<double> lnphi = state->ln_fugacity_coefficient();
            double residual = lnphi.at(static_cast<size_t>(record.target_index)) - record.target;
            raw_by_term[term_name].push_back(residual);
            scaled.push_back(record.scale * residual);
            return;
        }
        throw ValueError("Unknown native generic regression residual term.");
    } catch (const SolutionError&) {
        int penalty_count = 1;
        if (record.term == kGenericTermBinaryVLE) {
            penalty_count = 0;
            if (record.target_index >= 0) {
                ++penalty_count;
            }
            if (record.target_index_2 >= 0) {
                ++penalty_count;
            }
            penalty_count = std::max(1, penalty_count);
        }
        for (int i = 0; i < penalty_count; ++i) {
            double residual = 1.0e6;
            raw_by_term[term_name].push_back(residual);
            scaled.push_back(record.scale * residual);
        }
    }
}

GenericRegressionDebugResult evaluate_generic_residuals_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &theta
) {
    if (base_args_by_record.size() != records.size()) {
        throw ValueError("Native generic regression requires one base parameter payload per record.");
    }
    vector<double> scaled;
    std::map<std::string, vector<double>> raw_by_term;
    for (size_t i = 0; i < records.size(); ++i) {
        append_generic_residuals_cpp(
            base_args_by_record[i],
            records[i],
            target_kinds,
            target_indices,
            target_indices_2,
            theta,
            scaled,
            raw_by_term
        );
    }
    if (scaled.empty()) {
        throw ValueError("Native generic regression generated no residuals.");
    }
    GenericRegressionDebugResult out;
    out.residuals = std::move(scaled);
    out.cost = 0.0;
    for (double value : out.residuals) {
        out.cost += 0.5 * value * value;
    }
    out.residual_norm = std::sqrt(std::max(0.0, 2.0 * out.cost));
    out.jacobian_available = true;
    out.jacobian_backend = "finite_difference";
    out.jacobian_fallback_used = true;
    out.jacobian_fallback_reason =
        "Generic regression autodiff Jacobian is not implemented for all residual state calls yet.";
    out.finite_difference_fallback_count = static_cast<int>(target_kinds.size());
    for (const auto &item : raw_by_term) {
        out.metrics_by_term[item.first] = rms_metric_cpp(item.second);
    }
    return out;
}

GenericRegressionDebugResult evaluate_generic_residuals_with_jacobian_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &theta
) {
    GenericRegressionDebugResult out = evaluate_generic_residuals_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        theta
    );
    const Index rows = static_cast<Index>(out.residuals.size());
    const Index cols = static_cast<Index>(theta.size());
    ResidualJacobian jac = ResidualJacobian::Zero(rows, cols);
    for (Index j = 0; j < cols; ++j) {
        const double eps = 1.0e-6 * std::max(1.0, std::abs(theta[static_cast<size_t>(j)]));
        vector<double> xp = theta;
        vector<double> xm = theta;
        xp[static_cast<size_t>(j)] += eps;
        xm[static_cast<size_t>(j)] -= eps;
        GenericRegressionDebugResult fp = evaluate_generic_residuals_cpp(
            base_args_by_record, records, target_kinds, target_indices, target_indices_2, xp
        );
        GenericRegressionDebugResult fm = evaluate_generic_residuals_cpp(
            base_args_by_record, records, target_kinds, target_indices, target_indices_2, xm
        );
        for (Index i = 0; i < rows; ++i) {
            jac(i, j) = (fp.residuals[static_cast<size_t>(i)] - fm.residuals[static_cast<size_t>(i)]) / (2.0 * eps);
        }
    }
    out.jacobian_rows = static_cast<int>(rows);
    out.jacobian_cols = static_cast<int>(cols);
    out.jacobian_row_major.resize(static_cast<size_t>(rows * cols), 0.0);
    for (Index i = 0; i < rows; ++i) {
        for (Index j = 0; j < cols; ++j) {
            out.jacobian_row_major[static_cast<size_t>(i * cols + j)] = jac(i, j);
        }
    }
    out.jacobian_available = true;
    out.jacobian_backend = "finite_difference";
    out.jacobian_fallback_used = true;
    out.jacobian_fallback_reason =
        "Generic regression autodiff Jacobian is not implemented for all residual state calls yet.";
    out.finite_difference_fallback_count = static_cast<int>(cols);
    return out;
}

struct GenericBoundedTransformResult {
    vector<double> x;
};

GenericBoundedTransformResult generic_unconstrained_to_bounded_cpp(
    const LMInputVector &y,
    const vector<double> &lower,
    const vector<double> &upper
) {
    GenericBoundedTransformResult out;
    out.x.resize(static_cast<size_t>(y.size()), 0.0);
    for (Index i = 0; i < y.size(); ++i) {
        double lo = lower[static_cast<size_t>(i)];
        double hi = upper[static_cast<size_t>(i)];
        if (!std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo)) {
            throw ValueError("Native generic regression requires finite strictly increasing bounds.");
        }
        double sigma = logistic_cpp(y[i]);
        out.x[static_cast<size_t>(i)] = lo + (hi - lo) * sigma;
    }
    return out;
}

LMInputVector generic_bounded_to_unconstrained_cpp(
    const vector<double> &x,
    const vector<double> &lower,
    const vector<double> &upper
) {
    LMInputVector y(static_cast<Index>(x.size()));
    for (size_t i = 0; i < x.size(); ++i) {
        double lo = lower[i];
        double hi = upper[i];
        if (!std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo)) {
            throw ValueError("Native generic regression requires finite strictly increasing bounds.");
        }
        double clipped = clip_start_value_cpp(x[i], lo, hi);
        double p = (clipped - lo) / (hi - lo);
        p = std::min(1.0 - 1.0e-12, std::max(1.0e-12, p));
        y[static_cast<Index>(i)] = logit_cpp(p);
    }
    return y;
}

vector<vector<double>> generic_candidate_starts_cpp(
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart
) {
    if (x0.size() != lower.size() || x0.size() != upper.size()) {
        throw ValueError("Native generic regression starts and bounds must have matching lengths.");
    }
    vector<vector<double>> starts;
    vector<double> first = x0;
    for (size_t i = 0; i < first.size(); ++i) {
        first[i] = clip_start_value_cpp(first[i], lower[i], upper[i]);
    }
    starts.push_back(first);
    if (multistart <= 0) {
        return starts;
    }
    constexpr std::array<double, 5> fractions = {0.25, 0.5, 0.75, 0.1, 0.9};
    for (double fraction : fractions) {
        if (static_cast<int>(starts.size()) > multistart) {
            break;
        }
        vector<double> point = first;
        for (size_t i = 0; i < point.size(); ++i) {
            point[i] = lower[i] + fraction * (upper[i] - lower[i]);
        }
        append_start_if_distinct_cpp(starts, std::move(point));
    }
    return starts;
}

struct GenericLeastSquaresFunctor : Eigen::DenseFunctor<double> {
    static int residual_count_from_records(const vector<GenericRegressionRecord> &records) {
        int count = 0;
        for (const auto &record : records) {
            if (record.term == kGenericTermBinaryVLE) {
                if (record.target_index >= 0) {
                    ++count;
                }
                if (record.target_index_2 >= 0) {
                    ++count;
                }
            } else {
                ++count;
            }
        }
        return std::max(1, count);
    }

    GenericLeastSquaresFunctor(
        vector<add_args> base_args_by_record,
        vector<GenericRegressionRecord> records,
        vector<int> target_kinds,
        vector<int> target_indices,
        vector<int> target_indices_2,
        vector<double> lower,
        vector<double> upper
    )
        : Eigen::DenseFunctor<double>(
              static_cast<int>(target_kinds.size()),
              residual_count_from_records(records)
          ),
          base_args_by_record(std::move(base_args_by_record)),
          records(std::move(records)),
          target_kinds(std::move(target_kinds)),
          target_indices(std::move(target_indices)),
          target_indices_2(std::move(target_indices_2)),
          lower(std::move(lower)),
          upper(std::move(upper)) {}

    int operator()(const LMInputVector &y, ResidualVector &fvec) {
        try {
            GenericBoundedTransformResult transform = generic_unconstrained_to_bounded_cpp(y, lower, upper);
            GenericRegressionDebugResult eval = evaluate_generic_residuals_cpp(
                base_args_by_record,
                records,
                target_kinds,
                target_indices,
                target_indices_2,
                transform.x
            );
            for (Index i = 0; i < fvec.size(); ++i) {
                fvec[i] = eval.residuals[static_cast<size_t>(i)];
            }
            return 0;
        } catch (...) {
            return -1;
        }
    }

    int df(const LMInputVector &y, ResidualJacobian &fjac) {
        constexpr double eps = 1.0e-6;
        ResidualVector f0(values());
        if ((*this)(y, f0) != 0) {
            return -1;
        }
        for (Index j = 0; j < inputs(); ++j) {
            LMInputVector yp = y;
            LMInputVector ym = y;
            yp[j] += eps;
            ym[j] -= eps;
            ResidualVector fp(values());
            ResidualVector fm(values());
            if ((*this)(yp, fp) != 0 || (*this)(ym, fm) != 0) {
                return -1;
            }
            fjac.col(j) = (fp - fm) / (2.0 * eps);
        }
        return 0;
    }

    vector<add_args> base_args_by_record;
    vector<GenericRegressionRecord> records;
    vector<int> target_kinds;
    vector<int> target_indices;
    vector<int> target_indices_2;
    vector<double> lower;
    vector<double> upper;
};

GenericRegressionResult generic_result_from_eval_cpp(
    const GenericRegressionDebugResult &eval,
    const vector<double> &x,
    bool success,
    int status,
    const std::string &message,
    int nfev,
    int iterations
) {
    GenericRegressionResult out;
    out.x = x;
    out.cost = eval.cost;
    out.residual_norm = eval.residual_norm;
    out.metrics_by_term = eval.metrics_by_term;
    out.success = success;
    out.status = status;
    out.message = message;
    out.nfev = nfev;
    out.iterations = iterations;
    out.backend = "least_squares_native";
    out.jacobian_available = eval.jacobian_available;
    out.jacobian_backend = eval.jacobian_backend;
    out.jacobian_fallback_used = eval.jacobian_fallback_used;
    out.jacobian_fallback_reason = eval.jacobian_fallback_reason;
    out.finite_difference_fallback_count = eval.finite_difference_fallback_count;
    out.hessian_available = eval.hessian_available;
    out.hessian_backend = eval.hessian_backend;
    out.hessian_fallback_used = eval.hessian_fallback_used;
    out.hessian_fallback_reason = eval.hessian_fallback_reason;
    return out;
}

GenericRegressionResult solve_one_generic_start_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &start,
    const vector<double> &lower,
    const vector<double> &upper,
    int max_nfev
) {
    GenericRegressionDebugResult initial_eval = evaluate_generic_residuals_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        start
    );
    if (max_nfev <= 1) {
        GenericRegressionResult out = generic_result_from_eval_cpp(
            initial_eval,
            start,
            std::isfinite(initial_eval.residual_norm),
            0,
            "evaluated initial native generic residual without optimizer",
            1,
            0
        );
        out.initial_cost = initial_eval.cost;
        out.initial_residual_norm = initial_eval.residual_norm;
        return out;
    }
    GenericLeastSquaresFunctor functor(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        lower,
        upper
    );
    Eigen::LevenbergMarquardt<GenericLeastSquaresFunctor> lm(functor);
    lm.setFtol(1.0e-6);
    lm.setXtol(1.0e-6);
    lm.setGtol(0.0);
    lm.setFactor(10.0);
    lm.setMaxfev(max_nfev);

    LMInputVector y = generic_bounded_to_unconstrained_cpp(start, lower, upper);
    Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(y);
    GenericBoundedTransformResult final_transform = generic_unconstrained_to_bounded_cpp(y, lower, upper);
    GenericRegressionDebugResult final_eval = evaluate_generic_residuals_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        final_transform.x
    );
    GenericRegressionResult out = generic_result_from_eval_cpp(
        final_eval,
        final_transform.x,
        least_squares_status_success_cpp(status) || final_eval.residual_norm <= initial_eval.residual_norm + 1.0e-12,
        static_cast<int>(status),
        least_squares_status_message_cpp(status),
        static_cast<int>(lm.nfev() + lm.njev()),
        static_cast<int>(lm.iterations())
    );
    out.initial_cost = initial_eval.cost;
    out.initial_residual_norm = initial_eval.residual_norm;
    return out;
}

GenericRegressionResult choose_better_generic_result_cpp(
    bool have_result,
    const GenericRegressionResult &best,
    const GenericRegressionResult &candidate
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

}  // namespace regression_detail

using regression_detail::PureNeutralObjectiveEvaluation;
using regression_detail::PureNeutralResidualEvaluation;
using regression_detail::RegressionProfilingStats;
using regression_detail::candidate_starts_cpp;
using regression_detail::choose_better_result_cpp;
using regression_detail::evaluate_pure_neutral_objective_cpp;
using regression_detail::evaluate_residual_jacobian_cpp;
using regression_detail::kThetaSize;
using regression_detail::objective_from_residual_eval_cpp;
using regression_detail::solve_one_start_least_squares_cpp;
using regression_detail::validate_pure_neutral_base_args_cpp;
using regression_detail::choose_better_generic_result_cpp;
using regression_detail::evaluate_generic_residuals_cpp;
using regression_detail::evaluate_generic_residuals_with_jacobian_cpp;
using regression_detail::generic_candidate_starts_cpp;
using regression_detail::solve_one_generic_start_cpp;

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
    int starts_tried = 0;
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
            ++starts_tried;
            best = choose_better_result_cpp(have_result, best, candidate);
            have_result = true;
        } catch (...) {
        }
    }
    if (!have_result) {
        throw ValueError("Native pure-neutral least-squares regression did not generate any candidate starts.");
    }
    best.starts_tried = starts_tried;
    return best;
}

GenericRegressionDebugResult evaluate_generic_regression_debug_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &x
) {
    return evaluate_generic_residuals_with_jacobian_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        x
    );
}

GenericRegressionResult fit_generic_least_squares_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int multistart,
    int max_nfev
) {
    if (target_kinds.empty()) {
        throw ValueError("Native generic regression requires at least one optimization target.");
    }
    if (x0.size() != target_kinds.size() || lower.size() != target_kinds.size() || upper.size() != target_kinds.size()) {
        throw ValueError("Native generic regression target arrays must have matching lengths.");
    }
    vector<vector<double>> starts = generic_candidate_starts_cpp(x0, lower, upper, multistart);
    bool have_result = false;
    GenericRegressionResult best;
    int starts_tried = 0;
    int nfev_total = 0;
    for (const auto &start : starts) {
        GenericRegressionResult candidate = solve_one_generic_start_cpp(
            base_args_by_record,
            records,
            target_kinds,
            target_indices,
            target_indices_2,
            start,
            lower,
            upper,
            max_nfev
        );
        ++starts_tried;
        nfev_total += candidate.nfev;
        best = choose_better_generic_result_cpp(have_result, best, candidate);
        have_result = true;
    }
    if (!have_result) {
        throw ValueError("Native generic least-squares regression did not generate any candidate starts.");
    }
    best.starts_tried = starts_tried;
    best.nfev = nfev_total;
    best.backend = "least_squares_native";
    return best;
}

