#include "epcsaft_core_internal.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>


#ifndef EPCSAFT_HAS_CERES
#error "EPCSAFT_HAS_CERES must be defined; Ceres is a required native dependency."
#endif

#include <ceres/cost_function.h>
#include <ceres/problem.h>
#include <ceres/solver.h>

using Index = Eigen::Index;
using thermo_detail::kDispersionA0;
using thermo_detail::kDispersionA1;
using thermo_detail::kDispersionA2;
using thermo_detail::kDispersionB0;
using thermo_detail::kDispersionB1;
using thermo_detail::kDispersionB2;

namespace regression_detail {

constexpr double kRegressionGradientFloor = 1.0e-300;
constexpr int kThetaSize = 3;
constexpr int kSeedRho = 0;
constexpr int kSeedM = 1;
constexpr int kSeedS = 2;
constexpr int kSeedE = 3;
constexpr int kSeedCount = 4;

using ResidualVector = Eigen::VectorXd;
using ResidualJacobian = Eigen::MatrixXd;

double scalar_value(double x) {
    return x;
}

#ifdef EPCSAFT_HAS_CPPAD
double scalar_value(const CppADScalar &x) {
    return CppAD::Value(x);
}
#endif

bool scalar_domain_positive(double x) {
    return x > 0.0;
}

#ifdef EPCSAFT_HAS_CPPAD
bool scalar_domain_positive(const CppADScalar &) {
    return true;
}
#endif

double scalar_log(double x) {
    return std::log(x);
}

#ifdef EPCSAFT_HAS_CPPAD
CppADScalar scalar_log(const CppADScalar &x) {
    return CppAD::log(x);
}
#endif

double scalar_exp(double x) {
    return std::exp(x);
}

#ifdef EPCSAFT_HAS_CPPAD
CppADScalar scalar_exp(const CppADScalar &x) {
    return CppAD::exp(x);
}
#endif

double scalar_pow(double x, int exponent) {
    return std::pow(x, exponent);
}

#ifdef EPCSAFT_HAS_CPPAD
CppADScalar scalar_pow(const CppADScalar &x, int exponent) {
    return CppAD::pow(x, static_cast<double>(exponent));
}
#endif

double scalar_pow(double x, double exponent) {
    return std::pow(x, exponent);
}

#ifdef EPCSAFT_HAS_CPPAD
CppADScalar scalar_pow(const CppADScalar &x, double exponent) {
    return CppAD::pow(x, exponent);
}
#endif

template <typename Scalar>
Scalar regression_scalar_constant(double value) {
    return scalar_constant<Scalar>(value);
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
    Scalar mures = regression_scalar_constant<Scalar>(0.0);
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

double pure_neutral_gradient_norm_cpp(const PureNeutralObjectiveEvaluation &eval) {
    double norm_sq = 0.0;
    for (double value : eval.gradient) {
        norm_sq += value * value;
    }
    return std::sqrt(norm_sq);
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
    if (!scalar_domain_positive(state.Z)) {
        throw ValueError("Encountered non-positive compressibility factor during native regression evaluation.");
    }
    state.pressure = state.Z * kb * t * state.den * 1.0e30;
    state.mures = state.ares_total + state.zraw_total;
    state.lnfug = state.ares_total + state.zraw_total - scalar_log(state.Z);
    return state;
}

PureNeutralFusedState evaluate_fused_state_cpp(double t, double rho, const vector<double> &x) {
    std::vector<CppADScalar> ax(kSeedCount);
    ax[kSeedRho] = rho;
    ax[kSeedM] = x[0];
    ax[kSeedS] = x[1];
    ax[kSeedE] = x[2];
    CppAD::Independent(ax);

    auto state = pure_neutral_state_scalar_cpp<CppADScalar>(
        t,
        ax[kSeedRho],
        ax[kSeedM],
        ax[kSeedS],
        ax[kSeedE]
    );
    std::vector<CppADScalar> ay(3);
    ay[0] = state.pressure;
    ay[1] = state.lnfug;
    ay[2] = state.Z;

    CppAD::ADFun<double> function(ax, ay);
    std::vector<double> point = {rho, x[0], x[1], x[2]};
    std::vector<double> value = function.Forward(0, point);
    std::vector<double> jacobian = function.Jacobian(point);

    PureNeutralFusedState out;
    out.pressure = value[0];
    out.lnfug = value[1];
    out.Z = value[2];
    out.dpdrho = jacobian[0 * kSeedCount + kSeedRho];
    out.dlnfug_drho = jacobian[1 * kSeedCount + kSeedRho];
    for (int j = 0; j < kThetaSize; ++j) {
        out.dpdtheta[static_cast<size_t>(j)] = jacobian[0 * kSeedCount + j + 1];
        out.dlnfugdtheta[static_cast<size_t>(j)] = jacobian[1 * kSeedCount + j + 1];
    }
    return out;
}

epcsaft::native::cppad_support::CppADDerivativeResult cppad_pure_neutral_parameter_derivatives_cpp(
    double t,
    double rho,
    const add_args &base_args
) {
    if (base_args.m.size() != 1 || base_args.s.size() != 1 || base_args.e.size() != 1) {
        throw ValueError("unsupported: pure-neutral m/sigma/epsilon derivatives require exactly one component.");
    }
    if (!base_args.z.empty() && (base_args.z.size() != 1 || std::abs(base_args.z[0]) > 1.0e-12)) {
        throw ValueError("unsupported: pure-neutral m/sigma/epsilon derivatives support only neutral components.");
    }
    if (!base_args.assoc_num.empty() || !base_args.assoc_matrix.empty() || !base_args.k_hb.empty() || !base_args.e_assoc.empty() || !base_args.vol_a.empty()) {
        throw ValueError("unsupported: pure-neutral m/sigma/epsilon derivatives support only nonassociating components.");
    }
    using CppADScalar = CppAD::AD<double>;
    std::vector<CppADScalar> ax(kThetaSize);
    ax[0] = base_args.m[0];
    ax[1] = base_args.s[0];
    ax[2] = base_args.e[0];
    CppAD::Independent(ax);

    auto state = pure_neutral_state_scalar_cpp<CppADScalar>(
        t,
        CppADScalar(rho),
        ax[0],
        ax[1],
        ax[2]
    );
    std::vector<CppADScalar> ay(3);
    ay[0] = state.pressure;
    ay[1] = state.mures;
    ay[2] = state.lnfug;

    CppAD::ADFun<double> function(ax, ay);
    std::vector<double> point = {base_args.m[0], base_args.s[0], base_args.e[0]};
    auto value = function.Forward(0, point);
    auto jacobian = function.Jacobian(point);

    epcsaft::native::cppad_support::CppADDerivativeResult result;
    result.supported = true;
    result.backend = "cppad";
    result.message = "CppAD pure-neutral m/sigma/epsilon property derivatives available";
    result.value = std::move(value);
    result.jacobian_row_major = std::move(jacobian);
    result.outputs = {"pressure", "mu_res", "ln_fugacity"};
    result.variables = {"m", "sigma", "epsilon"};
    result.rows = 3;
    result.cols = kThetaSize;
    return result;
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

vector<double> canonical_start_cpp(
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper
) {
    if (x0.size() != lower.size() || x0.size() != upper.size()) {
        throw ValueError("Native regression starts and bounds must have matching lengths.");
    }
    vector<double> start = x0;
    for (size_t i = 0; i < start.size(); ++i) {
        start[i] = clip_start_value_cpp(start[i], lower[i], upper[i]);
    }
    return start;
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

#ifdef EPCSAFT_HAS_CERES
// AlgID: native_ceres_regression_adapter
ceres::Solver::Options ceres_regression_options_cpp(int max_num_iterations) {
    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = false;
    options.logging_type = ceres::SILENT;
    options.max_num_iterations = std::max(1, max_num_iterations);
    options.function_tolerance = 1.0e-10;
    options.gradient_tolerance = 1.0e-10;
    options.parameter_tolerance = 1.0e-10;
    return options;
}

// AlgID: pure_neutral_ceres_regression
class PureNeutralCeresCostFunction final : public ceres::CostFunction {
public:
    PureNeutralCeresCostFunction(
        add_args base_args,
        vector<PureNeutralRegressionDensityRecord> density_records,
        double density_scale,
        vector<PureNeutralRegressionVLERecord> pure_vle_records,
        double pure_vle_scale,
        RegressionProfilingStats *profiling
    )
        : base_args_(std::move(base_args)),
          density_records_(std::move(density_records)),
          density_scale_(density_scale),
          pure_vle_records_(std::move(pure_vle_records)),
          pure_vle_scale_(pure_vle_scale),
          profiling_(profiling)
    {
        set_num_residuals(static_cast<int>(density_records_.size() + pure_vle_records_.size()));
        mutable_parameter_block_sizes()->push_back(kThetaSize);
    }

    bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override {
        try {
            vector<double> theta(parameters[0], parameters[0] + kThetaSize);
            PureNeutralResidualEvaluation eval = evaluate_residual_jacobian_cpp(
                base_args_,
                density_records_,
                density_scale_,
                pure_vle_records_,
                pure_vle_scale_,
                theta,
                profiling_
            );
            for (Eigen::Index row = 0; row < eval.residuals.size(); ++row) {
                residuals[row] = eval.residuals[row];
            }
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                for (Eigen::Index row = 0; row < eval.jacobian.rows(); ++row) {
                    for (Eigen::Index col = 0; col < eval.jacobian.cols(); ++col) {
                        jacobians[0][row * kThetaSize + col] = eval.jacobian(row, col);
                    }
                }
            }
            return true;
        } catch (...) {
            return false;
        }
    }

private:
    add_args base_args_;
    vector<PureNeutralRegressionDensityRecord> density_records_;
    double density_scale_ = 1.0;
    vector<PureNeutralRegressionVLERecord> pure_vle_records_;
    double pure_vle_scale_ = 1.0;
    RegressionProfilingStats *profiling_ = nullptr;
};

PureNeutralRegressionResult solve_one_start_ceres_cpp(
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

    vector<double> theta = start;
    RegressionProfilingStats profiling;
    ceres::Problem problem;
    auto *cost = new PureNeutralCeresCostFunction(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        &profiling
    );
    problem.AddResidualBlock(cost, nullptr, theta.data());
    for (int j = 0; j < kThetaSize; ++j) {
        problem.SetParameterLowerBound(theta.data(), j, lower[static_cast<size_t>(j)]);
        problem.SetParameterUpperBound(theta.data(), j, upper[static_cast<size_t>(j)]);
    }

    ceres::Solver::Options options = ceres_regression_options_cpp(100);

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    PureNeutralResidualEvaluation final_eval = evaluate_residual_jacobian_cpp(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        theta,
        &profiling
    );
    PureNeutralObjectiveEvaluation final_objective = objective_from_residual_eval_cpp(final_eval);
    PureNeutralRegressionResult out;
    out.x = theta;
    out.initial_cost = initial_eval.objective;
    out.initial_density_metric = rms_metric_cpp(initial_eval.density_raw_residuals);
    out.initial_pure_vle_metric = rms_metric_cpp(initial_eval.pure_vle_raw_residuals);
    out.cost = final_objective.objective;
    out.residual_norm = std::sqrt(std::max(0.0, 2.0 * final_objective.objective));
    out.density_metric = rms_metric_cpp(final_objective.density_raw_residuals);
    out.pure_vle_metric = rms_metric_cpp(final_objective.pure_vle_raw_residuals);
    out.success = summary.IsSolutionUsable();
    out.status = static_cast<int>(summary.termination_type);
    out.message = summary.BriefReport();
    out.nfev = static_cast<int>(summary.num_residual_evaluations + summary.num_jacobian_evaluations);
    out.iterations = static_cast<int>(summary.iterations.size());
    out.objective_evaluations = static_cast<int>(summary.num_residual_evaluations);
    out.gradient_evaluations = static_cast<int>(summary.num_jacobian_evaluations);
    out.residual_evaluations = profiling.residual_evaluations;
    out.density_solves = profiling.density_solves;
    out.fused_state_evaluations = profiling.fused_state_evaluations;
    out.callback_wall_time_s = profiling.callback_wall_time_s;
    out.solve_wall_time_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - solve_start).count();
    out.backend = "ceres";
    out.optimizer_backend = "ceres";
    out.derivative_backend = "cppad_implicit";
    out.jacobian_backend = out.derivative_backend;
    out.gradient_norm = pure_neutral_gradient_norm_cpp(final_objective);
    double step_norm_sq = 0.0;
    for (int j = 0; j < kThetaSize; ++j) {
        double diff = theta[static_cast<size_t>(j)] - start[static_cast<size_t>(j)];
        step_norm_sq += diff * diff;
    }
    out.step_norm = std::sqrt(step_norm_sq);
    return out;
}
#endif

constexpr int kGenericTargetM = 0;
constexpr int kGenericTargetS = 1;
constexpr int kGenericTargetE = 2;
constexpr int kGenericTargetEAssoc = 3;
constexpr int kGenericTargetVolA = 4;
constexpr int kGenericTargetDBorn = 5;
constexpr int kGenericTargetKIJ = 6;
constexpr int kGenericTargetLIJ = 7;
constexpr int kGenericTargetKHB = 8;
constexpr int kGenericTargetFSolv = 9;
constexpr int kGenericTargetDielc = 10;

constexpr int kGenericTermDensity = 1;
constexpr int kGenericTermPureVLE = 2;
constexpr int kGenericTermOsmotic = 3;
constexpr int kGenericTermMIAC = 4;
constexpr int kGenericTermBinaryVLE = 5;
constexpr int kGenericTermComponentLnFugacity = 6;
constexpr int kGenericTermRelativePermittivity = 7;

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
        case kGenericTermRelativePermittivity:
            return "relative_permittivity";
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
            case kGenericTargetFSolv:
                set_vector_value_cpp(args.f_solv, index, value, "f_solv");
                break;
            case kGenericTargetDielc:
                set_vector_value_cpp(args.dielc, index, value, "dielc");
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
        if (record.term == kGenericTermRelativePermittivity) {
            auto state = mixture->state(record.t, record.x, record.phase, true, record.p, false, 0.0);
            vector<double> epsilon = state->relative_permittivity();
            double residual = relative_or_absolute_residual_cpp(epsilon.at(0), record.target);
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
    out.jacobian_available = false;
    for (const auto &item : raw_by_term) {
        out.metrics_by_term[item.first] = rms_metric_cpp(item.second);
    }
    return out;
}

struct BinaryKijResidualEvaluation {
    ResidualVector residuals;
    ResidualJacobian jacobian;
    double cost = 0.0;
    double residual_norm = 0.0;
    std::map<std::string, double> metrics_by_term;
};

struct PureIonBornResidualEvaluation {
    ResidualVector residuals;
    ResidualJacobian jacobian;
    double cost = 0.0;
    double residual_norm = 0.0;
    std::map<std::string, double> metrics_by_term;
};

int generic_residual_count_from_records_cpp(const vector<GenericRegressionRecord> &records) {
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

double regression_normalize_mw_cpp(double mw) {
    if (mw > 1.0) {
        return mw / 1000.0;
    }
    return mw;
}

void validate_pure_ion_born_ceres_problem_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &theta
) {
    if (base_args_by_record.size() != records.size()) {
        throw ValueError("Native Ceres pure-ion regression requires one base parameter payload per record.");
    }
    if (target_kinds.empty()) {
        throw ValueError("Native Ceres pure-ion regression requires at least one d_born target.");
    }
    if (target_kinds.size() != theta.size() || target_indices.size() != theta.size() || target_indices_2.size() != theta.size()) {
        throw ValueError("Native Ceres pure-ion regression target arrays must match theta length.");
    }
    for (size_t j = 0; j < target_kinds.size(); ++j) {
        if (target_kinds[j] != kGenericTargetS
            && target_kinds[j] != kGenericTargetE
            && target_kinds[j] != kGenericTargetDBorn
            && target_kinds[j] != kGenericTargetFSolv
            && target_kinds[j] != kGenericTargetDielc) {
            throw ValueError("unsupported: native Ceres pure-ion regression supports s, e, d_born, f_solv, and dielc targets only.");
        }
        (void)target_indices_2[j];
    }
    if (records.empty()) {
        throw ValueError("Native Ceres pure-ion regression requires osmotic or mean-ionic activity records.");
    }
    for (size_t r = 0; r < records.size(); ++r) {
        const auto &record = records[r];
        const auto &args = base_args_by_record[r];
        if (record.term != kGenericTermDensity
            && record.term != kGenericTermOsmotic
            && record.term != kGenericTermMIAC
            && record.term != kGenericTermRelativePermittivity) {
            throw ValueError("unsupported: native Ceres pure-ion regression supports density, osmotic, mean-ionic activity, and relative-permittivity rows only.");
        }
        if (!(record.t > 0.0) || !(record.p > 0.0) || record.phase != 0) {
            throw ValueError("Native Ceres pure-ion regression requires positive T/P liquid records.");
        }
        if (args.z.empty() || args.d_born.size() != args.z.size() || args.m.size() != args.z.size() || args.mw.size() != args.z.size()) {
            throw ValueError("Native Ceres pure-ion regression requires aligned ionic parameter payloads.");
        }
        bool has_ion = false;
        for (double charge : args.z) {
            has_ion = has_ion || std::abs(charge) > 1.0e-12;
        }
        if (!has_ion) {
            throw ValueError("Native Ceres pure-ion regression requires ionic species.");
        }
        for (size_t j = 0; j < target_indices.size(); ++j) {
            const int index = target_indices[j];
            if (index < 0 || static_cast<size_t>(index) >= args.d_born.size()) {
                throw ValueError("Native Ceres pure-ion target index is out of range.");
            }
            if (target_kinds[j] == kGenericTargetFSolv
                && (args.f_solv.size() != args.z.size() || static_cast<size_t>(index) >= args.f_solv.size())) {
                throw ValueError("Native Ceres pure-ion f_solv targets require aligned f_solv values.");
            }
            if (target_kinds[j] == kGenericTargetDielc
                && (args.dielc.size() != args.z.size() || static_cast<size_t>(index) >= args.dielc.size())) {
                throw ValueError("Native Ceres pure-ion dielc targets require aligned relative-permittivity values.");
            }
        }
    }
}

PureIonBornResidualEvaluation evaluate_pure_ion_born_residual_jacobian_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &theta
) {
    validate_pure_ion_born_ceres_problem_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        theta
    );
    PureIonBornResidualEvaluation out;
    out.residuals = ResidualVector::Zero(static_cast<int>(records.size()));
    out.jacobian = ResidualJacobian::Zero(static_cast<int>(records.size()), static_cast<int>(theta.size()));
    std::map<std::string, vector<double>> raw_by_term;
    for (size_t r = 0; r < records.size(); ++r) {
        add_args args = base_args_by_record[r];
        apply_generic_targets_cpp(args, target_kinds, target_indices, target_indices_2, theta);
        const auto &record = records[r];
        auto mixture = std::make_shared<ePCSAFTMixtureNative>(args);
        auto state = mixture->state(record.t, record.x, record.phase, true, record.p, false, 0.0);
        const double rho = state->density();
        double calc = 0.0;
        vector<double> dcalc(theta.size(), 0.0);
        if (record.term == kGenericTermDensity) {
            calc = record.density_kind == 1 ? generic_mass_density_cpp(*state, args, record.x) : state->density();
            double density_multiplier = 1.0;
            if (record.density_kind == 1) {
                if (args.mw.size() != record.x.size()) {
                    throw ValueError("Native Ceres pure-ion density derivative requires one MW value per component.");
                }
                density_multiplier = 0.0;
                for (size_t i = 0; i < record.x.size(); ++i) {
                    density_multiplier += record.x[i] * args.mw[i];
                }
            }
            for (size_t j = 0; j < theta.size(); ++j) {
                if ((target_kinds[j] == kGenericTargetDBorn || target_kinds[j] == kGenericTargetFSolv)
                    && args.born_model == 2) {
                    dcalc[j] = 0.0;
                    continue;
                }
                NeutralBinaryKijPhaseDerivatives phase = generic_component_parameter_phase_derivatives_cpp(
                    record.t,
                    rho,
                    record.x,
                    args,
                    target_kinds[j],
                    target_indices[j]
                );
                dcalc[j] = density_multiplier * phase.drhodk;
            }
        } else if (record.term == kGenericTermOsmotic) {
            ActivityCoefficientNative activity = state->activity_coefficient_native(true, record.solvent_index >= 0, record.solvent_index);
            const int solvent = activity.solvent_index;
            calc = activity.osmotic_coefficient;
            const double mw_solvent = regression_normalize_mw_cpp(args.mw.at(static_cast<size_t>(solvent)));
            double molality_sum = 0.0;
            for (size_t i = 0; i < record.x.size(); ++i) {
                if (static_cast<int>(i) == solvent) {
                    continue;
                }
                molality_sum += record.x[i] / (record.x[static_cast<size_t>(solvent)] * mw_solvent);
            }
            if (!(molality_sum > 0.0)) {
                throw ValueError("Native Ceres pure-ion osmotic record has zero total molality.");
            }
            vector<double> x0(record.x.size(), 0.0);
            x0[static_cast<size_t>(solvent)] = 1.0;
            const double rho_ref = mixture->solve_density(record.t, record.p, x0, record.phase);
            for (size_t j = 0; j < theta.size(); ++j) {
                const int target = target_indices[j];
                if ((target_kinds[j] == kGenericTargetDBorn || target_kinds[j] == kGenericTargetFSolv)
                    && args.born_model == 2) {
                    BornSSMDSDerivativeResult current = born_ssmds_liquid_derivatives_cpp(
                        record.t,
                        rho,
                        record.phase,
                        record.x,
                        args
                    );
                    if (!current.supported) {
                        throw ValueError(current.message);
                    }
                    const vector<double> &current_row_major = target_kinds[j] == kGenericTargetDBorn
                        ? current.lnfug_d_d_born_row_major
                        : current.lnfug_d_f_solv_row_major;
                    double reference_derivative = 0.0;
                    if (target >= 0 && static_cast<size_t>(target) < x0.size() && x0[static_cast<size_t>(target)] > 0.0) {
                        BornSSMDSDerivativeResult reference = born_ssmds_liquid_derivatives_cpp(
                            record.t,
                            rho_ref,
                            record.phase,
                            x0,
                            args
                        );
                        if (!reference.supported) {
                            throw ValueError(reference.message);
                        }
                        const vector<double> &reference_row_major = target_kinds[j] == kGenericTargetDBorn
                            ? reference.lnfug_d_d_born_row_major
                            : reference.lnfug_d_f_solv_row_major;
                        reference_derivative = reference_row_major.at(
                            static_cast<size_t>(target) * record.x.size() + static_cast<size_t>(solvent)
                        );
                    }
                    const double current_derivative = current_row_major.at(
                        static_cast<size_t>(target) * record.x.size() + static_cast<size_t>(solvent)
                    );
                    dcalc[j] = -(current_derivative - reference_derivative) / (mw_solvent * molality_sum);
                    continue;
                }
                NeutralBinaryKijPhaseDerivatives current = generic_component_parameter_phase_derivatives_cpp(
                    record.t,
                    rho,
                    record.x,
                    args,
                    target_kinds[j],
                    target
                );
                double reference_derivative = 0.0;
                if (target >= 0 && static_cast<size_t>(target) < x0.size() && x0[static_cast<size_t>(target)] > 0.0) {
                    NeutralBinaryKijPhaseDerivatives reference = generic_component_parameter_phase_derivatives_cpp(
                        record.t,
                        rho_ref,
                        x0,
                        args,
                        target_kinds[j],
                        target
                    );
                    reference_derivative = reference.dlnphi_dk_total.at(static_cast<size_t>(solvent));
                }
                const double dln_gamma = current.dlnphi_dk_total.at(static_cast<size_t>(solvent)) - reference_derivative;
                dcalc[j] = -dln_gamma / (mw_solvent * molality_sum);
            }
        } else if (record.term == kGenericTermMIAC) {
            ActivityCoefficientNative activity = state->activity_coefficient_native(false, record.solvent_index >= 0, record.solvent_index);
            int pair_index = -1;
            for (size_t k = 0; k < activity.pair_cation_indices.size(); ++k) {
                const bool cation_matches = record.target_index < 0 || activity.pair_cation_indices[k] == record.target_index;
                const bool anion_matches = record.target_index_2 < 0 || activity.pair_anion_indices[k] == record.target_index_2;
                if (cation_matches && anion_matches) {
                    pair_index = static_cast<int>(k);
                    break;
                }
            }
            if (pair_index < 0) {
                throw ValueError("Native Ceres pure-ion MIAC pair was not present in the activity coefficient state.");
            }
            calc = record.activity_basis == 1
                ? activity.mean_ionic_activity_coefficients_molality.at(static_cast<size_t>(pair_index))
                : activity.mean_ionic_activity_coefficients_mole_fraction.at(static_cast<size_t>(pair_index));
            const int ic = activity.pair_cation_indices.at(static_cast<size_t>(pair_index));
            const int ia = activity.pair_anion_indices.at(static_cast<size_t>(pair_index));
            const double nu_cat = static_cast<double>(activity.pair_nu_cation.at(static_cast<size_t>(pair_index)));
            const double nu_an = static_cast<double>(activity.pair_nu_anion.at(static_cast<size_t>(pair_index)));
            const double sum_nu = nu_cat + nu_an;
            for (size_t j = 0; j < theta.size(); ++j) {
                const int target = target_indices[j];
                if ((target_kinds[j] == kGenericTargetDBorn || target_kinds[j] == kGenericTargetFSolv)
                    && args.born_model == 2) {
                    BornSSMDSDerivativeResult current = born_ssmds_liquid_derivatives_cpp(
                        record.t,
                        rho,
                        record.phase,
                        record.x,
                        args
                    );
                    if (!current.supported) {
                        throw ValueError(current.message);
                    }
                    const vector<double> &current_row_major = target_kinds[j] == kGenericTargetDBorn
                        ? current.lnfug_d_d_born_row_major
                        : current.lnfug_d_f_solv_row_major;
                    const double dln_gamma_pm = (
                        nu_cat * current_row_major.at(
                            static_cast<size_t>(target) * record.x.size() + static_cast<size_t>(ic)
                        )
                        + nu_an * current_row_major.at(
                            static_cast<size_t>(target) * record.x.size() + static_cast<size_t>(ia)
                        )
                    ) / sum_nu;
                    dcalc[j] = calc * dln_gamma_pm;
                    continue;
                }
                NeutralBinaryKijPhaseDerivatives phase = generic_component_parameter_phase_derivatives_cpp(
                    record.t,
                    rho,
                    record.x,
                    args,
                    target_kinds[j],
                    target
                );
                const double dln_gamma_pm = (
                    nu_cat * phase.dlnphi_dk_total.at(static_cast<size_t>(ic))
                    + nu_an * phase.dlnphi_dk_total.at(static_cast<size_t>(ia))
                ) / sum_nu;
                dcalc[j] = calc * dln_gamma_pm;
            }
        } else if (record.term == kGenericTermRelativePermittivity) {
            vector<double> epsilon = state->relative_permittivity();
            calc = epsilon.at(0);
            for (size_t j = 0; j < theta.size(); ++j) {
                if (target_kinds[j] == kGenericTargetDielc && args.dielc_rule == 1) {
                    const int target = target_indices[j];
                    dcalc[j] = target >= 0 && static_cast<size_t>(target) < record.x.size()
                        ? record.x[static_cast<size_t>(target)]
                        : 0.0;
                } else {
                    dcalc[j] = 0.0;
                }
            }
        }
        const double residual = relative_or_absolute_residual_cpp(calc, record.target);
        raw_by_term[generic_term_name_cpp(record)].push_back(residual);
        out.residuals[static_cast<int>(r)] = record.scale * residual;
        const double residual_scale = std::abs(record.target) <= 1.0e-8 ? 1.0 : 1.0 / std::abs(record.target);
        for (size_t j = 0; j < theta.size(); ++j) {
            out.jacobian(static_cast<int>(r), static_cast<int>(j)) = record.scale * dcalc[j] * residual_scale;
        }
    }
    for (Eigen::Index idx = 0; idx < out.residuals.size(); ++idx) {
        out.cost += 0.5 * out.residuals[idx] * out.residuals[idx];
    }
    out.residual_norm = std::sqrt(std::max(0.0, 2.0 * out.cost));
    for (const auto &item : raw_by_term) {
        out.metrics_by_term[item.first] = rms_metric_cpp(item.second);
    }
    return out;
}

void validate_binary_kij_ceres_problem_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &theta
) {
    if (base_args_by_record.size() != records.size()) {
        throw ValueError("Native Ceres binary k_ij regression requires one base parameter payload per record.");
    }
    if (target_kinds.size() != 1 || target_indices.size() != 1 || target_indices_2.size() != 1 || theta.size() != 1) {
        throw ValueError("unsupported: native Ceres generic regression currently supports one binary k_ij target only.");
    }
    if (target_kinds[0] != kGenericTargetKIJ) {
        throw ValueError("unsupported: native Ceres generic regression currently supports binary k_ij targets only.");
    }
    if (records.empty()) {
        throw ValueError("Native Ceres binary k_ij regression requires at least one VLE record.");
    }
    for (size_t r = 0; r < records.size(); ++r) {
        const auto &record = records[r];
        const auto &args = base_args_by_record[r];
        if (record.term != kGenericTermBinaryVLE) {
            throw ValueError("unsupported: native Ceres binary k_ij regression supports binary VLE rows only.");
        }
        if (args.m.size() != 2 || args.s.size() != 2 || args.e.size() != 2) {
            throw ValueError("unsupported: native Ceres binary k_ij regression requires exactly two neutral components.");
        }
        if (!args.z.empty()) {
            for (double charge : args.z) {
                if (std::abs(charge) > 1.0e-12) {
                    throw ValueError("unsupported: native Ceres binary k_ij regression does not support ionic rows.");
                }
            }
        }
        if (args.k_ij.size() != 4) {
            throw ValueError("unsupported: native Ceres binary k_ij regression requires a dense 2x2 k_ij matrix.");
        }
        if (record.x.size() != 2 || record.y.size() != 2 || !(record.p > 0.0) || !(record.t > 0.0)) {
            throw ValueError("Native Ceres binary k_ij regression requires positive T/P and binary x/y records.");
        }
        for (int index : {record.target_index, record.target_index_2}) {
            if (index >= 0 && index >= 2) {
                throw ValueError("Native Ceres binary k_ij residual component index is out of range.");
            }
        }
    }
}

BinaryKijResidualEvaluation evaluate_binary_kij_residual_jacobian_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &theta
) {
    validate_binary_kij_ceres_problem_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        theta
    );
    const int nres = generic_residual_count_from_records_cpp(records);
    BinaryKijResidualEvaluation out;
    out.residuals = ResidualVector::Zero(nres);
    out.jacobian = ResidualJacobian::Zero(nres, 1);
    vector<double> raw;
    int row = 0;
    for (size_t r = 0; r < records.size(); ++r) {
        add_args args = base_args_by_record[r];
        apply_generic_targets_cpp(args, target_kinds, target_indices, target_indices_2, theta);
        const auto &record = records[r];
        const int i = target_indices[0];
        const int j = target_indices_2[0];
        const int k_index = i * 2 + j;
        const int reverse_k_index = j * 2 + i;
        const double rho_liq = den_cpp(record.t, record.p, record.x, 0, args);
        const double rho_vap = den_cpp(record.t, record.p, record.y, 1, args);
        NeutralBinaryKijPhaseDerivatives liq_forward = neutral_binary_kij_phase_derivatives_cpp(
            record.t,
            rho_liq,
            record.x,
            args,
            k_index
        );
        NeutralBinaryKijPhaseDerivatives vap_forward = neutral_binary_kij_phase_derivatives_cpp(
            record.t,
            rho_vap,
            record.y,
            args,
            k_index
        );
        NeutralBinaryKijPhaseDerivatives liq_reverse = neutral_binary_kij_phase_derivatives_cpp(
            record.t,
            rho_liq,
            record.x,
            args,
            reverse_k_index
        );
        NeutralBinaryKijPhaseDerivatives vap_reverse = neutral_binary_kij_phase_derivatives_cpp(
            record.t,
            rho_vap,
            record.y,
            args,
            reverse_k_index
        );
        for (int index : {record.target_index, record.target_index_2}) {
            if (index < 0) {
                continue;
            }
            const double xi = record.x.at(static_cast<size_t>(index));
            const double yi = record.y.at(static_cast<size_t>(index));
            if (!(xi > 0.0) || !(yi > 0.0)) {
                throw ValueError("Native Ceres binary k_ij VLE records require positive x/y values.");
            }
            const double residual = std::log(xi) + liq_forward.lnphi.at(static_cast<size_t>(index))
                - std::log(yi) - vap_forward.lnphi.at(static_cast<size_t>(index));
            const double derivative =
                liq_forward.dlnphi_dk_total.at(static_cast<size_t>(index))
                + liq_reverse.dlnphi_dk_total.at(static_cast<size_t>(index))
                - vap_forward.dlnphi_dk_total.at(static_cast<size_t>(index))
                - vap_reverse.dlnphi_dk_total.at(static_cast<size_t>(index));
            raw.push_back(residual);
            out.residuals[row] = record.scale * residual;
            out.jacobian(row, 0) = record.scale * derivative;
            ++row;
        }
    }
    out.cost = 0.0;
    for (Eigen::Index idx = 0; idx < out.residuals.size(); ++idx) {
        out.cost += 0.5 * out.residuals[idx] * out.residuals[idx];
    }
    out.residual_norm = std::sqrt(std::max(0.0, 2.0 * out.cost));
    out.metrics_by_term["binary_vle_fugacity_balance"] = rms_metric_cpp(raw);
    return out;
}

#ifdef EPCSAFT_HAS_CERES
// AlgID: pure_ion_ceres_regression
class PureIonCeresCostFunction final : public ceres::CostFunction {
public:
    PureIonCeresCostFunction(
        vector<add_args> base_args_by_record,
        vector<GenericRegressionRecord> records,
        vector<int> target_kinds,
        vector<int> target_indices,
        vector<int> target_indices_2
    )
        : base_args_by_record_(std::move(base_args_by_record)),
          records_(std::move(records)),
          target_kinds_(std::move(target_kinds)),
          target_indices_(std::move(target_indices)),
          target_indices_2_(std::move(target_indices_2))
    {
        set_num_residuals(static_cast<int>(records_.size()));
        mutable_parameter_block_sizes()->push_back(static_cast<int>(target_kinds_.size()));
    }

    bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override {
        try {
            vector<double> theta(parameters[0], parameters[0] + target_kinds_.size());
            PureIonBornResidualEvaluation eval = evaluate_pure_ion_born_residual_jacobian_cpp(
                base_args_by_record_,
                records_,
                target_kinds_,
                target_indices_,
                target_indices_2_,
                theta
            );
            for (Eigen::Index row = 0; row < eval.residuals.size(); ++row) {
                residuals[row] = eval.residuals[row];
            }
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                for (Eigen::Index row = 0; row < eval.jacobian.rows(); ++row) {
                    for (Eigen::Index col = 0; col < eval.jacobian.cols(); ++col) {
                        jacobians[0][row * eval.jacobian.cols() + col] = eval.jacobian(row, col);
                    }
                }
            }
            return true;
        } catch (...) {
            return false;
        }
    }

private:
    vector<add_args> base_args_by_record_;
    vector<GenericRegressionRecord> records_;
    vector<int> target_kinds_;
    vector<int> target_indices_;
    vector<int> target_indices_2_;
};

GenericRegressionResult solve_one_pure_ion_ceres_start_cpp(
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
    PureIonBornResidualEvaluation initial_eval = evaluate_pure_ion_born_residual_jacobian_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        start
    );
    vector<double> theta = start;
    {
        ceres::Problem problem;
        auto *cost = new PureIonCeresCostFunction(
            base_args_by_record,
            records,
            target_kinds,
            target_indices,
            target_indices_2
        );
        problem.AddResidualBlock(cost, nullptr, theta.data());
        for (int j = 0; j < static_cast<int>(theta.size()); ++j) {
            problem.SetParameterLowerBound(theta.data(), j, lower[static_cast<size_t>(j)]);
            problem.SetParameterUpperBound(theta.data(), j, upper[static_cast<size_t>(j)]);
        }

        ceres::Solver::Options options = ceres_regression_options_cpp(max_nfev);

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        PureIonBornResidualEvaluation final_eval = evaluate_pure_ion_born_residual_jacobian_cpp(
            base_args_by_record,
            records,
            target_kinds,
            target_indices,
            target_indices_2,
            theta
        );
        GenericRegressionResult out;
        out.x = theta;
        out.cost = final_eval.cost;
        out.residual_norm = final_eval.residual_norm;
        out.initial_cost = initial_eval.cost;
        out.initial_residual_norm = initial_eval.residual_norm;
        out.metrics_by_term = final_eval.metrics_by_term;
        out.success = summary.IsSolutionUsable() || final_eval.cost <= initial_eval.cost + 1.0e-12;
        out.status = static_cast<int>(summary.termination_type);
        out.message = summary.BriefReport();
        out.nfev = static_cast<int>(summary.num_residual_evaluations + summary.num_jacobian_evaluations);
        out.iterations = static_cast<int>(summary.iterations.size());
        out.backend = "ceres";
        out.optimizer_backend = "ceres";
        out.derivative_backend = "cppad_implicit";
        out.jacobian_available = true;
        out.jacobian_backend = "cppad_implicit";
        return out;
    }
}

// AlgID: binary_kij_ceres_regression
class BinaryKijCeresCostFunction final : public ceres::CostFunction {
public:
    BinaryKijCeresCostFunction(
        vector<add_args> base_args_by_record,
        vector<GenericRegressionRecord> records,
        vector<int> target_kinds,
        vector<int> target_indices,
        vector<int> target_indices_2
    )
        : base_args_by_record_(std::move(base_args_by_record)),
          records_(std::move(records)),
          target_kinds_(std::move(target_kinds)),
          target_indices_(std::move(target_indices)),
          target_indices_2_(std::move(target_indices_2))
    {
        set_num_residuals(generic_residual_count_from_records_cpp(records_));
        mutable_parameter_block_sizes()->push_back(1);
    }

    bool Evaluate(double const *const *parameters, double *residuals, double **jacobians) const override {
        try {
            vector<double> theta = {parameters[0][0]};
            BinaryKijResidualEvaluation eval = evaluate_binary_kij_residual_jacobian_cpp(
                base_args_by_record_,
                records_,
                target_kinds_,
                target_indices_,
                target_indices_2_,
                theta
            );
            for (Eigen::Index row = 0; row < eval.residuals.size(); ++row) {
                residuals[row] = eval.residuals[row];
            }
            if (jacobians != nullptr && jacobians[0] != nullptr) {
                for (Eigen::Index row = 0; row < eval.jacobian.rows(); ++row) {
                    jacobians[0][row] = eval.jacobian(row, 0);
                }
            }
            return true;
        } catch (...) {
            return false;
        }
    }

private:
    vector<add_args> base_args_by_record_;
    vector<GenericRegressionRecord> records_;
    vector<int> target_kinds_;
    vector<int> target_indices_;
    vector<int> target_indices_2_;
};

GenericRegressionResult solve_one_binary_kij_ceres_start_cpp(
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
    BinaryKijResidualEvaluation initial_eval = evaluate_binary_kij_residual_jacobian_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        start
    );
    vector<double> theta = start;
    {
        ceres::Problem problem;
        auto *cost = new BinaryKijCeresCostFunction(
            base_args_by_record,
            records,
            target_kinds,
            target_indices,
            target_indices_2
        );
        problem.AddResidualBlock(cost, nullptr, theta.data());
        problem.SetParameterLowerBound(theta.data(), 0, lower[0]);
        problem.SetParameterUpperBound(theta.data(), 0, upper[0]);

        ceres::Solver::Options options = ceres_regression_options_cpp(max_nfev);

        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        BinaryKijResidualEvaluation final_eval = evaluate_binary_kij_residual_jacobian_cpp(
            base_args_by_record,
            records,
            target_kinds,
            target_indices,
            target_indices_2,
            theta
        );
        GenericRegressionResult out;
        out.x = theta;
        out.cost = final_eval.cost;
        out.residual_norm = final_eval.residual_norm;
        out.initial_cost = initial_eval.cost;
        out.initial_residual_norm = initial_eval.residual_norm;
        out.metrics_by_term = final_eval.metrics_by_term;
        out.success = summary.IsSolutionUsable() || final_eval.cost <= initial_eval.cost + 1.0e-12;
        out.status = static_cast<int>(summary.termination_type);
        out.message = summary.BriefReport();
        out.nfev = static_cast<int>(summary.num_residual_evaluations + summary.num_jacobian_evaluations);
        out.iterations = static_cast<int>(summary.iterations.size());
        out.backend = "ceres";
        out.optimizer_backend = "ceres";
        out.derivative_backend = "cppad_implicit";
        out.jacobian_available = true;
        out.jacobian_backend = "cppad_implicit";
        return out;
    }
}
#endif

}  // namespace regression_detail

using regression_detail::PureNeutralObjectiveEvaluation;
using regression_detail::PureNeutralResidualEvaluation;
using regression_detail::RegressionProfilingStats;
using regression_detail::canonical_start_cpp;
using regression_detail::evaluate_pure_neutral_objective_cpp;
using regression_detail::evaluate_residual_jacobian_cpp;
using regression_detail::kThetaSize;
using regression_detail::objective_from_residual_eval_cpp;
using regression_detail::validate_pure_neutral_base_args_cpp;
#ifdef EPCSAFT_HAS_CERES
using regression_detail::solve_one_start_ceres_cpp;
#endif
using regression_detail::evaluate_generic_residuals_cpp;
using regression_detail::kGenericTargetDBorn;
using regression_detail::kGenericTargetDielc;
using regression_detail::kGenericTargetE;
using regression_detail::kGenericTargetFSolv;
using regression_detail::kGenericTargetKIJ;
using regression_detail::kGenericTargetS;
#ifdef EPCSAFT_HAS_CERES
using regression_detail::solve_one_binary_kij_ceres_start_cpp;
using regression_detail::solve_one_pure_ion_ceres_start_cpp;
#endif

epcsaft::native::cppad_support::CppADDerivativeResult cppad_pure_neutral_parameter_derivatives_cpp(
    double t,
    double rho,
    const add_args &base_args
) {
    return regression_detail::cppad_pure_neutral_parameter_derivatives_cpp(t, rho, base_args);
}

PureNeutralRegressionDebugResult evaluate_pure_neutral_objective_debug_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x
) {
    validate_pure_neutral_base_args_cpp(base_args);
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
    out.gradient.assign(eval.gradient.begin(), eval.gradient.end());
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

PureNeutralRegressionResult fit_pure_neutral_ceres_cpp(
    const add_args &base_args,
    const vector<PureNeutralRegressionDensityRecord> &density_records,
    double density_scale,
    const vector<PureNeutralRegressionVLERecord> &pure_vle_records,
    double pure_vle_scale,
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper
) {
    validate_pure_neutral_base_args_cpp(base_args);
    if (density_records.empty() || pure_vle_records.empty()) {
        throw ValueError("Native Ceres pure-neutral regression requires both density and pure-VLE record families.");
    }
    if (x0.size() != kThetaSize || lower.size() != kThetaSize || upper.size() != kThetaSize) {
        throw ValueError("Native Ceres pure-neutral regression requires 3-variable starts and bounds for m, s, and e.");
    }
    PureNeutralRegressionResult result = solve_one_start_ceres_cpp(
        base_args,
        density_records,
        density_scale,
        pure_vle_records,
        pure_vle_scale,
        canonical_start_cpp(x0, lower, upper),
        lower,
        upper
    );
    result.starts_tried = 1;
    result.backend = "ceres";
    result.optimizer_backend = "ceres";
    result.derivative_backend = "cppad_implicit";
    result.jacobian_backend = result.derivative_backend;
    return result;
}

GenericRegressionDebugResult evaluate_generic_regression_debug_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &x
) {
    return evaluate_generic_residuals_cpp(
        base_args_by_record,
        records,
        target_kinds,
        target_indices,
        target_indices_2,
        x
    );
}

GenericRegressionResult fit_generic_ceres_cpp(
    const vector<add_args> &base_args_by_record,
    const vector<GenericRegressionRecord> &records,
    const vector<int> &target_kinds,
    const vector<int> &target_indices,
    const vector<int> &target_indices_2,
    const vector<double> &x0,
    const vector<double> &lower,
    const vector<double> &upper,
    int max_nfev
) {
    if (target_kinds.empty()) {
        throw ValueError("Native Ceres generic regression requires at least one optimization target.");
    }
    if (x0.size() != target_kinds.size() || lower.size() != target_kinds.size() || upper.size() != target_kinds.size()) {
        throw ValueError("Native Ceres generic regression target arrays must have matching lengths.");
    }
    if (max_nfev < 1) {
        throw ValueError("Native Ceres generic regression requires max_nfev >= 1.");
    }
    const bool is_binary_kij = target_kinds.size() == 1 && target_kinds[0] == kGenericTargetKIJ;
    bool is_pure_ion_parameter_set = !target_kinds.empty();
    for (int kind : target_kinds) {
        if (kind != kGenericTargetS
            && kind != kGenericTargetE
            && kind != kGenericTargetDBorn
            && kind != kGenericTargetFSolv
            && kind != kGenericTargetDielc) {
            is_pure_ion_parameter_set = false;
            break;
        }
    }
    if (!is_binary_kij && !is_pure_ion_parameter_set) {
        throw ValueError("unsupported: native Ceres generic regression has no native analytic/CppAD/implicit derivative path for this target set.");
    }
    vector<double> start = canonical_start_cpp(x0, lower, upper);
    GenericRegressionResult result = is_binary_kij
        ? solve_one_binary_kij_ceres_start_cpp(
            base_args_by_record,
            records,
            target_kinds,
            target_indices,
            target_indices_2,
            start,
            lower,
            upper,
            max_nfev
        )
        : solve_one_pure_ion_ceres_start_cpp(
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
    result.starts_tried = 1;
    result.backend = "ceres";
    result.optimizer_backend = "ceres";
    result.derivative_backend = "cppad_implicit";
    result.jacobian_available = true;
    result.jacobian_backend = "cppad_implicit";
    return result;
}

