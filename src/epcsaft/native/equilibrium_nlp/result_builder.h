#pragma once

#include "route_builders.h"

#include <string>
#include <vector>

struct add_args;

namespace epcsaft::native::equilibrium_nlp {

struct NeutralTwoPhaseEosPhasePayload {
    std::string label;
    std::vector<double> composition;
    std::vector<double> ln_fugacity_coefficient;
    std::vector<double> fugacity_coefficient;
    double density = 0.0;
    double temperature = 0.0;
    double pressure = 0.0;
    double phase_fraction = 0.0;
    double amount_total = 0.0;
    double volume = 0.0;
    double eos_pressure = 0.0;
    double pressure_consistency_residual = 0.0;
    double compressibility_factor = 0.0;
};

struct NeutralTwoPhaseEosResultPayload {
    bool accepted = false;
    bool stable = false;
    bool split_detected = false;
    std::string backend = "native_equilibrium_nlp";
    std::string problem_kind = "neutral_two_phase_eos";
    std::string derivative_backend;
    std::string rejection_reason;
    double objective = 0.0;
    double material_balance_norm = 0.0;
    double pressure_consistency_norm = 0.0;
    double chemical_potential_consistency_norm = 0.0;
    double phase_distance = 0.0;
    std::vector<double> feed_amounts;
    std::vector<double> constraints;
    std::vector<NeutralTwoPhaseEosPhasePayload> phases;
};

NeutralTwoPhaseEosResultPayload build_neutral_two_phase_eos_result(
    const add_args& args,
    double temperature,
    double target_pressure,
    const std::vector<std::vector<double>>& phase_amounts,
    const std::vector<double>& volumes,
    const std::vector<double>& feed_amounts,
    double material_tolerance,
    double pressure_tolerance,
    double chemical_potential_tolerance,
    double phase_distance_tolerance
);

}  // namespace epcsaft::native::equilibrium_nlp
