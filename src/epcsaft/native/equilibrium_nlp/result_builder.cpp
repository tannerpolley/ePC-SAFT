#include "result_builder.h"

#include "eos_phase_block.h"
#include "epcsaft_core_internal.h"

#include <cmath>
#include <numeric>

namespace epcsaft::native::equilibrium_nlp {

namespace {

void require_positive_finite(double value, const std::string& label) {
    if (std::isfinite(value) && value > 0.0) {
        return;
    }
    throw ValueError(label + " must be positive and finite.");
}

std::vector<double> ln_fugacity_coefficients_for_block(
    const add_args& args,
    const EosPhaseBlockResult& block
) {
    FugacityContributionResult fugacity = fugacity_coefficient_result_cpp(
        block.temperature,
        block.density,
        block.composition,
        args
    );
    if (fugacity.lnfugcoef.total.size() != block.composition.size()) {
        throw ValueError("Neutral EOS result-builder fugacity payload size did not match composition.");
    }
    return fugacity.lnfugcoef.total;
}

NeutralTwoPhaseEosPhasePayload phase_payload(
    const add_args& args,
    const EosPhaseBlockResult& block,
    std::size_t phase_index,
    double total_amount
) {
    require_positive_finite(total_amount, "Neutral EOS result-builder total amount");
    NeutralTwoPhaseEosPhasePayload out;
    out.label = "phase_" + std::to_string(phase_index);
    out.composition = block.composition;
    out.ln_fugacity_coefficient = ln_fugacity_coefficients_for_block(args, block);
    out.fugacity_coefficient = out.ln_fugacity_coefficient;
    out.density = block.density;
    out.temperature = block.temperature;
    out.pressure = block.target_pressure;
    out.phase_fraction = block.total_amount / total_amount;
    out.amount_total = block.total_amount;
    out.volume = block.volume;
    out.eos_pressure = block.eos_pressure;
    out.pressure_consistency_residual = block.pressure_consistency_residual;
    out.compressibility_factor = block.compressibility_factor;
    return out;
}

}  // namespace

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
) {
    NeutralTwoPhaseEosPostsolve postsolve = evaluate_neutral_two_phase_eos_postsolve(
        args,
        temperature,
        target_pressure,
        phase_amounts,
        volumes,
        feed_amounts,
        material_tolerance,
        pressure_tolerance,
        chemical_potential_tolerance,
        phase_distance_tolerance
    );
    NeutralTwoPhaseEosResultPayload out;
    out.accepted = postsolve.accepted;
    out.split_detected = postsolve.accepted;
    out.derivative_backend = postsolve.derivative_backend;
    out.rejection_reason = postsolve.rejection_reason;
    out.objective = postsolve.objective;
    out.material_balance_norm = postsolve.material_balance_norm;
    out.pressure_consistency_norm = postsolve.pressure_consistency_norm;
    out.chemical_potential_consistency_norm = postsolve.chemical_potential_consistency_norm;
    out.phase_distance = postsolve.phase_distance;
    out.feed_amounts = feed_amounts;
    out.constraints = postsolve.constraints;

    if (!postsolve.accepted) {
        return out;
    }
    const EosPhaseSystemResult system = evaluate_eos_phase_system(
        args,
        temperature,
        target_pressure,
        phase_amounts,
        volumes,
        feed_amounts,
        {}
    );
    if (system.phase_count != 2 || system.phase_blocks.size() != 2) {
        throw ValueError("Neutral EOS result-builder currently requires exactly two phase blocks.");
    }
    const double total_amount = std::accumulate(
        postsolve.phase_amount_totals.begin(),
        postsolve.phase_amount_totals.end(),
        0.0
    );
    out.phases.reserve(system.phase_blocks.size());
    for (std::size_t index = 0; index < system.phase_blocks.size(); ++index) {
        out.phases.push_back(phase_payload(args, system.phase_blocks[index], index, total_amount));
    }
    return out;
}

}  // namespace epcsaft::native::equilibrium_nlp
