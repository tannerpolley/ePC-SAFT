using Clapeyron
using ForwardDiff
using JSON

const R_GAS = 8.31446261815324
const T_REF = 298.15
const P_REF = 1.0e5
const EPS = 1.0e-8
const EPS_INF = 1.0e-12
const FD_H = 1.0e-8

const SCRIPT_DIR = dirname(@__FILE__)
const PACKAGE_DIR = dirname(SCRIPT_DIR)
const DIELECTRIC_OVERRIDE = joinpath(PACKAGE_DIR, "clapeyron_dielectric_overrides.csv")

function require_args()
    if length(ARGS) < 2
        error("Usage: julia --project=<clapeyron_root> extract_clapeyron_na_component_breakdown.jl <output_json> <clapeyron_root>")
    end
    return ARGS[1], ARGS[2]
end

function neutral_userlocations(clapeyron_root::String)
    return [
        joinpath(clapeyron_root, "database", "SAFT", "PCSAFT", "ePCSAFTAdv", "ePCSAFTAdv_like.csv"),
        joinpath(clapeyron_root, "database", "SAFT", "PCSAFT", "ePCSAFTAdv", "ePCSAFTAdv_unlike.csv"),
        joinpath(clapeyron_root, "database", "SAFT", "PCSAFT", "ePCSAFT", "ePCSAFT_like.csv"),
        joinpath(clapeyron_root, "database", "SAFT", "PCSAFT", "ePCSAFT", "ePCSAFT_unlike.csv"),
    ]
end

function build_model(clapeyron_root::String)
    return ESElectrolyte(
        ["water"],
        ["sodium", "chloride"];
        neutralmodel = pharmaPCSAFT,
        ionmodel = DHBorn,
        RSPmodel = LinMixRSP,
        neutralmodel_userlocations = neutral_userlocations(clapeyron_root),
        RSPmodel_userlocations = [DIELECTRIC_OVERRIDE],
        verbose = false,
    )
end

function state_tuple(model)
    z_bulk = [1.0 - 2.0 * EPS, EPS, EPS]
    V_bulk = Clapeyron.volume(model, P_REF, T_REF, z_bulk, phase = :l)
    z_ref = [1.0, 0.0, 0.0]
    p_ref = Clapeyron.pressure(model, V_bulk, T_REF, z_ref)
    z_inf = copy(z_ref)
    z_inf[2] = EPS_INF
    z_inf[1] = 1.0 - EPS_INF
    V_inf = Clapeyron.volume(model, p_ref, T_REF, z_inf, phase = :l)
    return z_inf, Float64(p_ref), Float64(V_inf)
end

function molar_a_hc(model, V, T, z)
    return R_GAS * T * Clapeyron.a_hc(model.neutralmodel, V, T, z) / 1000.0
end

function molar_a_disp(model, V, T, z)
    return R_GAS * T * Clapeyron.a_disp(model.neutralmodel, V, T, z) / 1000.0
end

function molar_a_assoc(model, V, T, z)
    return R_GAS * T * Clapeyron.a_assoc(model.neutralmodel, V, T, z) / 1000.0
end

function ion_data(model, V, T, z)
    return Clapeyron.iondata(model, V, T, z)
end

function molar_a_dh(model, V, T, z)
    return R_GAS * T * Clapeyron.a_dh(model.ionmodel, V, T, z, ion_data(model, V, T, z)) / 1000.0
end

function molar_a_born(model, V, T, z)
    return R_GAS * T * Clapeyron.a_born(model.ionmodel, V, T, z, ion_data(model, V, T, z)) / 1000.0
end

function eos_from_molar(fn)
    return (model, V, T, z) -> 1000.0 * fn(model, V, T, z) * sum(z)
end

function pressure_from_molar_fn(model, V, T, z, fn)
    energy_fn = eos_from_molar(fn)
    return -ForwardDiff.derivative(v -> energy_fn(model, v, T, z), V)
end

function simplex_fd(model, fn, V, T, z, idx, bal_idx; h = FD_H)
    z2 = copy(z)
    z2[idx] += h
    z2[bal_idx] -= h
    return (fn(model, V, T, z2) - fn(model, V, T, z)) / h
end

function diag_values(values)
    if ndims(values) == 1
        return [Float64(values[i]) for i in eachindex(values)]
    end
    return [Float64(values[i, i]) for i in 1:size(values, 1)]
end

function contribution_payload(model, z, p_ref, V_inf, name, fn)
    mu = Clapeyron.VT_molar_gradient(model, V_inf, T_REF, z, eos_from_molar(fn)) ./ 1000.0
    a = fn(model, V_inf, T_REF, z)
    p_alpha = pressure_from_molar_fn(model, V_inf, T_REF, z, fn)
    z_total = p_ref * V_inf / (sum(z) * R_GAS * T_REF)
    z_residual = z_total - 1.0
    z_alpha = p_alpha * V_inf / (sum(z) * R_GAS * T_REF)
    z_term = (R_GAS * T_REF / 1000.0) * z_alpha
    z_correction = (R_GAS * T_REF / 1000.0) * (-(z_alpha / z_residual) * log(z_total))
    return Dict(
        "name" => name,
        "a_kj_mol" => Float64(a),
        "z_kj_mol" => Float64(z_term),
        "a_plus_z_kj_mol" => Float64(a + z_term),
        "z_correction_kj_mol" => Float64(z_correction),
        "mu_components_kj_mol" => [Float64(value) for value in mu],
        "lnfug_target_kj_mol" => Float64(mu[2] + z_correction),
        "weighted_mu_kj_mol" => Float64(sum(z .* mu)),
        "mu_minus_a_z_components_kj_mol" => [Float64(value - (a + z_term)) for value in mu],
        "simplex_fd_target_kj_mol" => Float64(simplex_fd(model, fn, V_inf, T_REF, z, 2, 1)),
        "simplex_fd_counter_kj_mol" => Float64(simplex_fd(model, fn, V_inf, T_REF, z, 3, 1)),
        "simplex_fd_water_kj_mol" => Float64(simplex_fd(model, fn, V_inf, T_REF, z, 1, 2)),
    )
end

function assoc_payload(model, z, p_ref, V_inf)
    payload = contribution_payload(model, z, p_ref, V_inf, "assoc", molar_a_assoc)
    neutral = model.neutralmodel
    neutral_data = Clapeyron.data(neutral, V_inf, T_REF, z)
    X, Δ = Clapeyron.X_and_Δ(neutral, V_inf, T_REF, z, neutral_data)
    K = Clapeyron.assoc_site_matrix(neutral, V_inf, T_REF, z, neutral_data, Δ)
    payload["assoc_X_values"] = [Float64(value) for value in X.v]
    payload["assoc_delta_values"] = [Float64(value) for value in Δ.values]
    payload["assoc_site_matrix"] = [[Float64(K[i, j]) for j in axes(K, 2)] for i in axes(K, 1)]
    return payload
end

function sanitize_json_value(x)
    if x isa AbstractFloat
        return isfinite(x) ? Float64(x) : nothing
    elseif x isa AbstractVector
        return [sanitize_json_value(v) for v in x]
    elseif x isa AbstractMatrix
        return [[sanitize_json_value(x[i, j]) for j in axes(x, 2)] for i in axes(x, 1)]
    elseif x isa Dict
        return Dict(String(k) => sanitize_json_value(v) for (k, v) in x)
    else
        return x
    end
end

function main()
    output_json, clapeyron_root = require_args()
    model = build_model(clapeyron_root)
    z, p_ref, V_inf = state_tuple(model)
    payload = Dict(
        "package" => "Clapeyron.jl",
        "config" => "mixed_current",
        "components" => model.components,
        "target_component" => "sodium",
        "counter_component" => "chloride",
        "z" => [Float64(value) for value in z],
        "density_mol_m3" => Float64(sum(z) / V_inf),
        "reference_pressure_pa" => p_ref,
        "state_volume_m3" => V_inf,
        "compressibility_factor" => Float64(p_ref * V_inf / (sum(z) * R_GAS * T_REF)),
        "segment_diag" => diag_values(model.neutralmodel.params.segment.values),
        "sigma_diag_m" => diag_values(model.neutralmodel.params.sigma.values),
        "epsilon_diag_k" => diag_values(model.neutralmodel.params.epsilon.values),
        "k_matrix" => [[Float64(model.neutralmodel.params.k.values[i, j]) for j in 1:3] for i in 1:3],
        "sigma_born_m" => [Float64(value) for value in model.ionmodel.params.sigma_born.values],
        "dielectric_loaded" => [Float64(value) for value in model.ionmodel.RSPmodel.params.dielectric_constant.values],
        "neutralmodel_userlocations" => neutral_userlocations(clapeyron_root),
        "RSPmodel_userlocations" => [DIELECTRIC_OVERRIDE],
        "contributions" => Dict(
            "hc" => contribution_payload(model, z, p_ref, V_inf, "hc", molar_a_hc),
            "disp" => contribution_payload(model, z, p_ref, V_inf, "disp", molar_a_disp),
            "assoc" => assoc_payload(model, z, p_ref, V_inf),
            "dh" => contribution_payload(model, z, p_ref, V_inf, "dh", molar_a_dh),
            "born" => contribution_payload(model, z, p_ref, V_inf, "born", molar_a_born),
        ),
    )
    mkpath(dirname(output_json))
    open(output_json, "w") do io
        JSON.print(io, sanitize_json_value(payload), 2)
    end
    println("Wrote $(output_json)")
end

main()
