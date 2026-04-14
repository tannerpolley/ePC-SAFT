using Clapeyron
using JSON

const R_GAS = 8.31446261815324
const T_REF = 298.15
const P_REF = 1.0e5
const EPS = 1.0e-8
const EPS_INF = 1.0e-12

const SCRIPT_DIR = dirname(@__FILE__)
const PACKAGE_DIR = dirname(SCRIPT_DIR)
const REPO_OVERRIDE = joinpath(PACKAGE_DIR, "clapeyron_repo_unlike_overrides.csv")
const DIELECTRIC_OVERRIDE = joinpath(PACKAGE_DIR, "clapeyron_dielectric_overrides.csv")

const ION_SETUPS = Dict(
    "Li+" => Dict("ions" => ["lithium", "chloride"], "target_index" => 2, "counter_index" => 3),
    "Na+" => Dict("ions" => ["sodium", "chloride"], "target_index" => 2, "counter_index" => 3),
    "K+" => Dict("ions" => ["potassium", "chloride"], "target_index" => 2, "counter_index" => 3),
    "F-" => Dict("ions" => ["sodium", "fluoride"], "target_index" => 3, "counter_index" => 2),
    "Cl-" => Dict("ions" => ["sodium", "chloride"], "target_index" => 3, "counter_index" => 2),
    "Br-" => Dict("ions" => ["sodium", "bromide"], "target_index" => 3, "counter_index" => 2),
    "I-" => Dict("ions" => ["sodium", "iodide"], "target_index" => 3, "counter_index" => 2),
)

function require_args()
    if length(ARGS) < 2
        error("Usage: julia --project=<clapeyron_root> scan_clapeyron_dispersion_configs.jl <output_json> <clapeyron_root>")
    end
    return ARGS[1], ARGS[2]
end

function config_files(clapeyron_root::String)
    adv_like = joinpath(clapeyron_root, "database", "SAFT", "ePCSAFT", "ePCSAFTAdv", "ePCSAFTAdv_like.csv")
    adv_unlike = joinpath(clapeyron_root, "database", "SAFT", "ePCSAFT", "ePCSAFTAdv", "ePCSAFTAdv_unlike.csv")
    rev_like = joinpath(clapeyron_root, "database", "SAFT", "ePCSAFT", "ePCSAFT", "ePCSAFT_like.csv")
    rev_unlike = joinpath(clapeyron_root, "database", "SAFT", "ePCSAFT", "ePCSAFT", "ePCSAFT_unlike.csv")
    return Dict(
        "mixed_current" => [adv_like, adv_unlike, rev_like, rev_unlike],
        "advanced_only" => [adv_like, adv_unlike],
        "repo_unlike_override" => [adv_like, adv_unlike, rev_like, rev_unlike, REPO_OVERRIDE],
    )
end

function build_model(ions::Vector{String}, neutral_files::Vector{String})
    return ESElectrolyte(
        ["water"],
        ions;
        neutralmodel = pharmaePCSAFT,
        ionmodel = DHBorn,
        RSPmodel = LinMixRSP,
        neutralmodel_userlocations = neutral_files,
        RSPmodel_userlocations = [DIELECTRIC_OVERRIDE],
        verbose = false,
    )
end

function eos_disp(model, V, T, z)
    return Clapeyron.Rgas(model) * T * sum(z) * Clapeyron.a_disp(model.neutralmodel, V, T, z)
end

function dispersion_mu(model, target_index::Int)
    z_bulk = [1.0 - 2.0 * EPS, EPS, EPS]
    V_bulk = Clapeyron.volume(model, P_REF, T_REF, z_bulk, phase = :l)
    z_ref = [1.0, 0.0, 0.0]
    p_ref = Clapeyron.pressure(model, V_bulk, T_REF, z_ref)
    z_inf = copy(z_ref)
    z_inf[target_index] = EPS_INF
    z_inf[1] = 1.0 - EPS_INF
    V_inf = Clapeyron.volume(model, p_ref, T_REF, z_inf, phase = :l)
    grad = Clapeyron.VT_molar_gradient(model, V_inf, T_REF, z_inf, eos_disp)
    return Float64(grad[target_index] / 1000.0), Float64(p_ref), Float64(V_inf)
end

function diag_values(matrix)
    n = size(matrix, 1)
    return [Float64(matrix[i, i]) for i in 1:n]
end

function collect_one(ion::String, setup::Dict{String, Any}, neutral_files::Vector{String})
    ions = Vector{String}(setup["ions"])
    target_index = Int(setup["target_index"])
    counter_index = Int(setup["counter_index"])
    model = build_model(ions, neutral_files)
    disp_kj_mol, p_ref, V_inf = dispersion_mu(model, target_index)
    k = model.neutralmodel.params.k.values
    return Dict(
        "components" => model.components,
        "disp_kj_mol" => disp_kj_mol,
        "reference_pressure_pa" => p_ref,
        "state_volume_m3" => V_inf,
        "k_matrix" => [[Float64(k[i, j]) for j in axes(k, 2)] for i in axes(k, 1)],
        "sigma_diag_m" => diag_values(model.neutralmodel.params.sigma.values),
        "epsilon_diag_k" => diag_values(model.neutralmodel.params.epsilon.values),
        "water_counter_k" => Float64(k[1, counter_index]),
        "water_target_k" => Float64(k[1, target_index]),
        "counter_target_k" => Float64(k[counter_index, target_index]),
    )
end

function main()
    output_json, clapeyron_root = require_args()
    results = Dict{String, Any}()
    for (config_name, neutral_files) in config_files(clapeyron_root)
        config_results = Dict{String, Any}()
        for ion in ("Li+", "Na+", "K+", "F-", "Cl-", "Br-", "I-")
            try
                config_results[ion] = collect_one(ion, ION_SETUPS[ion], neutral_files)
            catch err
                config_results[ion] = Dict(
                    "error" => sprint(showerror, err),
                    "neutralmodel_userlocations" => neutral_files,
                )
            end
        end
        results[config_name] = Dict(
            "neutralmodel_userlocations" => neutral_files,
            "results" => config_results,
        )
    end
    payload = Dict(
        "package" => "Clapeyron.jl",
        "RSPmodel_userlocations" => [DIELECTRIC_OVERRIDE],
        "configs" => results,
    )
    mkpath(dirname(output_json))
    open(output_json, "w") do io
        JSON.print(io, payload, 2)
    end
    println("Wrote $(output_json)")
end

main()

