using Clapeyron
using ForwardDiff
using JSON

const R_GAS = 8.31446261815324
const T_REF = 298.15
const P_REF = 1.0e5
const EPS = 1.0e-8
const EPS_INF = 1.0e-12
const TERM_KEYS = ("hc", "disp", "assoc", "dh", "born")

const ION_SETUPS = Dict(
    "Li+" => Dict("ions" => ["lithium", "chloride"], "target_index" => 2),
    "Na+" => Dict("ions" => ["sodium", "chloride"], "target_index" => 2),
    "K+" => Dict("ions" => ["potassium", "chloride"], "target_index" => 2),
    "F-" => Dict("ions" => ["sodium", "fluoride"], "target_index" => 3),
    "Cl-" => Dict("ions" => ["sodium", "chloride"], "target_index" => 3),
    "Br-" => Dict("ions" => ["sodium", "bromide"], "target_index" => 3),
    "I-" => Dict("ions" => ["sodium", "iodide"], "target_index" => 3),
)

function require_args()
    if length(ARGS) < 2
        error("Usage: julia --project=<clapeyron_root> extract_clapeyron_figure3.jl <output_json> <clapeyron_root>")
    end
    return ARGS[1], ARGS[2]
end

function neutral_userlocations(clapeyron_root::String)
    return [
        joinpath(clapeyron_root, "database", "SAFT", "ePCSAFT", "ePCSAFTAdv", "ePCSAFTAdv_like.csv"),
        joinpath(clapeyron_root, "database", "SAFT", "ePCSAFT", "ePCSAFTAdv", "ePCSAFTAdv_unlike.csv"),
        joinpath(clapeyron_root, "database", "SAFT", "ePCSAFT", "ePCSAFT", "ePCSAFT_like.csv"),
        joinpath(clapeyron_root, "database", "SAFT", "ePCSAFT", "ePCSAFT", "ePCSAFT_unlike.csv"),
    ]
end

function rsp_userlocations()
    return [joinpath(dirname(@__FILE__), "clapeyron_dielectric_overrides.csv")]
end

function build_model(ions::Vector{String}, clapeyron_root::String)
    return ESElectrolyte(
        ["water"],
        ions;
        neutralmodel = pharmaePCSAFT,
        ionmodel = DHBorn,
        RSPmodel = LinMixRSP,
        neutralmodel_userlocations = neutral_userlocations(clapeyron_root),
        RSPmodel_userlocations = rsp_userlocations(),
        verbose = false,
    )
end

function eos_hc(model, V, T, z)
    return Clapeyron.Rgas(model) * T * sum(z) * Clapeyron.a_hc(model.neutralmodel, V, T, z)
end

function eos_disp(model, V, T, z)
    return Clapeyron.Rgas(model) * T * sum(z) * Clapeyron.a_disp(model.neutralmodel, V, T, z)
end

function eos_assoc(model, V, T, z)
    return Clapeyron.Rgas(model) * T * sum(z) * Clapeyron.a_assoc(model.neutralmodel, V, T, z)
end

function eos_dh(model, V, T, z)
    return Clapeyron.Rgas(model) * T * sum(z) * Clapeyron.a_dh(model.ionmodel, V, T, z, Clapeyron.iondata(model, V, T, z))
end

function eos_born(model, V, T, z)
    return Clapeyron.Rgas(model) * T * sum(z) * Clapeyron.a_born(model.ionmodel, V, T, z, Clapeyron.iondata(model, V, T, z))
end

function branch_functions()
    return Dict(
        "hc" => eos_hc,
        "disp" => eos_disp,
        "assoc" => eos_assoc,
        "dh" => eos_dh,
        "born" => eos_born,
    )
end

function pressure_from_energy_fn(model, V, T, z, fn)
    return -ForwardDiff.derivative(v -> fn(model, v, T, z), V)
end

function compute_one(ion::String, setup::Dict{String, Any}, clapeyron_root::String)
    ions = Vector{String}(setup["ions"])
    target_index = Int(setup["target_index"])
    model = build_model(ions, clapeyron_root)

    z_bulk = [1.0 - 2.0 * EPS, EPS, EPS]
    V_bulk = Clapeyron.volume(model, P_REF, T_REF, z_bulk, phase = :l)

    z_ref = [1.0, 0.0, 0.0]
    p_ref = Clapeyron.pressure(model, V_bulk, T_REF, z_ref)

    z_inf = copy(z_ref)
    z_inf[target_index] = EPS_INF
    z_inf[1] = 1.0 - EPS_INF
    V_inf = Clapeyron.volume(model, p_ref, T_REF, z_inf, phase = :l)

    mu_total = Clapeyron.VT_chemical_potential_res(model, V_inf, T_REF, z_inf)[target_index] / 1000.0
    compressibility = p_ref * V_inf / (sum(z_inf) * R_GAS * T_REF)
    total = mu_total - (R_GAS * T_REF * log(compressibility) / 1000.0)
    z_residual = compressibility - 1.0

    terms = Dict{String, Float64}()
    z_corrections = Dict{String, Float64}()
    lnfug_terms = Dict{String, Float64}()
    for (term, fn) in branch_functions()
        values = Clapeyron.VT_molar_gradient(model, V_inf, T_REF, z_inf, fn)
        terms[term] = Float64(values[target_index] / 1000.0)
        p_alpha = pressure_from_energy_fn(model, V_inf, T_REF, z_inf, fn)
        z_alpha = p_alpha * V_inf / (sum(z_inf) * R_GAS * T_REF)
        z_corrections[term] = Float64((R_GAS * T_REF / 1000.0) * (-(z_alpha / z_residual) * log(compressibility)))
        lnfug_terms[term] = Float64(terms[term] + z_corrections[term])
    end

    mu_sum = sum(terms[term] for term in TERM_KEYS)
    lnfug_sum = sum(lnfug_terms[term] for term in TERM_KEYS)
    if !isfinite(total) || !isfinite(mu_total) || !isfinite(mu_sum)
        error("Non-finite Clapeyron result for $(ion).")
    end
    if !isfinite(lnfug_sum)
        error("Non-finite Clapeyron lnfug-sum result for $(ion).")
    end
    if abs(mu_total - mu_sum) > 1.0e-6
        error("Clapeyron mu_total and mapped mu_sum do not close for $(ion): $(mu_total) vs $(mu_sum)")
    end

    return Dict(
        "hc" => terms["hc"],
        "disp" => terms["disp"],
        "assoc" => terms["assoc"],
        "dh" => terms["dh"],
        "born" => terms["born"],
        "hc_z_correction_kj_mol" => z_corrections["hc"],
        "disp_z_correction_kj_mol" => z_corrections["disp"],
        "assoc_z_correction_kj_mol" => z_corrections["assoc"],
        "dh_z_correction_kj_mol" => z_corrections["dh"],
        "born_z_correction_kj_mol" => z_corrections["born"],
        "hc_lnfug_kj_mol" => lnfug_terms["hc"],
        "disp_lnfug_kj_mol" => lnfug_terms["disp"],
        "assoc_lnfug_kj_mol" => lnfug_terms["assoc"],
        "dh_lnfug_kj_mol" => lnfug_terms["dh"],
        "born_lnfug_kj_mol" => lnfug_terms["born"],
        "mu_total_kj_mol" => mu_total,
        "mu_sum_kj_mol" => mu_sum,
        "lnfug_sum_kj_mol" => lnfug_sum,
        "lnfug_gap_kj_mol" => total - lnfug_sum,
        "total_kj_mol" => total,
        "reference_pressure_pa" => Float64(p_ref),
        "state_volume_m3" => Float64(V_inf),
        "compressibility_factor" => Float64(compressibility),
    )
end

function compute_all(clapeyron_root::String)
    results = Dict{String, Any}()
    for ion in ("Li+", "Na+", "K+", "F-", "Cl-", "Br-", "I-")
        results[ion] = compute_one(ion, ION_SETUPS[ion], clapeyron_root)
    end
    return results
end

function main()
    output_json, clapeyron_root = require_args()
    payload = Dict(
        "package" => "Clapeyron.jl",
        "neutralmodel" => "pharmaePCSAFT",
        "ionmodel" => "DHBorn",
        "RSPmodel" => "LinMixRSP",
        "clapeyron_root" => clapeyron_root,
        "neutralmodel_userlocations" => neutral_userlocations(clapeyron_root),
        "RSPmodel_userlocations" => rsp_userlocations(),
        "results" => compute_all(clapeyron_root),
    )
    mkpath(dirname(output_json))
    open(output_json, "w") do io
        JSON.print(io, payload, 2)
    end
    println("Wrote $(output_json)")
end

main()

