using Clapeyron
using JSON

function require_args()
    if length(ARGS) != 2
        error("Usage: julia external_package_case2_clapeyron.jl <input_json> <output_json>")
    end
    return ARGS[1], ARGS[2]
end

function err_string(err)
    io = IOBuffer()
    showerror(io, err)
    return String(take!(io))
end

function as_matrix_rows(x)
    return [collect(row) for row in eachrow(Matrix(x))]
end

function as_vector(x)
    return collect(x)
end

function build_model(neutrals, ions; charges=nothing)
    if charges === nothing
        return ePCSAFT(neutrals, ions)
    end
    return ePCSAFT(neutrals, ions, charge=charges)
end

function build_result(neutrals, ions; charges=nothing)
    try
        model = build_model(neutrals, ions; charges=charges)
        return Dict(
            "success" => true,
            "components" => collect(model.components),
        )
    catch err
        return Dict(
            "success" => false,
            "error" => err_string(err),
        )
    end
end

function build_flash_attempt(label, K0, model, p, T, feed, alcohol_index, charges)
    method = MichelsenTPFlash(equilibrium=:lle, K0=K0)
    x, nph, G = tp_flash(model, p, T, feed, method)

    xmat = Matrix(x)
    nmat = Matrix(nph)
    phase_totals = vec(sum(nmat; dims=2))
    betas = phase_totals ./ sum(phase_totals)
    organic_row = xmat[1, alcohol_index] >= xmat[2, alcohol_index] ? 1 : 2
    aqueous_row = organic_row == 1 ? 2 : 1
    org_x = vec(xmat[organic_row, :])
    aq_x = vec(xmat[aqueous_row, :])
    org_phi = as_vector(fugacity_coefficient(model, p, T, org_x; phase=:l))
    aq_phi = as_vector(fugacity_coefficient(model, p, T, aq_x; phase=:l))
    log_p_bar = log(p / 1.0e5)
    org_lnfug_bar = log.(org_phi) .+ log.(org_x) .+ log_p_bar
    aq_lnfug_bar = log.(aq_phi) .+ log.(aq_x) .+ log_p_bar
    charge_residuals = vec(xmat * charges)

    return Dict(
        "success" => true,
        "method" => label,
        "K0" => collect(K0),
        "phase_molefracs" => as_matrix_rows(xmat),
        "phase_component_moles" => as_matrix_rows(nmat),
        "phase_total_moles" => collect(phase_totals),
        "phase_betas" => collect(betas),
        "organic_phase_index" => organic_row,
        "aqueous_phase_index" => aqueous_row,
        "organic_lnfug_bar" => collect(org_lnfug_bar),
        "aqueous_lnfug_bar" => collect(aq_lnfug_bar),
        "organic_lnphi" => collect(log.(org_phi)),
        "aqueous_lnphi" => collect(log.(aq_phi)),
        "phase_charge_residuals" => collect(charge_residuals),
        "gibbs_energy" => G,
    )
end

function flash_result(neutrals, ions, charges, p, T, feed, method_configs)
    try
        model = build_model(neutrals, ions; charges=charges)
        attempts = Any[]
        alcohol_index = 2
        for method_config in method_configs
            label = String(method_config["label"])
            K0 = Float64.(method_config["K0"])
            try
                payload = build_flash_attempt(label, K0, model, p, T, feed, alcohol_index, charges)
                payload["components"] = collect(model.components)
                return payload
            catch err
                push!(attempts, Dict("method" => label, "error" => err_string(err)))
            end
        end
        return Dict(
            "success" => false,
            "components" => collect(model.components),
            "attempts" => attempts,
        )
    catch err
        return Dict(
            "success" => false,
            "error" => err_string(err),
        )
    end
end

function main()
    input_json, output_json = require_args()
    payload = JSON.parsefile(input_json)
    neutrals = Vector{String}(payload["neutral_components"])
    ions = Vector{String}(payload["ion_components"])
    charges = Int.(payload["charges"])
    result = Dict(
        "case_key" => payload["case_key"],
        "label" => payload["label"],
        "build_without_charge" => build_result(neutrals, ions),
        "build_with_charge" => build_result(neutrals, ions; charges=charges),
    )
    if payload["action"] == "flash"
        feed = Float64.(payload["feed_moles"])
        result["flash"] = flash_result(
            neutrals,
            ions,
            charges,
            Float64(payload["pressure_pa"]),
            Float64(payload["temperature_k"]),
            feed,
            payload["method_configs"],
        )
    end

    mkpath(dirname(output_json))
    open(output_json, "w") do io
        JSON.print(io, result, 2)
    end
end

main()
