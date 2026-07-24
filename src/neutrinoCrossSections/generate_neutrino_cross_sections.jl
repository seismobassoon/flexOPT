"""
    CubicSpline

Knots and polynomial coefficients of a SciPy-generated cubic spline.
"""
struct CubicSpline
    k::Vector{Float64}
    c::Matrix{Float64}

    function CubicSpline(k::Vector{Float64}, c::Matrix{Float64})
        length(k) >= 2 ||
            throw(ArgumentError("a cubic spline needs at least two knots"))
        issorted(k) ||
            throw(ArgumentError("cubic-spline knots must be sorted"))
        size(c, 1) == 4 ||
            throw(DimensionMismatch("cubic-spline coefficients need four rows"))
        size(c, 2) >= length(k) - 1 ||
            throw(DimensionMismatch(
                "not enough coefficient columns for $(length(k)) knots",
            ))
        new(k, c)
    end
end

const CrossSectionSplines = NamedTuple{
    (:νe_CC, :νμ_CC, :ντ_CC, :antiνe_CC, :antiνμ_CC, :antiντ_CC,
     :ν_NC, :antiν_NC, :source),
}

const NeutrinoCrossSectionTable = NamedTuple{
    (:energies, :νe_CC, :νμ_CC, :ντ_CC, :antiνe_CC, :antiνμ_CC,
     :antiντ_CC, :ν_NC, :antiν_NC, :units, :source),
}

const _CROSS_SECTION_JSON_KEYS = (
    νe_CC="nue_CC_logE",
    νμ_CC="num_CC_logE",
    ντ_CC="nut_CC_logE",
    antiνe_CC="anue_CC_logE",
    antiνμ_CC="anum_CC_logE",
    antiντ_CC="anut_CC_logE",
    ν_NC="nu_NC_logE",
    antiν_NC="anu_NC_logE",
)

function _resolve_cross_section_path(path::AbstractString)
    candidates = String[path]
    endswith(lowercase(path), ".json") || push!(candidates, path * ".json")
    data_dir = normpath(
        joinpath(@__DIR__, "..", "..", "dataInput", "neutrinoCrossSections"),
    )
    append!(candidates, joinpath.(Ref(data_dir), basename.(copy(candidates))))
    match = findfirst(isfile, candidates)
    isnothing(match) &&
        throw(ArgumentError(
            "cross-section JSON not found; tried $(join(candidates, ", "))",
        ))
    return abspath(candidates[match])
end

"""
    read_neutrino_cross_sections_info(path="cross_sections.json")

Read all charged-current and neutral-current splines from the cross-section
JSON file. The result is a named tuple; no positional unpacking is required.
"""
function read_neutrino_cross_sections_info(
    path::AbstractString="cross_sections.json",
)
    data_path = _resolve_cross_section_path(path)
    endswith(lowercase(data_path), ".json") ||
        throw(ArgumentError("unsupported cross-section format: $data_path"))
    data = JSON.parsefile(data_path)

    parsed_splines = map(Base.values(_CROSS_SECTION_JSON_KEYS)) do json_key
        haskey(data, json_key) ||
            throw(ArgumentError("cross-section JSON is missing key $json_key"))
        entry = data[json_key]
        knots = Float64.(entry["k"])
        coefficients = permutedims(Float64.(hcat(entry["c"]...)))
        CubicSpline(knots, coefficients)
    end
    return (
        νe_CC=parsed_splines[1],
        νμ_CC=parsed_splines[2],
        ντ_CC=parsed_splines[3],
        antiνe_CC=parsed_splines[4],
        antiνμ_CC=parsed_splines[5],
        antiντ_CC=parsed_splines[6],
        ν_NC=parsed_splines[7],
        antiν_NC=parsed_splines[8],
        source=data_path,
    )
end

"""
    evaluate_cubicspline(spline, x; extrapolate=true)

Evaluate a spline at `x` using Horner's method. With `extrapolate=false`, an
out-of-domain coordinate raises an informative error.
"""
function evaluate_cubicspline(
    spline::CubicSpline,
    x::Real;
    extrapolate::Bool=true,
)
    coordinate = Float64(x)
    if !extrapolate && !(first(spline.k) <= coordinate <= last(spline.k))
        throw(DomainError(
            coordinate,
            "spline coordinate must lie in " *
            "[$(first(spline.k)), $(last(spline.k))]",
        ))
    end
    interval = clamp(
        searchsortedlast(spline.k, coordinate),
        1,
        length(spline.k) - 1,
    )
    dx = coordinate - spline.k[interval]
    coefficients = spline.c
    return (
        (
            coefficients[1, interval] * dx +
            coefficients[2, interval]
        ) * dx + coefficients[3, interval]
    ) * dx + coefficients[4, interval]
end

"""
    evaluate_neutrino_cross_sections(energies;
        path="cross_sections.json", extrapolate=false, clamp_nonnegative=true)

Evaluate every available neutrino cross section on an energy grid in GeV.
Returned arrays follow the same energy axis as `neutrinoFluxTable.energies`.
The tabulated quantity is `(σ/E) × 10⁴² m²/GeV` per nucleon.
Small negative cubic-spline overshoots near interaction thresholds are clamped
to zero by default.
"""
function evaluate_neutrino_cross_sections(
    energies::AbstractVector{<:Real};
    path::AbstractString="cross_sections.json",
    extrapolate::Bool=false,
    clamp_nonnegative::Bool=true,
)
    energy_grid = Float64.(collect(energies))
    isempty(energy_grid) &&
        throw(ArgumentError("cross-section energy grid cannot be empty"))
    all(isfinite, energy_grid) ||
        throw(ArgumentError("cross-section energies must be finite"))
    all(>(0), energy_grid) ||
        throw(ArgumentError("cross-section energies must be positive"))
    all(diff(energy_grid) .> 0) ||
        throw(ArgumentError(
            "cross-section energies must be strictly increasing",
        ))

    splines = read_neutrino_cross_sections_info(path)
    log_energies = log10.(energy_grid)
    evaluated = map(
        (
            splines.νe_CC,
            splines.νμ_CC,
            splines.ντ_CC,
            splines.antiνe_CC,
            splines.antiνμ_CC,
            splines.antiντ_CC,
            splines.ν_NC,
            splines.antiν_NC,
        ),
    ) do spline
        values = evaluate_cubicspline.(
            Ref(spline),
            log_energies;
            extrapolate=extrapolate,
        )
        return clamp_nonnegative ? max.(values, 0.0) : values
    end
    return (
        energies=energy_grid,
        νe_CC=evaluated[1],
        νμ_CC=evaluated[2],
        ντ_CC=evaluated[3],
        antiνe_CC=evaluated[4],
        antiνμ_CC=evaluated[5],
        antiντ_CC=evaluated[6],
        ν_NC=evaluated[7],
        antiν_NC=evaluated[8],
        units=Symbol("(σ/E)_1e-42_m2_per_GeV_per_nucleon"),
        source=splines.source,
    )
end

"""
    plot_neutrino_cross_sections(table; channels=:CC)

Create a Makie plot of cross sections per nucleon. `channels` may be `:CC`,
`:NC`, or `:all`. The returned named tuple contains `figure` and `axis`.
"""
function plot_neutrino_cross_sections(
    table::NeutrinoCrossSectionTable;
    channels::Symbol=:CC,
)
    channels in (:CC, :NC, :all) ||
        throw(ArgumentError("channels must be :CC, :NC or :all"))
    selected = if channels === :CC
        (
            (:νe_CC, "νₑ CC"),
            (:νμ_CC, "νμ CC"),
            (:ντ_CC, "ντ CC"),
            (:antiνe_CC, "ν̄ₑ CC"),
            (:antiνμ_CC, "ν̄μ CC"),
            (:antiντ_CC, "ν̄τ CC"),
        )
    elseif channels === :NC
        ((:ν_NC, "ν NC"), (:antiν_NC, "ν̄ NC"))
    else
        (
            (:νe_CC, "νₑ CC"),
            (:νμ_CC, "νμ CC"),
            (:ντ_CC, "ντ CC"),
            (:antiνe_CC, "ν̄ₑ CC"),
            (:antiνμ_CC, "ν̄μ CC"),
            (:antiντ_CC, "ν̄τ CC"),
            (:ν_NC, "ν NC"),
            (:antiν_NC, "ν̄ NC"),
        )
    end

    figure = Figure()
    axis = Axis(
        figure[1, 1];
        xscale=log10,
        xlabel="E / GeV",
        ylabel="(σ/E) × 10⁻⁴² m²/GeV per nucleon",
    )
    for (field, label) in selected
        lines!(
            axis,
            table.energies,
            getproperty(table, field);
            label=label,
        )
    end
    axislegend(axis; position=:lt)
    return (; figure, axis)
end
