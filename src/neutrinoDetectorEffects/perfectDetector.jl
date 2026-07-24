"""
    PerfectDetector

Ideal detector exposure used before efficiency, resolution and reconstruction
effects are applied. Masses are in kg and live time is in seconds.
"""
Base.@kwdef struct PerfectDetector
    mass_kg::Float64 = 1.0e9
    livetime_s::Float64 = 365.25 * 24.0 * 60.0 * 60.0
    nucleon_mass_kg::Float64 = 1.67262e-27
end

const PerfectDetectorEventTable = NamedTuple{
    (:energies, :cosθgrid, :energy_bin_widths, :cosθ_bin_widths,
     :solid_angle_bin_widths, :detector, :neutrino, :antineutrino,
     :combined, :units),
}

function _validated_centers(centers::AbstractVector{<:Real}, name::AbstractString)
    values = Float64.(collect(centers))
    length(values) >= 2 ||
        throw(ArgumentError("$name needs at least two bin centers"))
    all(isfinite, values) ||
        throw(ArgumentError("$name must contain only finite values"))
    all(diff(values) .> 0) ||
        throw(ArgumentError("$name must be strictly increasing"))
    return values
end

"""
    bin_edges_from_centers(centers; scale=:linear, bounds=nothing)

Construct bin edges from ordered centers. For logarithmic coordinates, interior
edges are geometric means and end edges continue the neighboring log spacing.
Optional physical `bounds` clamp the two outer edges.
"""
function bin_edges_from_centers(
    centers::AbstractVector{<:Real};
    scale::Symbol=:linear,
    bounds::Union{Nothing,Tuple{<:Real,<:Real}}=nothing,
)
    values = _validated_centers(centers, "bin centers")
    scale in (:linear, :log) ||
        throw(ArgumentError("scale must be :linear or :log"))

    edges = if scale === :linear
        interior = (values[1:end-1] .+ values[2:end]) ./ 2
        [
            values[1] - (interior[1] - values[1]);
            interior;
            values[end] + (values[end] - interior[end])
        ]
    else
        all(>(0), values) ||
            throw(ArgumentError("logarithmic bin centers must be positive"))
        log_values = log.(values)
        log_interior = (log_values[1:end-1] .+ log_values[2:end]) ./ 2
        exp.([
            log_values[1] - (log_interior[1] - log_values[1]);
            log_interior;
            log_values[end] + (log_values[end] - log_interior[end])
        ])
    end

    if !isnothing(bounds)
        lower, upper = Float64.(bounds)
        lower < upper ||
            throw(ArgumentError("bin-edge bounds must be increasing"))
        edges[1] = max(edges[1], lower)
        edges[end] = min(edges[end], upper)
        edges[1] < edges[2] && edges[end-1] < edges[end] ||
            throw(ArgumentError("bounds exclude one or more bin centers"))
    end
    return edges
end

bin_widths_from_centers(centers; kwargs...) =
    diff(bin_edges_from_centers(centers; kwargs...))

function _raw_probability_array(probabilities)
    array = parent(probabilities)
    ndims(array) == 4 ||
        throw(DimensionMismatch(
            "oscillation probabilities must have dimensions " *
            "(energy, cosθ/path, initial flavor, final flavor)",
        ))
    size(array, 3) >= 3 && size(array, 4) >= 3 ||
        throw(DimensionMismatch(
            "oscillation probabilities need at least three flavors",
        ))
    return array
end

"""
    oscillation_channel_table(neutrino, antineutrino)

Expose oscillation channels as nested named tuples without transposing or
copying the underlying `(energy, cosθ, initial, final)` arrays.
"""
function oscillation_channel_table(neutrino, antineutrino)
    ν = _raw_probability_array(neutrino)
    antiν = _raw_probability_array(antineutrino)
    size(ν) == size(antiν) ||
        throw(DimensionMismatch(
            "neutrino and antineutrino probability arrays must have equal size",
        ))

    channels(array) = (
        e_to_e=@view(array[:, :, 1, 1]),
        e_to_μ=@view(array[:, :, 1, 2]),
        e_to_τ=@view(array[:, :, 1, 3]),
        μ_to_e=@view(array[:, :, 2, 1]),
        μ_to_μ=@view(array[:, :, 2, 2]),
        μ_to_τ=@view(array[:, :, 2, 3]),
        τ_to_e=@view(array[:, :, 3, 1]),
        τ_to_μ=@view(array[:, :, 3, 2]),
        τ_to_τ=@view(array[:, :, 3, 3]),
    )
    return (
        neutrino=channels(ν),
        antineutrino=channels(antiν),
        size=size(ν),
    )
end

function _validate_detector(detector::PerfectDetector)
    detector.mass_kg > 0 ||
        throw(ArgumentError("detector mass must be positive"))
    detector.livetime_s > 0 ||
        throw(ArgumentError("detector live time must be positive"))
    detector.nucleon_mass_kg > 0 ||
        throw(ArgumentError("nucleon mass must be positive"))
    return detector
end

function _validate_interaction_inputs(flux, cross_sections, channels)
    energies = _validated_centers(flux.energies, "flux energies")
    cosines = _validated_centers(flux.cosθgrid, "flux cosθ grid")
    cross_sections.energies == energies ||
        throw(DimensionMismatch(
            "cross sections and flux must use the same energy grid",
        ))
    channels.size[1:2] == (length(energies), length(cosines)) ||
        throw(DimensionMismatch(
            "probability dimensions $(channels.size[1:2]) do not match " *
            "the flux grid $((length(energies), length(cosines)))",
        ))
    expected = (length(energies), length(cosines))
    for field in (:νe, :νμ, :antiνe, :antiνμ)
        size(getproperty(flux, field)) == expected ||
            throw(DimensionMismatch(
                "flux.$field must have size $expected",
            ))
    end
    return (; energies, cosines)
end

function _flavor_events(
    normalization,
    flux_e,
    flux_μ,
    probability_e,
    probability_μ,
    cross_section_over_energy,
    energies,
)
    differential_flux = (
        flux_e .* probability_e .+
        flux_μ .* probability_μ
    )
    # JSON values are (σ/E) × 10⁻⁴² m²/GeV per nucleon.
    cross_section_m² = cross_section_over_energy .* energies .* 1.0e-42
    return normalization .* reshape(cross_section_m², :, 1) .* differential_flux
end

"""
    interacting_events(flux, cross_sections, neutrino_probabilities,
        antineutrino_probabilities; detector=PerfectDetector())

Compute ideal binned charged-current event counts. Flux matrices and
probabilities must already use `(energy, cosθ)` ordering. The calculation uses

`N = (M/mₙ) T σ(E) Φ(E,cosθ) P ΔE 2πΔcosθ`.

No empirical factor, transpose, or extra energy power is applied.
"""
function interacting_events(
    flux,
    cross_sections,
    neutrino_probabilities,
    antineutrino_probabilities;
    detector::PerfectDetector=PerfectDetector(),
)
    _validate_detector(detector)
    channels = oscillation_channel_table(
        neutrino_probabilities,
        antineutrino_probabilities,
    )
    (; energies, cosines) = _validate_interaction_inputs(
        flux,
        cross_sections,
        channels,
    )

    energy_widths = bin_widths_from_centers(energies; scale=:log)
    cosine_widths = bin_widths_from_centers(
        cosines;
        scale=:linear,
        # The propagation grid includes its requested angular endpoints
        # (for example -1 and 0 for the upgoing hemisphere).
        bounds=(first(cosines), last(cosines)),
    )
    solid_angle_widths = 2π .* cosine_widths
    exposure_nucleon_seconds =
        detector.livetime_s * detector.mass_kg / detector.nucleon_mass_kg
    normalization =
        exposure_nucleon_seconds .*
        reshape(energy_widths, :, 1) .*
        reshape(solid_angle_widths, 1, :)

    ν = channels.neutrino
    antiν = channels.antineutrino
    neutrino_events = (
        νe=_flavor_events(
            normalization,
            flux.νe,
            flux.νμ,
            ν.e_to_e,
            ν.μ_to_e,
            cross_sections.νe_CC,
            energies,
        ),
        νμ=_flavor_events(
            normalization,
            flux.νe,
            flux.νμ,
            ν.e_to_μ,
            ν.μ_to_μ,
            cross_sections.νμ_CC,
            energies,
        ),
        ντ=_flavor_events(
            normalization,
            flux.νe,
            flux.νμ,
            ν.e_to_τ,
            ν.μ_to_τ,
            cross_sections.ντ_CC,
            energies,
        ),
    )
    antineutrino_events = (
        νe=_flavor_events(
            normalization,
            flux.antiνe,
            flux.antiνμ,
            antiν.e_to_e,
            antiν.μ_to_e,
            cross_sections.antiνe_CC,
            energies,
        ),
        νμ=_flavor_events(
            normalization,
            flux.antiνe,
            flux.antiνμ,
            antiν.e_to_μ,
            antiν.μ_to_μ,
            cross_sections.antiνμ_CC,
            energies,
        ),
        ντ=_flavor_events(
            normalization,
            flux.antiνe,
            flux.antiνμ,
            antiν.e_to_τ,
            antiν.μ_to_τ,
            cross_sections.antiντ_CC,
            energies,
        ),
    )
    combined = (
        νe=neutrino_events.νe .+ antineutrino_events.νe,
        νμ=neutrino_events.νμ .+ antineutrino_events.νμ,
        ντ=neutrino_events.ντ .+ antineutrino_events.ντ,
    )
    combined = merge(
        combined,
        (total=combined.νe .+ combined.νμ .+ combined.ντ,),
    )

    return (
        energies=energies,
        cosθgrid=cosines,
        energy_bin_widths=energy_widths,
        cosθ_bin_widths=cosine_widths,
        solid_angle_bin_widths=solid_angle_widths,
        detector=detector,
        neutrino=neutrino_events,
        antineutrino=antineutrino_events,
        combined=combined,
        units=:events_per_bin,
    )
end

"""
    plot_interacting_events(events; flavor=:νe, component=:combined,
        logscale=true)

Create a Makie heatmap of ideal interacting-event counts. `component` may be
`:neutrino`, `:antineutrino`, or `:combined`; combined also supports
`flavor=:total`.
"""
function plot_interacting_events(
    events::PerfectDetectorEventTable;
    flavor::Symbol=:νe,
    component::Symbol=:combined,
    logscale::Bool=true,
)
    component in (:neutrino, :antineutrino, :combined) ||
        throw(ArgumentError(
            "component must be :neutrino, :antineutrino or :combined",
        ))
    selected = getproperty(events, component)
    hasproperty(selected, flavor) ||
        throw(ArgumentError(
            "flavor $flavor is unavailable for component $component",
        ))
    counts = getproperty(selected, flavor)
    plotted = if logscale
        log10.(max.(counts, floatmin(Float64)))
    else
        counts
    end

    figure = Figure()
    axis = Axis(
        figure[1, 1];
        xlabel="log₁₀(E / GeV)",
        ylabel="cosθ",
        title="$(component) $(flavor)",
    )
    heatmap_plot = heatmap!(
        axis,
        log10.(events.energies),
        events.cosθgrid,
        plotted,
    )
    colorbar = Colorbar(
        figure[1, 2],
        heatmap_plot;
        label=logscale ? "log₁₀(events/bin)" : "events/bin",
    )
    return (; figure, axis, plot=heatmap_plot, colorbar)
end
