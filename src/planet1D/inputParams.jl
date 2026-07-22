using YAML, StructTypes

"""Runtime parameters for `planet1D`.

The module is self-contained: loading it no longer requires a `Main.ParamFile`.
Call [`configure_input!`](@ref) to change frequency or time-window parameters.
"""
input = Input()
metaInfo = MetaInfo()
dsm1Dconfig = planet1Dconfig()

StructTypes.StructType(::Type{MetaInfo}) = StructTypes.Mutable()
StructTypes.StructType(::Type{planet1Dconfig}) = StructTypes.Mutable()

function _load_dsm_config()
    metadata_path = normpath(joinpath(@__DIR__, "../../dataInput/config.yaml"))
    isfile(metadata_path) || return metaInfo, dsm1Dconfig

    metadata = YAML.load_file(metadata_path; dicttype=Dict{Symbol,Any})
    metadata_values = metadata[:classicDSM][1]
    loaded_meta = StructTypes.constructfrom(MetaInfo, metadata_values)

    config_path = normpath(joinpath(@__DIR__, "..", loaded_meta.planet1DconfigFile))
    isfile(config_path) || return loaded_meta, dsm1Dconfig
    config = YAML.load_file(config_path; dicttype=Dict{Symbol,Any})
    return loaded_meta, StructTypes.constructfrom(planet1Dconfig, config[:global][1])
end

metaInfo, dsm1Dconfig = _load_dsm_config()

"""
    configure_input!(; averagedPlanetRadius=0, timeWindow=nothing,
                       timeWindowMinimum=nothing, maxFrequencyMin=2,
                       minFrequencyMax=0, GUIoption=false)

Configure the global `planet1D.input` without relying on a CSV parameter file.
Frequencies are expressed in Hz and time windows in seconds. If only
`timeWindowMinimum` is supplied, the window is rounded up to `0.1 * 2^n`.
"""
function configure_input!(; averagedPlanetRadius::Real=0.0,
                          timeWindow::Union{Nothing,Real}=nothing,
                          timeWindowMinimum::Union{Nothing,Real}=nothing,
                          maxFrequencyMin::Real=2.0,
                          minFrequencyMax::Real=0.0,
                          GUIoption::Bool=false)
    maxFrequencyMin > 0 || throw(ArgumentError("maxFrequencyMin must be positive"))
    minFrequencyMax >= 0 || throw(ArgumentError("minFrequencyMax must be non-negative"))

    selected_window = if timeWindow !== nothing
        Float64(timeWindow)
    elseif timeWindowMinimum !== nothing
        minimum_window = Float64(timeWindowMinimum)
        minimum_window > 0 || throw(ArgumentError("timeWindowMinimum must be positive"))
        0.1 * 2.0^ceil(Int, log2(minimum_window / 0.1))
    else
        3276.8
    end
    selected_window > 0 || throw(ArgumentError("timeWindow must be positive"))

    i_end = 2^floor(Int, log2(Float64(maxFrequencyMin) * selected_window))
    i_start = minFrequencyMax > 0 ?
        2^floor(Int, log2(Float64(minFrequencyMax) * selected_window)) : 1
    i_start <= i_end || throw(ArgumentError("minFrequencyMax must not exceed maxFrequencyMin"))

    input.averagedPlanetRadius = Float64(averagedPlanetRadius)
    input.timeWindow = selected_window
    input.ωᵢ = -log(dsm1Dconfig.omegaiTlen) / selected_window
    dω = 2π / selected_window
    input.ωᵣ = collect(range(dω * i_start; step=dω, length=i_end - i_start + 1))
    input.GUIoption = GUIoption
    return input
end

configure_input!()
