
# For now, the default neutrino flux is the same as in EarthProbe
const DEFAULT_NUFLUX_PARAMS = Ref(Dict{String, Any}(
    "model"                  => :Honda,
    "location_hf"            => "gran_sasso",      
    "season_hf"              => "all_year",
    "angles_separation_hf"   => "averaged_ϕ",       
    "mountain_overburden_hf" => "without",
    "solar_activity_hf"      => "minimum"   
))

const DEFAULT_NUFLUX = Ref{Any}(nothing)

function __init__()
    DEFAULT_NUFLUX[] = produce_neutrino_flux(DEFAULT_NUFLUX_PARAMS[])
end

"""
    set_default_neutrino_flux!(nuflux_params::Dict{String, Any})

Updates the default configuration parameters and recalculates the global 
default neutrino flux data using the new settings.
"""
function set_default_neutrino_flux!(nuflux_params::Dict{String, Any})
    DEFAULT_NUFLUX_PARAMS[] = nuflux_params
    DEFAULT_NUFLUX[]        = produce_neutrino_flux(nuflux_params)
end

const LOCATION_MAPPING = Dict{Symbol, String}(
    :kamioka     => "kam",
    :gran_sasso  => "grn",
    :sudbury     => "sno",
    :frejus      => "frj",
    :ino         => "ino",
    :south_pole  => "spl", 
    :pythasalmi  => "pyh",
    :homestake   => "hms",
    :juno        => "juno"
)

const SEASON_MAPPING = Dict{Symbol, String}(
    :all_year  => "ally",
    :march_may => "0305",
    :june_aug  => "0608",
    :sept_nov  => "0911",
    :dec_feb   => "1202"
)

const ANGLES_MAPPING = Dict{Symbol, String}(
    :variable_ϕ   => "20-12", 
    :averaged_ϕ   => "20-01",
    :averaged_ϕθ  => "01-01" 
)

const MOUNTAIN_MAPPING = Dict{Symbol, String}(
    :with    => "-mtn",  
    :without => "" 
)

const SOLAR_MAPPING = Dict{Symbol, String}(
    :minimum => "solmin",
    :maximum => "solmax"
)
