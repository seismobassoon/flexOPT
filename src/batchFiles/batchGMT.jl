using GMT

# Nobuaki Fuji October 2025

#include("../src/batchRevise.jl")

if @isdefined DEFAULT_PLANET
    #
else
    include("../src/GeoPoints.jl")
end


function getTopoViaGMT(params::Dict)
    @unpack precision, region = params
    resolution=String(chopprefix(precision, "@earth_relief_"))
    if DEFAULT_PLANET[] === :Mars
        G = gmtread(remotegrid("mars", res = resolution))
    elseif DEFAULT_PLANET[] === :Earth
        # do nothing
    else
        @show "not yet coded for your planet"
    end

    topo = GMT.grdcut(precision, region=region)
    return @strdict(topo)
end