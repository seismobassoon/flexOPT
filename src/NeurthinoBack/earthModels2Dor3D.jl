# Estelle Salomé (stagiaire ENSG 2ème année) and Nobuaki Fuji (IPGP/UPC/IUF)
# extends PREM.jl in Neurthino.jl package for a half circle disk


function earthpath(zenith::T, zposition; samples=100, discrete_densities=nothing,EARTH_RADIUS=6371.0) where {T <: Number}
"""


# Arguments
- `zenith::Float64`: Zenith angle of the path with respect to the detector frame [rad]
- `zposition::Float64` Detector depth (positive value) [km] 
- `samples` The number of steps with equal distance
"""
    trklen = tracklength(zenith, zposition)
    x = Array(range(0.0; stop=trklen, length=samples))
    sections = (x[2:end] - x[1:end-1])
    total_pathlen = 0.5 * (x[2:end] + x[1:end-1])
    zprime = EARTH_RADIUS - zposition
    radii = map(x -> sqrt(zprime^2 + x^2 + 2*zprime*x*cos(zenith)), total_pathlen)
    @show densities = PREM.(radii)
    if !isnothing(discrete_densities)
        idx = map(d->searchsortedfirst(discrete_densities,d), densities)
        densities = map(i->discrete_densities[i], idx)
    end
    Path(densities, sections)
end
