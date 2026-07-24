#= these functions here are the beautiful fruits of Estelle Salomé's 2-month internship in 2025 =#
# and Isabel & Yael 2026

"""
    set_oscillation_parameters(ordering=:normal)

Initializes a 3-flavor neutrino oscillation object using the `Neurthino` framework
and the latest global best-fit values from NuFIT v6.1.
Accepts either `:normal` or `:inverted` for the `ordering` keyword.

### Examples:
set_oscillation_parameters()
set_oscillation_parameters(:inverted)
"""
function set_oscillation_parameters(ordering::Symbol=:normal)
 
    # initialize a 3-flavor oscillation object
    osc = OscillationParameters(3)
    # NuFIT v6.1 normal ordering with SK atmospheric data
    if ordering === :normal
        # mixing angles (radians)
        setθ!(osc, 1=>2, 0.589)   # θ12 ≈ 33.76°
        setθ!(osc, 1=>3, 0.150)   # θ13 ≈ 8.62°
        setθ!(osc, 2=>3, 0.756)   # θ23 ≈ 43.29° 
        # CP-violating phase (radians)
        setδ!(osc, 1=>3, 3.700)   # δCP ≈ 212°
        # Mass-squared splittings (eV²)
        setΔm²!(osc, 2=>1, 7.537e-5)   
        setΔm²!(osc, 3=>2, 2.511e-3)  
        
    elseif ordering === :inverted
        # Mixing angles (radians)
        setθ!(osc, 1=>2, 0.589)   # θ12 ≈ 33.76°
        setθ!(osc, 1=>3, 0.151)   # θ13 ≈ 8.65°
        setθ!(osc, 2=>3, 0.836)   # θ23 ≈ 47.90° 
        # CP-violating phase (radians)
        setδ!(osc, 1=>3, 4.782)   # δCP ≈ 274°
        # Mass-squared splittings (eV²)
        setΔm²!(osc, 2=>1,  7.537e-5)  
        setΔm²!(osc, 3=>2, -2.483e-3) 
        
    else
        error("Unknown mass ordering: :\$ordering. Use :normal or :inverted.")

    end
    
    return osc
end

"""
    set_oscillation_parameters(θ12, θ13, θ23, δCP, m12, m23)

Initializes a 3-flavor neutrino oscillation object using the `Neurthino` framework
and user-defined oscillation parameters.

### Note on units:
* Mixing angles (`θ12`, `θ13`, `θ23`) & CP-violating phase (`δCP`): must be provided in degrees. 
* Mass-squared splittings (`m12`, `m23`): must be provided in eV^2.

### Example:
set_oscillation_parameters(34.0, 9.0, 43.0, 212.0, 7.4e-5, 2.5e-3)
"""
function set_oscillation_parameters(θ12, θ13, θ23, δCP, m12, m23)
 
    # initialize a 3-flavor oscillation object
    osc = OscillationParameters(3)

    # mixing angles (radians)
    setθ!(osc, 1=>2, θ12*pi/180)   
    setθ!(osc, 1=>3, θ13*pi/180)   
    setθ!(osc, 2=>3, θ23*pi/180)  
    # CP-violating phase (radians)
    setδ!(osc, 1=>3, δCP)  
    # Mass-squared splittings (eV²)
    setΔm²!(osc, 2=>1, m12)   
    setΔm²!(osc, 3=>2, m23)  

    return osc
end

"""
    set_oscillation_parameters(; ordering=:normal, kwargs...)

Configure a 3-flavor neutrino oscillation object using default oscillation parameters
and doing some manual modifications.

### Arguments:
* `ordering::Symbol`: The baseline mass ordering. Options are `:normal` (default) or `:inverted`.

### Keyword arguments (optional modifications):
* `θ12`, `θ13`, `θ23`: in degrees.
* `δCP`: CP-violating phase in degrees.
* `m12`, `m23`: Mass-squared splittings in eV^2.

### Example:
# Use inverted ordering as the baseline, but modify just the CP-violating phase:
set_oscillation_parameters(ordering=:inverted, δCP=1.5)
"""
function set_oscillation_parameters(; ordering::Symbol=:normal,
    θ12 = nothing,
    θ13 = nothing,
    θ23 = nothing,
    δCP = nothing,
    m12 = nothing,
    m23 = nothing)

    # initialize a 3-flavor oscillation object 
    # and assign it the oscillation parameters of the chosen ordering
    osc = set_oscillation_parameters(ordering)

    # modify the parameters of interest
    if !isnothing(θ12) setθ!(osc, 1=>2, θ12*pi/180) end
    if !isnothing(θ13) setθ!(osc, 1=>3, θ13*pi/180) end 
    if !isnothing(θ23) setθ!(osc, 2=>3, θ23*pi/180) end 
    if !isnothing(δCP) setδ!(osc, 1=>3, δCP) end 
    if !isnothing(m12) setΔm²!(osc, 2=>1, m12) end
    if !isnothing(m23) setΔm²!(osc, 3=>2, m23) end 

    return osc
end



function solveQuadraticEquation(a,b,c)
    Δ = b^2 - 4*a*c

    if Δ>0
        x1 = ((-b - sqrt(Δ))/(2*a))
        x2 = ((-b + sqrt(Δ))/(2*a))
    else
        x1 = -b/(2 *a)
        x2 = x1
    end
    return x1, x2

end

function posOrNeg(cos_θ, sign = :positive)
    #to choose if we want a positive or a negative angle
    if sign == :positive
        θ = acos.(cos_θ)
    else
        θ = .- acos.(cos_θ)
    end
    return θ
end


function roundExt(x,step)
    #to round the values of H and U
    if x isa Complex
        xreal = round(real(x)/step)*step
        ximag = round(imag(x)/step)*step
        return complex(xreal, ximag)
    else
        return round(x/step)*step
    end
end

"""
    pathsFromElectronDensity(path_sampling; reference_z_over_a=0.5)

Convert source-to-detector profiles produced by `GeoPoints.creationPaths` from
the electron-density proxy `nₑ = ρ * (Z/A)` into `Path` objects accepted by
`Pνν`.

The matter solver uses `density * zoa`. Consequently, each path stores the
equivalent density

    effective_density = nₑ / reference_z_over_a

and `Pνν` must be called with the same `zoa=reference_z_over_a`. With the
conventional reference value 0.5, this is the explicit conversion
`effective_density = 2nₑ`.

The returned NamedTuple exposes the scale and convention so this physical
assumption remains visible to notebook users.
"""
function pathsFromElectronDensity(
    path_sampling;
    reference_z_over_a::Real=0.5,
)
    0 < reference_z_over_a <= 1 ||
        throw(ArgumentError("reference_z_over_a must lie in (0, 1]"))
    hasproperty(path_sampling, :profiles) ||
        throw(ArgumentError("path_sampling must contain GeoPoints path profiles"))

    profiles = path_sampling.profiles
    isempty(profiles) && throw(ArgumentError("path_sampling.profiles cannot be empty"))
    effective_density_scale = inv(Float64(reference_z_over_a))
    paths = map(profiles) do profile
        hasproperty(profile, :values) && hasproperty(profile, :baselines) ||
            throw(ArgumentError("each profile must contain values and baselines"))
        Path(profile.values .* effective_density_scale, profile.baselines)
    end

    return (;
        paths,
        reference_z_over_a=Float64(reference_z_over_a),
        effective_density_scale,
        cosθgrid=Float64.(collect(path_sampling.cosθgrid)),
        convention=:electron_density_proxy_to_effective_mass_density,
    )
end
