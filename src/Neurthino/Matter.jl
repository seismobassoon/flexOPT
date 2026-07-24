struct Path
    density::Vector{Float64}
    baseline::Vector{Float64}

    function Path(density::AbstractVector, baseline::AbstractVector)
        length(density) == length(baseline) ||
            throw(DimensionMismatch("density and baseline must have equal lengths"))
        isempty(density) && throw(ArgumentError("a Path cannot be empty"))
        all(isfinite, density) ||
            throw(ArgumentError("Path density values must be finite"))
        all(isfinite, baseline) ||
            throw(ArgumentError("Path baselines must be finite"))
        all(>=(0), density) ||
            throw(ArgumentError("Path density values must be non-negative"))
        all(>(0), baseline) ||
            throw(ArgumentError("Path baselines must be positive"))
        new(Float64.(density), Float64.(baseline))
    end
end

Path(density::Number, baseline::Number) = Path([density],[baseline])




Base.iterate(p::Path, state=1) = state > length(p.density) ? nothing : ( (p.density[state], p.baseline[state]),  state+1)

Base.length(p::Path) = length(p.density)

"""


Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `P`: Vacuum PMNS Matrix
- `H`: Vacuum Hamiltonian
- `density`: Matter density [g*cm^-3] 
- `energy`: Neutrino energy [GeV]
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function MatterOscillationMatrices(U, H, energy, density; zoa=0.5, anti=false)
    H_eff = convert(Array{ComplexF64}, U * Diagonal{Complex}(H) * adjoint(U))
    MatterOscillationMatrices(H_eff, energy, density; zoa=zoa, anti=anti)
end


"""


Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `H_eff`: Effective Matter Hamiltonian
- `density`: Matter density [g*cm^-3] 
- `energy`: Neutrino energy [GeV]
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function MatterOscillationMatrices(H_eff, energy, density; zoa=0.5, anti=false)
    A = sqrt(2) * G_F * N_A * density
    if anti
        H_eff[1,1] -= A * zoa * 2 * energy * 1e9
    else
        H_eff[1,1] += A * zoa * 2 * energy * 1e9
    end
    # Subtract Vz (= - A*Nn*E*1e9) for any sterile flavours
    if size(H_eff)[1] > 3
        for i in 4:size(H_eff)[1]
            if anti
                H_eff[i,i] -= A * (1 - zoa) * energy * 1e9
            else
                H_eff[i,i] += A * (1 - zoa) * energy * 1e9
            end
        end
    end
    tmp = eigen(H_eff)
    return tmp.vectors, tmp.values
end

"""


Create modified oscillation parameters for neutrino propagation through matter

# Arguments
- `osc_vacuum::OscillationParameters`: Oscillation parameters in vacuum
- `energy`: Neutrino energy [GeV]
- `density`: Matter density in g*cm^-3 
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino

"""
function MatterOscillationMatrices(osc_vacuum::OscillationParameters, energy, density; zoa=0.5, anti=false)
    H_vacuum = Diagonal(Hamiltonian(osc_vacuum)) 
    U_vacuum = PMNSMatrix(osc_vacuum; anti=anti)
    H_eff = convert(Array{ComplexF64}, U_vacuum * Diagonal{ComplexF64}(H_vacuum) * adjoint(U_vacuum))
    return MatterOscillationMatrices(H_eff, energy, density; zoa=zoa, anti=anti)
end

"""


# Arguments
- `U`: Vacuum PMNS Matrix
- `H`: Vacuum Hamiltonian
- `energy`: Neutrino energy [GeV]
- `path::Vector{Path}`: Neutrino path
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function oscprob(
    U,
    H,
    energy::AbstractVector{T},
    path::AbstractVector{<:Path};
    zoa=0.5,
    anti=false,
) where {T <: Real}
    energies = Float64.(collect(energy))
    isempty(path) && throw(ArgumentError("path collection cannot be empty"))
    if anti
        H_eff = conj.(U) * Diagonal{ComplexF64}(H) * adjoint(conj.(U))
    else
        H_eff = U * Diagonal{ComplexF64}(H) * adjoint(U)
    end
    A = zeros(ComplexF64, length(energies), length(path), size(U)...)
    cache_size = length(energies) * sum(map(x->length(x.density), path))
    lru = LRU{Tuple{Float64, Float64},
              Tuple{Array{ComplexF64,2}, Vector{ComplexF64}}}(maxsize=cache_size)
    for k in eachindex(energies)
        @inbounds E = energies[k]
        for (l, p) in enumerate(path)
            tmp = Matrix{ComplexF64}(1I, size(U))
            for (m,b) in enumerate(p.baseline)
                @inbounds ρ = p.density[m]

                U_mat, H_mat = get!(lru, (E, ρ)) do
                    MatterOscillationMatrices(copy(H_eff), E, ρ; zoa=zoa, anti=anti)
                end

                # Profiles are ordered source → detector. For column-state
                # evolution, each later segment multiplies on the left.
                segment_evolution = _oscprobampl(U_mat, H_mat, E, b)
                tmp = segment_evolution * tmp
            end
            # Store the public probability axes as (initial, final).
            @inbounds A[k, l, :, :] = transpose(tmp)
        end
    end
    P = map(x -> abs.(x) .^ 2, A)
    flavrange = _make_flavour_range(first(size(U)))
    AxisArray(P; Energy=energies, Path=collect(path), InitFlav=flavrange, FinalFlav=flavrange)
end

#const oscprob(U, H, energy::T, path::Vector{Path}; zoa=0.5, anti=false) where {T <: Real} = oscprob(U, H, [energy], path; zoa=zoa, anti=anti)

#const oscprob(U, H, energy, path::Path; zoa=0.5, anti=false) = oscprob(U, H, energy, [path]; zoa=zoa, anti=anti)

"""


# Arguments
- `osc_vacuum::OscillationParameters`: Vacuum oscillation parameters
- `energy`: Neutrino energy [GeV]
- `path`: Neutrino path
- `zoa`: Proton nucleon ratio (Z/A)
- `anti`: Is anti neutrino
"""
function oscprob(
    osc_vacuum::OscillationParameters,
    energy,
    path::Union{Path,AbstractVector{<:Path}};
    zoa=0.5,
    anti=false,
)
    # TODO: attach U_vac and H_vac to the oscillation parameters, so that it's
    # only calculated once and invalidated when any of the oscillation parameters
    # are changed
    U_vac = PMNSMatrix(osc_vacuum; anti=anti)
    H_vac = Hamiltonian(osc_vacuum)
    paths = path isa Path ? [path] : path
    oscprob(U_vac, H_vac, energy, paths; zoa=zoa, anti=anti)
end
