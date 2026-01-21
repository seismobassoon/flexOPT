

module Neurthino

    using LinearAlgebra
    using SparseArrays
    using Graphs
    #using StaticArrays
    #using Polynomials
    #using DocStringExtensions
    using LRUCache
    #using LightGraphs
    using AxisArrays

    import Base

    export oscprob, Pνν, Pνν, OscillationParameters, PMNSMatrix, Hamiltonian, MatterOscillationMatrices
    export masssquareddiff!, setΔm²!, cpphase!, setδ!, mixingangle!, setθ!
    export cpphases, mixingangles
    export Path
    export NeutrinoFlavour, Electron, Muon, Tau

    const N_A = 6.022e23 #[mol^-1]
    const G_F = 8.961877245622253e-38 #[eV*cm^3]

    # Julia 1.0 compatibility - normally obsolete
    #isnothing(::Any) = false
    #isnothing(::Nothing) = true

    #specific sources

    include("Neurthino/Oscillation.jl")
    include("Neurthino/Matter.jl")
    include("Neurthino/earthModels2Dor3D.jl")
    include("Neurthino/Earthmodels.jl")
    include("Neurthino/premFunctions.jl")
    include("Neurthino/usefulFunctionsToPlot.jl")
    include("Neurthino/NeurthinoRelated.jl")


    
    export N_A, G_F
    export oscprob




end # module
