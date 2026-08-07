module flexOPT

    using Symbolics,UnPack,Tullio
    using Dates
    using LinearAlgebra
    using KernelAbstractions
    using ..commonBatchs
    using Adapt

    # Use the backend selected by batchFiles/batchGPU.jl when the caller has
    # loaded it.  A module import must nevertheless remain valid on its own:
    # recipe construction then falls back to the KernelAbstractions CPU.
    const backend = isdefined(Main, :backend) ?
        getfield(Main, :backend) : KernelAbstractions.CPU()
    const makeGPUarray = isdefined(Main, :makeGPUarray) ?
        getfield(Main, :makeGPUarray) :
        ((selected_backend, array) -> Adapt.adapt(selected_backend, array))


    # wrapping functions

    include("../src/motorsOPT/makeOPTrecette.jl")
    export TaylorOptions,makeOPTsemiSymbolic
    
    # deriving symbolic and semi-symbolic OPT operators
    include("../src/compactSymbolicFunctions/compactFunctionsArray.jl")
    include("../src/compactSymbolicFunctions/BsplineHelpers.jl")
    include("../src/CompactSymbolicFunctions/TaylorExpansionHelpers.jl")
    include("../src/CompactSymbolicFunctions/integralWYYKK.jl")


    # famous equations etc.

    include("../src/motorsOPT/others.jl")
    export timeDimensionString
    include("../src/motorsOPT/famousEquations.jl")
    include("../src/motorsOPT/famousSourceFunctions.jl")
    include("../src/motorsOPT/famousBoundaryConditions.jl")
    export CerjanBoundarySpec, BoundaryConditionSet, boundary_geometry
    export cerjan_padding, famousBoundaryConditions

    # semi-symbolics operators to fully numerical operators
    #include("fullyNumericalOPT/makeCostFunctions.jl")
    include("fullyNumericalOPT/makeCostFunctions.jl")
    include("fullyNumericalOPT/others.jl")
    export numericalOperatorConstruction,getModelPoints

    include("../src/numSolvers/linearSolution.jl")
    include("../src/numSolvers/timeMarchingSchemes.jl")
    include("../src/numSolvers/diffTools.jl")
    export prepareNumericalOperators,timeMarchingSchemePrepared,prepareLinearSystem,evaluateLinearSystem!,evaluateLinearSystem,timeMarchingSchemeLinear,propagateLinearSystem,overlapBoundaryLinearSystem
    export prepareLumpedFuturePredictor,applyLumpedFuturePredictor!


    export quasiNumericalOperatorConstruction,constructingNumericalDiscretisedEquations

    export famousEquations, Ricker

end
