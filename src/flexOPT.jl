module flexOPT

    using Symbolics,UnPack,Tullio
    using KernelAbstractions
    using ..commonBatchs

    # GPU backend (batchFiles/batchGPU.jl should be called in Main)
    backend = Main.backend
    makeGPUarray = Main.makeGPUarray

    # batch files

    include("batchFiles/batchSymbolics.jl")
    export x,y,z,t
    export ∂x,∂y,∂z,∂t
    export ∂x²,∂y²,∂z²,∂t² 
    export Num2Float64,usefulPartials,myCoeff,mySolvefor,mySimplify,integrateTaylorPolynomials


    # Semi symbolic OPT motors
    include("motorsOPT/WYYKK.jl")
    export getIntegralWYYKK
    include("motorsOPT/IntegrateBsplineAndPolynomials.jl")
    export BsplineTimesPolynomialsIntegrated
    include("motorsOPT/others.jl")
    #export findCartesianDependency, makeMixPartials,varMmaker,PDECoefFinder,timeDimensionString
    export timeDimensionString
    include("motorsOPT/OPTEngines.jl")
    export OPTobj
    
    # famous equations
    include("motorsOPT/famousEquations.jl")
    include("motorsOPT/famousSourceFunctions.jl")
    export famousEquations, Ricker

end
