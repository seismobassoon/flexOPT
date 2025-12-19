module flexOPT

    using Symbolics,UnPack,Tullio

    # batch files

    include("batchFiles/batchSymbolics.jl")
    export x,y,z,t
    export ∂x,∂y,∂z,∂t
    export ∂x²,∂y²,∂z²,∂t² 
    export Num2Float64,usefulPartials,myCoeff,mySolvefor,mySimplify,integrateTaylorPolynomials


    # Semi symbolic OPT motors
    include("motorsOPT/WYYKK.jl")
    export getIntegralWYYKK
    include("motorsOPT/others.jl")
    export findCartesianDependency, makeMixPartials,varMmaker,PDECoefFinder
    

    
    # famous equations
    include("motorsOPT/famousEquations.jl")
    include("motorsOPT/famousSourceFunctions.jl")
    export famousEquations, Ricker

end
