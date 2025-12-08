# some useful structures and functions to deal with planetary 3D coordinates

module GeoPoints


    using Geodesy, StaticArrays,LinearAlgebra
    using Interpolations
    import Base: +,-,/,*


    # specific functions

    include("GeoPoints/GeoPoints.jl")
    include("GeoPoints/getSeismicParamTopo.jl")
    include("batchFiles/batchDrWatson.jl")
    

    # batch files

    include("batchFiles/batchGMT.jl")
    include("batchFiles/batchDrWatson.jl")

end