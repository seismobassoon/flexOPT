# some useful structures and functions to deal with planetary 3D coordinates

module GeoPoints
    # be careful, this module needs planet1D

    using Geodesy, StaticArrays,LinearAlgebra,GMT, UnPack
    using Interpolations
    import Base: +,-,/,*


    using ..planet1D
    
    # batch files

  
    include("batchFiles/batchDrWatson.jl")

    # specific functions
    include("GeoPoints/GeoPoints.jl")
    include("GeoPoints/getSeismicParamTopo.jl")
    

    export DEFAULT_PLANET,DEFAULT_ELLIPSOID, planet_ellipsoid, set_default_planet!, localCoord2D, localCoord3D
    export GeoPoint, +, -, /, *
    #export effectiveRadius, makeLocalCoordinateUnitVectors, makeLocalCoordinateUnitVectors, makeLocalCoordinateUnitVectors
    export constructLocalBox
    export getParamsAndTopo, interpolate_params, GMTprecision

end