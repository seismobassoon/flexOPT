# some useful structures and functions to deal with planetary 3D coordinates

module GeoPoints
    # be careful, this module needs planet1D


    using ..planet1D
    using ..commonBatchs
    using Geodesy, StaticArrays,LinearAlgebra,GMT, UnPack
    using Interpolations
    using FileIO,CairoMakie
    using Colors, Images
    import Base: +,-,/,*


    



    # specific functions: former GeoPoints
    include("GeoPoints/GeoPoints.jl")
    include("GeoPoints/getSeismicParamTopo.jl")
    
    # specific functions: former imageReder
    include("imageReader/imageReader.jl")
    include("imageReader/modelMaker.jl")


    # functions GeoPoints
    export DEFAULT_PLANET,DEFAULT_ELLIPSOID, planet_ellipsoid, set_default_planet!, localCoord2D, localCoord3D
    export GeoPoint, +, -, /, *
    #export effectiveRadius, makeLocalCoordinateUnitVectors, makeLocalCoordinateUnitVectors, makeLocalCoordinateUnitVectors
    export constructLocalBox
    export getParamsAndTopo,getParamsWithoutTopo #, interpolate_params, GMTprecision
    export makeAdHocSeismicModel,initiateSeismicModel # soon we need to make a structure of seismicModel (1D-3D)
    


    # functions imageReader
    export defineModel
    export read2DimageModel
    export adjustArray
end