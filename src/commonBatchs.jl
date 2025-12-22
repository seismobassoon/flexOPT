module commonBatchs

    # batch files

  
    include("batchFiles/batchDrWatson.jl")
    include("batchFiles/batchImages.jl")
    include("batchFiles/batchUseful.jl")


    export myProduceOrLoad,lazyProduceOrLoad, @strdict
    export getColorPalette, regenerataionColorMap, color2Float
    export safeget
    export var2string,reinterpolateArrayMembers,expandVectors,string_as_varname
    export car2vec,carDropDim,carAddDim,vec2car
    export flatten,deep_flatten,is_all_less_than_or_equal,distance2_point_to_box
    export myInv

end