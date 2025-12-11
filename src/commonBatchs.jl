module commonBatchs

    # batch files

  
    include("batchFiles/batchDrWatson.jl")
    include("batchFiles/batchImages.jl")


    export myProduceOrLoad, lazyProduceOrLoad, @strdict
    export getColorPalette, regenerataionColorMap, color2Float

end