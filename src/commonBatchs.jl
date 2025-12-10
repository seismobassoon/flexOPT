module commonBatchs

    # batch files

  
    include("batchFiles/batchDrWatson.jl")
    include("batchFiles/batchImages.jl")


    export myProduceOrLoad, lazyProduceOrLoad

end