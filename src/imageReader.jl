module imageReader
    using FileIO,CairoMakie

    include("batchFiles/batchImages.jl")
    include("imageReader/imageReader.jl")
    include("imageReader/modelMaker.jl")


    export defineModel
    export read2DimageModel


end