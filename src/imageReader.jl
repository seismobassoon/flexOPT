module imageReader
    using FileIO,CairoMakie

    include("batchFiles/batchImages.jl")

    include("imageReader/imageReader.jl")

    export read2DimageModel, defineModel

end