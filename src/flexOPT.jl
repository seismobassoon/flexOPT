module flexOPT

    #include("imageReader/imageReader.jl")
    include("../batchFiles/batchGPU.jl") # now indispensable for >3D problem for Coef computation

end
