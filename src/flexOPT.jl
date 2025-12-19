module flexOPT
    using KernelAbstractions
    using Adapt

    # GPU motors
    include("batchFiles/batchGPU.jl") # now indispensable for >3D problem for Coef computation
    export backend, makeGPUarray, recoverCPUarray


    # Semi symbolic OPT motors
    

end
