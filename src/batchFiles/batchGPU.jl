


# ---------------------------------------------------------------------
# Robust backend detection (CUDA → Metal → CPU)
function detect_backend()
    if @isdefined(CUDA) && CUDA.has_cuda()
        println("→ Using CUDA backend")
        return CUDABackend()
    elseif @isdefined(Metal)
        try
            devs = Metal.devices()
            if !isempty(devs)
                println("→ Using Metal backend (", length(devs), " device(s))")
                return MetalBackend()
            end
        catch err
            @warn "Metal available but cannot query devices: $err"
        end
    end
    println("→ Using CPU backend (no GPU detected)")
    return CPU()
end

@kernel function piecewiseProduct!(a,b,result)
    i = @index(Global)
    result[i]=a[i]*b[i]
end

@kernel function piecewiseProduct!(a,b,c,result)
    i = @index(Global)
    result[i]=a[i]*b[i]*c[i]
end

@kernel function piecewiseProduct!(a,b,c,d,result)
    i = @index(Global)
    result[i]=a[i]*b[i]*c[i]*d[i]
end

@kernel function piecewiseProduct!(a,b,c,d,e,result)
    i = @index(Global)
    result[i]=a[i]*b[i]*c[i]*d[i]*e[i]
end

@kernel function piecewiseProductAll!(A, B, C, result, ::Val{N}) where {N}
    i = @index(Global)
    temp = one(eltype(result))
    for n in 1:N
        temp *= A[n, i]
    end
    temp *= B[i]*C[i]
    result[i] = temp
end

function makeGPUarray(backend,A)
    A_gpu= Adapt.adapt(backend, Float32.(A))
    return A_gpu
end
function recoverCPUarray(A_gpu)
    A_cpu = Array(A_gpu)         # preserves shape
    return Float64.(A_cpu)       # convert elements
end

function GPUsum(gpuArray)
    return KernelAbstractions.sum(gpuArray)
end

const backend = detect_backend()
println("Selected backend type: ", typeof(backend))

