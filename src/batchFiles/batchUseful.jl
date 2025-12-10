# In this code, I collect some useful macros and functions 
using LinearAlgebra, JLD2

#good old safeget
safeget(A, inds...; default=0) = checkbounds(Bool, A, inds...) ? A[inds...] : default


# macro to convert a variable name to a string 
# https://discourse.julialang.org/t/can-i-use-a-function-variable-name-to-generate-a-string-to-set-it-equal-to/73659
macro var2string(var)
    :($(esc(var)) = $(String(var)))
end



function reinterpolateArrayMembers(Ncolor,nSegments)
    nColors = Array{Int}(undef,nSegments)
    smallNumberMembers = Ncolor ÷ nSegments
    oneMoreMemberSegment = Ncolor % nSegments


    nColors .=smallNumberMembers
    nColors[1:oneMoreMemberSegment] .= smallNumberMembers+1
    
    intervalIndice=Array{Int,2}(undef,(2,nSegments))

    for i in 1:nSegments
        if i == 1
            intervalIndice[1,i] = 1
            intervalIndice[2,i] = nColors[1]
        else
            intervalIndice[1,i] = intervalIndice[1,i-1]+nColors[i-1]
            intervalIndice[2,i] = intervalIndice[1,i]+nColors[i]-1
        end
    end

    return nColors,intervalIndice
end

function expandVectors(vector,onesvector;fillinTheGap=1,Type=Int)
    # this function expands a <N dimension vector on N dimension onesvector 
    vector=collect(vector)
    onesvector=collect(onesvector)
    if size(vector)[1] !== sum(onesvector)
        @show vector,onesvector
        @error "the input vector dimension is not what onesvector can help"
    end
    Ndimension=size(onesvector)[1]
    NewVector=fillinTheGap .* ones(Type,Ndimension)
    iCoord = 0
    
    for jCoord in 1:Ndimension # here I cannot use eachindex since I need to be quite sure about the order of coordinates
        if onesvector[jCoord] === 1
            iCoord +=1
            NewVector[jCoord]=vector[iCoord]
        end
    end
    return NewVector
end

function string_as_varname(s::AbstractString,v::Any)
    s=Symbol(s)
    @eval (($s) = ($v))
end

car2vec(x::CartesianIndex) = collect(Tuple(x))
carDropDim(x::CartesianIndex) = CartesianIndex(Tuple(car2vec(x)[1:end-1]))
carAddDim(x::CartesianIndex,n::Int) = CartesianIndex(Tuple([car2vec(x);n]))
vec2car(x::Array{Int})=CartesianIndex(Tuple(vec(x))) 

function flatten(x)
    if isa(x, AbstractVector)
        return vcat(flatten.(x)...)
    else
        return [x]
    end
end

function deep_flatten(x)
    if isa(x, Number)       # base case
        return [x]
    elseif isa(x, AbstractArray)
        return vec(x)       # convert any array to 1D vector
    elseif isa(x, AbstractVector)
        return vcat(deep_flatten.(x)...)
    else
        error("Unsupported type: $(typeof(x))")
    end
end


function vec2car(x::Vector{Int}, Ndimension::Int)
    @assert length(x) % Ndimension == 0 "Flat array length not divisible by dimension"

    M = reshape(x, Ndimension, :)
    return [CartesianIndex(Tuple(M[:,i])) for i in axes(M, 2)]
end

function is_all_less_than_or_equal(c1::CartesianIndex, c2::CartesianIndex)
    all(x -> x[1] ≤ x[2], zip(Tuple(c1), Tuple(c2)))
end

function distance2_point_to_box(p::CartesianIndex, c1::CartesianIndex, c2::CartesianIndex)
    lower = min.(Tuple(c1), Tuple(c2))
    upper = max.(Tuple(c1), Tuple(c2))
    point = Tuple(p)

    # For each dimension, compute distance to box surface
    distsq = 0.0
    for i in eachindex(point)
        if point[i] < lower[i]
            distsq += (lower[i] - point[i])^2
        elseif point[i] > upper[i]
            distsq += (point[i] - upper[i])^2
        else
            # inside the slab: distance = 0 for this dimension
        end
    end

    return distsq
end


function myInv(a;method="LU")
    ainv=nothing
    if method==="macGPU"# this does not work for the moment
        a=convert(Array{Float32},a)
        a_gpu = MtlArray(a)  # Move A to GPU
        ainv = inv(a_gpu) 
        ainv=convert(Array{Float64},ainv)
    elseif method==="LU"
        F=lu(a)
        ainv=inv(F)
    elseif method==="basic"
        ainv=inv(a)
    end
    return ainv
end

