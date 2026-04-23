using Symbolics
using CairoMakie
using Base
import Base: ==,+, -, *, /, getindex, size

include("../batchFiles/batchSymbolics.jl")

# again, this is not completely handmade but it works nicely... (April 2026 Nobuaki Fuji)

struct CompactSymbolicFunctions{N,A<:AbstractArray{Num}}
    numberNodes::Int
    numberFunctions::Int
    auxDims::NTuple{N,Int}
    variables::Vector{Num}
    nodes::Vector{Num}
    data::A
end

function CompactSymbolicFunctions(
    nodes::AbstractVector,
    numberFunctions::Int;
    auxDims::Tuple=(),
    variables::AbstractVector=Num[],
    init=nothing,
)
    N = length(auxDims)
    dims = (length(nodes), numberFunctions, auxDims...)
    data = isnothing(init) ? zeros(Num, dims) : collect(init)

    #ndims(data) == 2 + N || error("data must have $(2+N) dimensions")
    size(data) == dims || error("data must have size $dims")

    return CompactSymbolicFunctions{N,typeof(data)}(
        length(nodes),
        numberFunctions,
        auxDims,
        collect(variables),
        collect(Num.(nodes)),
        data,
    )
end

function Base.:(==)(a::CompactSymbolicFunctions, b::CompactSymbolicFunctions)
    a.numberNodes == b.numberNodes &&
    a.numberFunctions == b.numberFunctions &&
    a.auxDims == b.auxDims &&
    size(a.data) == size(b.data) &&
    same_vector(a.nodes, b.nodes) &&
    same_vector(a.variables, b.variables) &&
    same_vector(vec(a.data), vec(b.data))
end

function Base.:(+)(a::CompactSymbolicFunctions{N}, b::CompactSymbolicFunctions{N}) where {N}
    check_compatible(a, b)
    CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=a.auxDims,
        variables=a.variables,
        init=a.data .+ b.data,
    )
end

function Base.:(-)(a::CompactSymbolicFunctions{N}, b::CompactSymbolicFunctions{N}) where {N}
    check_compatible(a, b)
    CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=a.auxDims,
        variables=a.variables,
        init=a.data .- b.data,
    )
end


function Base.:(*)(c, a::CompactSymbolicFunctions)
    CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=a.auxDims,
        variables=a.variables,
        init=c .* a.data,
    )
end

Base.:(*)(a::CompactSymbolicFunctions, c) = c * a

function Base.:(/)(a::CompactSymbolicFunctions, c)
    CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=a.auxDims,
        variables=a.variables,
        init=a.data ./ c,
    )
end

function Base.:(*)(a::CompactSymbolicFunctions{N}, b::CompactSymbolicFunctions{N}) where {N}
    check_compatible(a, b)
    CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=a.auxDims,
        variables=a.variables,
        init=a.data .* b.data,
    )
end

function Base.:(/)(a::CompactSymbolicFunctions{N}, b::CompactSymbolicFunctions{N}) where {N}
    check_compatible(a, b)
    CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=a.auxDims,
        variables=a.variables,
        init=a.data ./ b.data,
    )
end

function differentiate(a::CompactSymbolicFunctions, var, maxIter; simplify_expr=mySimplify)
    D = Differential(var)
    newAuxDims = (maxIter + 1, a.auxDims...)
    newCSF = CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=newAuxDims,
        variables=a.variables,
    )

    # store the original object in derivative slot 1
    selectdim(newCSF.data, 3, 1) .= a.data

    current = copy(a.data)
    for i in 1:maxIter
        current = simplify_expr.(D.(current))
        selectdim(newCSF.data, 3, i + 1) .= current
    end

    return newCSF
end


function integrate(a::CompactSymbolicFunctions, var, maxIter; simplify_expr=identity)
    newAuxDims = (maxIter + 1, a.auxDims...)
    newCSF = CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=newAuxDims,
        variables=a.variables,
    )

    selectdim(newCSF.data, 3, 1) .= a.data

    current = copy(a.data)
    for i in 1:maxIter
        current = simplify_expr.(integrateTaylorPolynomials.(current, Ref(var)))
        selectdim(newCSF.data, 3, i + 1) .= current
    end

    return newCSF
end


function differentiate(a::CompactSymbolicFunctions, var; simplify_expr=mySimplify)
    D = Differential(var)
    CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=a.auxDims,
        variables=a.variables,
        init=simplify_expr.(D.(a.data)),
    )
end

function integrate(a::CompactSymbolicFunctions, var; simplify_expr=mySimplify)
    CompactSymbolicFunctions(
        a.nodes,
        a.numberFunctions;
        auxDims=a.auxDims,
        variables=a.variables,
        #init=simplify_expr.(Symbolics.integrate.(a.data, Ref(var))),
        init=simplify_expr.(integrateTaylorPolynomials.(a.data, Ref(var)))
    )
end

function same_vector(a, b)
    length(a) == length(b) || return false
    all(isequal.(a, b))
end

function check_compatible(a::CompactSymbolicFunctions, b::CompactSymbolicFunctions)
    a.numberNodes == b.numberNodes || error("numberNodes mismatch")
    a.numberFunctions == b.numberFunctions || error("numberFunctions mismatch")
    a.auxDims == b.auxDims || error("auxDims mismatch")
    same_vector(a.nodes, b.nodes) || error("nodes mismatch")
    same_vector(a.variables, b.variables) || error("variables mismatch")
end


Base.getindex(a::CompactSymbolicFunctions, I...) = a.data[I...]
Base.size(a::CompactSymbolicFunctions) = size(a.data)


function plotCompactSymbolicFunctions(csf; order=0, derivOrder=0, N=10, Δxval=1.0)
    allNodes = Symbolics.value.(csf.nodes)
    x = csf.variables[1]
    Δx = csf.variables[2]

    orderSlot = order + 1
    derivSlot = derivOrder + 1
    ξgrid = segment_grid(allNodes; N=N)

    fig = Figure(size=(900, 500))
    ax = Axis(fig[1, 1], xlabel="x / Δx", ylabel="Derivative order $derivOrder")

    for idx in 1:csf.numberFunctions
        vals = [
            begin
                seg = find_segment(ξ, allNodes)
                expr = csf.data[seg, idx, derivSlot, orderSlot]
                Symbolics.value(Symbolics.substitute(expr, Dict(x => ξ * Δxval, Δx => Δxval)))
            end
            for ξ in ξgrid
        ]
        lines!(ax, ξgrid, vals, label="B_$idx")
    end

    vlines!(ax, allNodes; color=:gray70, linestyle=:dash)
    axislegend(ax)
    fig
end
