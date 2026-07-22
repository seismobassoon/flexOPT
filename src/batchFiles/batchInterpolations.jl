# Utilities for DIVAnd interpolation on regular and deformed query grids.

correlationLengthDefault = (20e3, 20e3, 20e3)
epsilon2Default = 1.0

_axis_spec(values::NamedTuple) = values
function _axis_spec(values::AbstractVector)
    lo, hi = extrema(values)
    return (min=lo, max=hi, n=length(values))
end

giveMaskPX(Xs, Ys) = giveMaskPX(Xs, Ys, nothing)

function giveMaskPX(Xs, Ys, Zs)
    xs = _axis_spec(Xs)
    ys = _axis_spec(Ys)
    minX = hasproperty(xs, :minX) ? xs.minX : xs.min
    maxX = hasproperty(xs, :maxX) ? xs.maxX : xs.max
    nX = hasproperty(xs, :nX) ? xs.nX : xs.n
    minY = hasproperty(ys, :minY) ? ys.minY : ys.min
    maxY = hasproperty(ys, :maxY) ? ys.maxY : ys.max
    nY = hasproperty(ys, :nY) ? ys.nY : ys.n

    if Zs === nothing
        mask, (pm, pn), (xi, yi) = DIVAnd_rectdom(
            range(minX, maxX; length=nX),
            range(minY, maxY; length=nY),
        )
        return (mask=mask, Ps=(pm, pn), XIs=(xi, yi))
    end

    zs = _axis_spec(Zs)
    minZ = hasproperty(zs, :minZ) ? zs.minZ : zs.min
    maxZ = hasproperty(zs, :maxZ) ? zs.maxZ : zs.max
    nZ = hasproperty(zs, :nZ) ? zs.nZ : zs.n
    mask, (pm, pn, po), (xi, yi, zi) = DIVAnd_rectdom(
        range(minX, maxX; length=nX),
        range(minY, maxY; length=nY),
        range(minZ, maxZ; length=nZ),
    )
    return (mask=mask, Ps=(pm, pn, po), XIs=(xi, yi, zi))
end

interpolateField(rawField, Xs, Ys, Xnode, Ynode;
    correlationLength=(correlationLengthDefault[1], correlationLengthDefault[3]),
    epsilon2=epsilon2Default,
) = interpolateField(
    rawField, Xs, Ys, nothing, Xnode, Ynode, nothing;
    correlationLength=correlationLength,
    epsilon2=epsilon2,
)

function interpolateField(
    rawField, Xs, Ys, Zs, Xnode, Ynode, Znode;
    correlationLength=correlationLengthDefault,
    epsilon2=epsilon2Default,
)
    masks = giveMaskPX(Xs, Ys, Zs)
    XXnodes = Znode === nothing ? (Xnode, Ynode) : (Xnode, Ynode, Znode)
    fieldCartesian, _ = DIVAndrun(
        masks.mask,
        masks.Ps,
        masks.XIs,
        XXnodes,
        vec(rawField),
        correlationLength,
        epsilon2,
    )
    return fieldCartesian
end

function _fast_analysis_size(Xquery, Yquery, correlationLength, analysis_size, max_analysis_points)
    analysis_size !== nothing && return analysis_size
    minX, maxX = extrema(Xquery)
    minY, maxY = extrema(Yquery)
    nX = max(2, ceil(Int, (maxX - minX) / correlationLength[1]) + 1)
    nY = max(2, ceil(Int, (maxY - minY) / correlationLength[2]) + 1)
    if nX * nY > max_analysis_points
        scale = sqrt(max_analysis_points / (nX * nY))
        nX = max(2, floor(Int, nX * scale))
        nY = max(2, floor(Int, nY * scale))
    end
    return (nX, nY)
end

"""
    interpolateField(rawField, Xquery::AbstractMatrix, Yquery::AbstractMatrix,
                     Xnode, Ynode; ...)

Run DIVAnd on a bounded regular analysis grid and then sample that result at
non-separable point-wise query coordinates. This avoids constructing and
solving the DIVAnd system at the full deformed-grid resolution.
"""
function interpolateField(
    rawField,
    Xquery::AbstractMatrix{<:Real},
    Yquery::AbstractMatrix{<:Real},
    Xnode,
    Ynode;
    correlationLength=(correlationLengthDefault[1], correlationLengthDefault[3]),
    epsilon2=epsilon2Default,
    analysis_size=nothing,
    max_analysis_points=500_000,
)
    size(Xquery) == size(Yquery) || throw(DimensionMismatch(
        "Xquery and Yquery must have the same size",
    ))
    nX, nY = _fast_analysis_size(
        Xquery, Yquery, correlationLength, analysis_size, max_analysis_points,
    )
    minX, maxX = extrema(Xquery)
    minY, maxY = extrema(Yquery)
    Xs = (minX=minX, maxX=maxX, nX=nX)
    Ys = (minY=minY, maxY=maxY, nY=nY)
    analysis = interpolateField(
        rawField, Xs, Ys, Xnode, Ynode;
        correlationLength=correlationLength,
        epsilon2=epsilon2,
    )

    xaxis = range(minX, maxX; length=nX)
    yaxis = range(minY, maxY; length=nY)
    gridded = interpolate((xaxis, yaxis), analysis, Gridded(Linear()))
    sampler = extrapolate(gridded, Flat())
    result = Array{eltype(analysis)}(undef, size(Xquery))
    Base.Threads.@threads for index in eachindex(result)
        result[index] = sampler(Xquery[index], Yquery[index])
    end
    return result
end
