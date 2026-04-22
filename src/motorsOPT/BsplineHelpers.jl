using Symbolics
using CairoMakie

# this source file was truly a fruit of some weeks of optimisations with codex (chatGPT) March-April 2026
# Nobuaki Fuji, IPGP/UPC/IUF

"""
    constructBsplineFamily(params; simplify_expr=Symbolics.simplify, boundary_mode=:clamped)

Build a limited 1D B-spline family on the node grid described by
`(maximumOrder, allNodes, idx_nodesNum, idx_refPoints)`.

The construction:

- augments the requested reference points with numerical left/right closures
  when needed,
- pads the resulting knot list with `maximumOrder` auxiliary knots on each side,
- constructs the full auxiliary family,
- returns only the centre functions corresponding to the original
  `idx_refPoints`.

The returned `b` array stores piecewise symbolic expressions with layout
`(segment, function, order_slot)`, where `order_slot = p + 2` corresponds to
degree `p`.
"""
function constructBsplineFamily(params; simplify_expr=mySimplify, boundary_mode=:ghost, correction_truncation = true)
    

    maximumOrder = params.maximumOrder
    allNodes = collect(params.allNodes)
    idx_nodesNum = collect(params.idx_nodesNum) # an ordinary consecutive integer increment from 1 (the numerical nodes with Δy)
    idx_refPoints_original = collect(params.idx_refPoints) # supporting nodes to construct the Bspline family
    idx_selectedPoints = collect(params.idx_selectedPoints) # the node addresses that user needs to take 

    if !issubset(Set(idx_selectedPoints),Set(idx_refPoints_original))
        @error "idx_selectedPoints should be a subset of idx_refPoints_original"
    end

    selected_positions = indexin(Int.(idx_selectedPoints), Int.(idx_refPoints_original))

    @variables x Δx
    ∂x = Differential(x)

    numberNodes = length(allNodes)
    NorderShiftedByTwo = maximumOrder + 2

    core_idx_refPoints = copy(idx_refPoints_original)
    leftPlusNumberFunctions = 0
    rightPlusNumberFunctions = 0


    if allNodes[core_idx_refPoints[1]] > allNodes[idx_nodesNum[1]]
        core_idx_refPoints = vcat(idx_nodesNum[1], core_idx_refPoints)
        leftPlusNumberFunctions = 1
    end

    if allNodes[core_idx_refPoints[end]] < allNodes[idx_nodesNum[end]]
        core_idx_refPoints = vcat(core_idx_refPoints, idx_nodesNum[end])
        rightPlusNumberFunctions = 1
    end

    coreKnots = allNodes[core_idx_refPoints]
    refKnots = extend_knots(coreKnots, maximumOrder; mode=boundary_mode)
    refKnotsSymbolic = Δx .* refKnots

    numberFunctionsOriginal = length(idx_refPoints_original)
    numberFunctions = length(refKnots) -1
    denominators = zeros(Num, 2, numberFunctions, maximumOrder)
    for p in 1:maximumOrder
        for idx in 1:numberFunctions-p
            denominators[1, idx, p] = refKnotsSymbolic[idx + p] - refKnotsSymbolic[idx]
            denominators[2, idx, p] = refKnotsSymbolic[idx + p + 1] - refKnotsSymbolic[idx + 1]
        end
    end

    

    b = zeros(Num, numberNodes, numberFunctions, NorderShiftedByTwo)
    b_to_give = zeros(Num, numberNodes, numberFunctionsOriginal, NorderShiftedByTwo)
    b_deriv = zeros(Num, numberNodes, numberFunctions, NorderShiftedByTwo, NorderShiftedByTwo)
    b_deriv_to_give = zeros(Num, numberNodes, numberFunctionsOriginal, NorderShiftedByTwo, NorderShiftedByTwo)
    
    # -1 order indicator on the numerical support.
    slot = -1 + 2
    b[idx_nodesNum[1]:idx_nodesNum[end], :, 1] .= 1
    firstCentral = 1 + leftPlusNumberFunctions + maximumOrder
    lastCentral = firstCentral + numberFunctionsOriginal - 1
    centre_range = firstCentral:lastCentral
    b_to_give[:,:,slot] = b[:,centre_range,slot]
    b_deriv_to_give[:,:,:,slot] = b_deriv[:,centre_range,:,slot]



    incrementDueToCentreShift = 0
    for p in 0:maximumOrder
        slot = p + 2

        if p == 0
            for idx in 1:numberFunctions
                leftKnot = refKnots[idx]
                rightKnot = refKnots[idx + 1]
                support_segments = support_segment_indices(leftKnot, rightKnot, allNodes)
                for seg in support_segments
                    b[seg, idx, slot] = 1
                end
            end
        else
            b[:, :, slot] .= 0
            for idx in 1:numberFunctions-p
                up = denominators[1, idx, p]
                down = denominators[2, idx, p]
                leftKnot = refKnotsSymbolic[idx]
                rightKnot = refKnotsSymbolic[idx + p + 1]

                active_segments_left = support_segment_indices(refKnots[idx], refKnots[idx + p], allNodes)
                active_segments_right = support_segment_indices(refKnots[idx + 1], refKnots[idx + p + 1], allNodes)

                for seg in active_segments_left
                    if !iszero_denominator(up)
                        b[seg, idx, slot] += simplify_expr((x - leftKnot) / up * b[seg, idx, slot - 1])
                    end
                end

                for seg in active_segments_right
                    if !iszero_denominator(down)
                        b[seg, idx, slot] += simplify_expr((rightKnot - x) / down * b[seg, idx + 1, slot - 1])
                    end
                end
            end

            if boundary_mode == :clamped
                rebuild_clamped_residuals!(b, allNodes, refKnots, numberFunctions, slot, p, simplify_expr)
            end

        end

        for pp in 0:1:p
            if pp == 0
                b_deriv[:,:,pp+1,slot] = b[:,:,slot]
            else
                b_deriv[:,:,pp+1,slot] = mySimplify.(∂x.(b_deriv[:,:,pp,slot]))
            
            end
        end

    
        incrementDueToCentreShift -= Int((1 - (-1)^(p)) / 2)
        
        @show firstCentral = leftPlusNumberFunctions + maximumOrder + 1 + incrementDueToCentreShift  
        @show lastCentral = firstCentral + numberFunctionsOriginal - 1
        centre_range = firstCentral:lastCentral


        b_to_give[:,:,slot] = b[:,centre_range,slot]
        b_deriv_to_give[:,:,:,slot] = b_deriv[:,centre_range,:,slot]


        if correction_truncation && p > 0
            left_support = support_segment_indices(refKnots[firstCentral], refKnots[firstCentral + p + 1], allNodes)
            right_support = support_segment_indices(refKnots[lastCentral], refKnots[lastCentral + p + 1], allNodes)

            b_to_give[:, 1, slot] .= 0
            b_to_give[:, end, slot] .= 0

            for seg in left_support
                s = zero(Num)
                for j in 2:numberFunctionsOriginal
                    s += b_to_give[seg, j, slot]
                end
                b_to_give[seg, 1, slot] = simplify_expr(1 - s)
            end

            for seg in right_support
                s = zero(Num)
                for j in 1:numberFunctionsOriginal-1
                    s += b_to_give[seg, j, slot]
                end
                b_to_give[seg, end, slot] = simplify_expr(1 - s)
            end


            for pp in 0:1:p
                if pp == 0
                    b_deriv_to_give[:,1,pp+1,slot] = b_to_give[:,1,slot]
                    b_deriv_to_give[:,end,pp+1,slot] = b_to_give[:,end,slot]
                else
                    b_deriv_to_give[:,1,pp+1,slot] = mySimplify.(∂x.(b_deriv_to_give[:,1,pp,slot]))
                    b_deriv_to_give[:,end,pp+1,slot] = mySimplify.(∂x.(b_deriv_to_give[:,end,pp,slot]))
                end
            end

        end







    end

    #b_selected = b_to_give[:, selected_positions, :]
    b_deriv_selected = b_deriv_to_give[:, selected_positions, :, :]



    return (
        #b = b_to_give,
        #b_full = b,
        b = b_deriv_selected,
        b_full = b_deriv,
        b_support = b_deriv_to_give,
        x = x,
        Δx = Δx,
        allNodes = allNodes,
        core_idx_refPoints = core_idx_refPoints,
        coreKnots = coreKnots,
        refKnots = refKnots,
        boundary_mode = boundary_mode,
    )
end



function extend_knots(coreKnots, maximumOrder; mode=:clamped)
    if mode == :clamped
        return vcat(
            fill(coreKnots[1], maximumOrder),
            coreKnots,
            fill(coreKnots[end], maximumOrder),
        )
    elseif mode == :ghost
        left_step, right_step = edge_steps(coreKnots)
        left_ghosts = [coreKnots[1] - k * left_step for k in maximumOrder:-1:1]
        right_ghosts = [coreKnots[end] + k * right_step for k in 1:maximumOrder]
        return vcat(left_ghosts, coreKnots, right_ghosts)
    else
        throw(ArgumentError("unknown boundary_mode=$mode; use :clamped or :ghost"))
    end
end

function edge_steps(coreKnots)
    if length(coreKnots) == 1
        return (1.0, 1.0)
    end

    return (coreKnots[2] - coreKnots[1], coreKnots[end] - coreKnots[end - 1])
end

function support_segment_indices(leftKnot, rightKnot, allNodes)
    if rightKnot <= leftKnot
        return Int[]
    end

    segments = Int[]
    for seg in 1:length(allNodes)-1
        leftNode = allNodes[seg]
        rightNode = allNodes[seg + 1]
        if max(leftNode, leftKnot) < min(rightNode, rightKnot)
            push!(segments, seg)
        end
    end
    return segments
end


function iszero_denominator(expr)
    try
        return Bool(Symbolics.value(iszero(expr)))
    catch
        return false
    end
end

function rebuild_clamped_residuals!(b, allNodes, refKnots, numberFunctions, slot, p, simplify_expr)
    left_support = support_segment_indices(1, refKnots[p + 2], allNodes)
    for seg in left_support
        s = zero(Num)
        for idx in 2:numberFunctions
            s += b[seg, idx, slot]
        end
        b[seg, 1, slot] = simplify_expr(1 - s)
    end

    rightExtreme = length(allNodes)

    right_support = support_segment_indices(refKnots[numberFunctions-p-1], rightExtreme, allNodes)
    for seg in right_support
        s = zero(Num)
        for idx in 1:numberFunctions-1
            s += b[seg, idx, slot]
        end
        b[seg, numberFunctions, slot] = simplify_expr(1 - s)
    end

    return b
end

"""
    evaluate_bspline_piecewise(b, idx, orderSlot, ξ, allNodes, x, Δx; Δxval=1.0)

Evaluate the piecewise symbolic basis function stored in `b[:, idx, orderSlot]`
at the normalized coordinate `ξ = x / Δx`.
"""
function evaluate_bspline_piecewise(b, idx, orderSlot, ξ, allNodes, x, Δx; Δxval=1.0)
    segment = find_segment(ξ, allNodes)
    expr = b[segment, idx, orderSlot]
    value = Symbolics.substitute(expr, Dict(x => ξ * Δxval, Δx => Δxval))
    return Symbolics.value(value)
end

function find_segment(ξ, allNodes)
    for k in 1:length(allNodes)-1
        if allNodes[k] <= ξ < allNodes[k + 1]
            return k
        end
    end
    return length(allNodes)
end

"""
    plot_bspline_family(result; order=result_order(result), N=10, Δxval=1.0, show_full=false)

Sample and plot the symbolic B-spline family returned by `construct_bspline_family`.
The x-axis is the normalized coordinate `ξ = x / Δx`, with `N` samples per unit span.
"""
function segment_grid(allNodes; N=10)
    ξgrid = Float64[]
    for k in 1:length(allNodes)-1
        a = allNodes[k]
        b = allNodes[k+1]
        localgrid = collect(range(a, b; length=N+1))
        if k > 1
            localgrid = localgrid[2:end]  # avoid duplicating left endpoint
        end
        append!(ξgrid, localgrid)
    end
    return ξgrid
end
function plot_bspline_family(
    result;
    derivOrder=0,
    N=10,
    Δxval=1.0,
    show_full=false,
)
    order = size(result.b, 4) - 2
    arr = show_full ? result.b_full : result.b
    ylabel_text = derivOrder == 0 ? "B-spline value" : "Derivative order $derivOrder"

    allNodes = result.allNodes
    x = result.x
    Δx = result.Δx

    orderSlot = order + 2
    derivSlot = derivOrder + 1

    ξgrid = segment_grid(allNodes; N=N)

    fig = Figure(size=(900, 500))
    ax = Axis(fig[1, 1], xlabel="x / Δx", ylabel=ylabel_text)

    for idx in 1:size(arr, 2)
        vals = [
            evaluate_bspline_piecewise_deriv(arr, idx, derivSlot, orderSlot, ξ, allNodes, x, Δx; Δxval=Δxval)
            for ξ in ξgrid
        ]
        lines!(ax, ξgrid, vals, label="B_$idx")
    end

    vlines!(ax, allNodes; color=:gray70, linestyle=:dash)
    axislegend(ax)
    return fig
end



function evaluate_bspline_piecewise_deriv(b_deriv, idx, derivSlot, orderSlot, ξ, allNodes, x, Δx; Δxval=1.0)
    segment = find_segment(ξ, allNodes)
    expr = b_deriv[segment, idx, derivSlot, orderSlot]
    value = Symbolics.substitute(expr, Dict(x => ξ * Δxval, Δx => Δxval))
    return Symbolics.value(value)
end

