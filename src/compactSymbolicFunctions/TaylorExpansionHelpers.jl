using Symbolics


function constructTaylorExpansions(params;simplify_expr=mySimplify,)
    # this function constructs the Taylor expansion expressions from given μ points ("selectedPoints)

    maxL_MINUS_N= params.maxL_MINUS_N
    allNodes = collect(params.allNodes)
    idx_nodesNum = collect(params.idx_nodesNum)
    idx_selectedPoints = collect(params.idx_selectedPoints)

    numberFunctions = length(idx_selectedPoints)
    numberL_MINUS_N = 1+maxL_MINUS_N

    @variables x Δx
    variables = Num[x, Δx]
    rangeSegments = idx_nodesNum[1]:idx_nodesNum[end]
    segmentNodes = allNodes[rangeSegments]

    allNodesSymbolic = Δx .* allNodes


    taylorKernel = CompactSymbolicFunctions(
        segmentNodes,
        numberFunctions;
        auxDims=(numberL_MINUS_N,),
        variables=variables,
    )

    for l_n in 0:maxL_MINUS_N
        slot = l_n + 1
        denominator=1/factorial(BigInt(l_n))
        for iμ in eachindex(idx_selectedPoints)
            pointCoord = allNodesSymbolic[iμ]
            taylorKernel.data[:,iμ,slot].=(x-pointCoord)^l_n*denominator
        end
    end
    return (k=taylorKernel,)
end

function plotTaylorExpansions(csf::CompactSymbolicFunctions;l_n=0,N=10, Δxval=1.0) 
    slots=(l_n+1,)
    ylabel="Taylor expansion of the degree of $l_n"
    return plotCompactSymbolicFunctions(csf, slots::NTuple; N=N, Δxval=Δxval,ylabel=ylabel)
end