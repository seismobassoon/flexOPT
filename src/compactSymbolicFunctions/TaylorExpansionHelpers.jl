using Symbolics


function constructTaylorExpansions(params;simplify_expr=mySimplify,)
    # this function constructs the Taylor expansion expressions from given μ points ("selectedPoints)

    maxL_MINUS_N= params.maxL_MINUS_N
    allNodes = collect(params.allNodes)
    idx_nodesNum = collect(params.idx_nodesNum)
    idx_selectedPoints = collect(params.idx_selectedPoints)

    @variables x Δx
    variables = Num[x, Δx]
    rangeSegments = idx_nodesNum[1]:idx_nodesNum[end]
    segmentNodes = allNodes[rangeSegments]

    allNodesSymbolic = Δx .* allNodes

    for iμ in idx_selectedPoints
        pointCoord = allNodesSymbolic[iμ]
        
    end

end