# look at out_of_src/Bspline_with_centralPoint_Symbolics.jl for some experiments
# that motivated to make this code



function BsplineTimesPolynomialsIntegrated(params::Dict)

    @unpack maximumOrder, numberNodes = params
    @variables x Δx ξ
   
    extFns = Symbolics.variables(:extFns,1:2,1:numberNodes,1:maximumOrder+1)
    # maximum order of B-spline
    
    maximumOrder = maximumOrder + 1

    

    ∂x = Differential(x)
    
    # input parameters
    
  
    
    # the left and right indices
    
    νₗ = 1 # we cannot touch this anymore
    νᵣ = numberNodes # like in the middle we know that the expression is ok 

    
    # let's start
    
    
    nodeIndices = collect(νₗ:1:νᵣ)
    nodesSymbolic = Δx .* nodeIndices
    append!(nodesSymbolic,nodesSymbolic[end])


    
    # b-splines
    
    b = zeros(Num, numberNodes, numberNodes, maximumOrder)
    bX = zeros(Num,numberNodes, 2, maximumOrder)
    
    # b-splines derivatives
    
    b_deriv = zeros(Num,numberNodes, numberNodes,maximumOrder,maximumOrder)
   
    for ι in 0:1:maximumOrder-1
        
        local neighbour = Int((1 - (-1)^ι) / 2)
        floor = Int((ι - neighbour) / 2)
        ceiling = Int((ι + neighbour) / 2)
        if ι === 0
            for ν in nodeIndices # this will run for all the ν related to nodes
                tmpν = ν - νₗ + 1
            
                b[tmpν, tmpν, ι+1] = 1
            end
        else
    
            for ν in nodeIndices # this will run for all the ν related to nodes
                tmpν = ν - νₗ + 1
                # the denominator for the ι>0
    
                denominator = BigInt(ι) * Δx
    
                # for the upgoing part
    
                if ν - neighbour <= νᵣ && ν - neighbour >= νₗ
                    rightlimit = minimum((numberNodes, tmpν + floor))
                    leftlimit = maximum((1, tmpν - ceiling))
                    for νSegment in leftlimit:1:rightlimit
                        tmpνSegment = νSegment #- νₗ + 1
                        numerator = x - BigInt(ν - ceiling) * Δx
                        b[tmpνSegment, tmpν, ι+1] += mySimplify(numerator / denominator * b[tmpνSegment, tmpν-neighbour, ι])
                    end
                end
    
                #for the downgoing part
    
                if ν - neighbour + 1 <= νᵣ && ν - neighbour + 1 >= νₗ
                    rightlimit = minimum((numberNodes, tmpν + floor + 1))
                    leftlimit = maximum((1, tmpν - ceiling + 1))
                    for νSegment in leftlimit:1:rightlimit
                        tmpνSegment = νSegment #- νₗ + 1
                        numerator = BigInt(ν + floor + 1) * Δx - x
                        b[tmpνSegment, tmpν, ι+1] += mySimplify(numerator / denominator * b[tmpνSegment, tmpν-neighbour+1, ι])
                    end
    
                end
    
            end
        end
    
        # unfortunately the νₗ and νᵣ related functions should recover the 'lost' amplitudes due to the trunctation 
        # we need to re-build based on the residual
    
        for ν in nodeIndices
            tmpν = ν - νₗ + 1
            
            #b[:,tmpν,ι+1] .= 0
            if ν===νₗ  # this loop runs only for the extremeties!!!!!
                # the nodes involved at most are: νₗ, νₗ + 1, ⋯, νₗ + ceiling 
                bX[1:ceiling+1,1,ι+1] .=1
                for νSegment in 1:1:ceiling+1
                    for νFunction in 2:1:numberNodes
                        bX[νSegment,1,ι+1] -= b[νSegment,νFunction,ι+1]
                    end
                end
            elseif ν===νᵣ
                bX[numberNodes-ceiling-1:numberNodes,2,ι+1] .=1
                for νSegment in numberNodes-ceiling-1:1:numberNodes
                    for νFunction in 1:1:numberNodes-1
                        bX[νSegment,2,ι+1] -= b[νSegment,νFunction,ι+1]
                    end
                end
            end
    
        end
    

        b[:,1,ι+1]=bX[:,1,ι+1]
        b[:,end,ι+1]=bX[:,end,ι+1]
        
          
        
        # computing the derivatives 
    
        for i in 0:1:maximumOrder-1 
            if i === 0
                b_deriv[:,:,i+1,ι+1] = b[:,:,ι+1]
            else
                b_deriv[:,:,i+1,ι+1] = mySimplify.(∂x.(b_deriv[:,:,i,ι+1]))
            
            end
    
        end
    
    end
        
    
    # re-calculated centre
    
    modμ=zeros(Num,3,numberNodes,maximumOrder) # min, max, mid for each b-splines this can be a multiple of Δx / 2 

    for ι in 0:1:maximumOrder-1
        for μ in nodeIndices
            leftLimit = nodesSymbolic[findfirst(x->x!==0,b[:,μ,ι+1])] /Δx
            rightLimit = nodesSymbolic[findlast(x->x!==0,b[:,μ,ι+1])+1] /Δx
            midNode = mySimplify((leftLimit+rightLimit)/2)
            modμ[1,μ,ι+1]=leftLimit 
            modμ[2,μ,ι+1]=rightLimit
            modμ[3,μ,ι+1]=midNode
        end
    end
    

    integral_b= zeros(Num,numberNodes,maximumOrder)

    #gvec = (@variables (g(x))[1:maximumOrder])[1]
    
    

    for ι in 0:1:maximumOrder-1
        for i in 0:1:maximumOrder-1
            for ν in nodeIndices # this will run for all the ν related to nodes
                tmpν = ν - νₗ + 1
                for νSegment in nodeIndices
                    tmpνSegment = νSegment - νₗ + 1
                    
                    integral_b[tmpν,ι+1] -= (-1)^(i)*extFns[1,tmpνSegment,i+1]*substitute(b_deriv[tmpνSegment,tmpν,i+1,ι+1],Dict(x=>nodesSymbolic[tmpνSegment]))
    
                    integral_b[tmpν,ι+1] += (-1)^(i)*extFns[2,tmpνSegment,i+1]*substitute(b_deriv[tmpνSegment,tmpν,i+1,ι+1],Dict(x=>nodesSymbolic[tmpνSegment+1]))
                    
                end
            end
        end
    end
    
    
    integral_b=mySimplify.(integral_b)
    
    BsplineIntegraters=(nodeIndices=nodeIndices,nodesSymbolic=nodesSymbolic,b_deriv=b_deriv,integral_b=integral_b,Δx=Δx,extFns=extFns,x=x,modμ=modμ)
    return @strdict(BsplineIntegraters)

end
    



function BsplineTimesPolynomialsIntegrated_deprecated(params::Dict)

    @unpack maximumOrder = params
    @variables x Δx ξ
    ∂x = Differential(x)
    
    # input parameters
    
    # maximum order of B-spline
    
    maximumOrder = maximumOrder + 1
    
    # the left and right indices
    νₗ = 0
    νᵣ = (maximumOrder+1)*2 + 3 # like in the middle we know that the expression is ok 
    νₘ = (νₗ+νᵣ) ÷ 2
    
    # let's start
    
    numberNodes = νᵣ - νₗ + 1 # [νₗ,νₗ + 1), ⋯, [νᵣ - 1 to νᵣ), {νᵣ}
    
    nodeIndices = collect(νₗ:1:νᵣ)
    
    nodesSymbolic = Δx .* nodeIndices
    append!(nodesSymbolic,nodesSymbolic[end])
    
    # b-splines
    
    b = zeros(Num, numberNodes, numberNodes, maximumOrder)
    bX = zeros(Num,numberNodes, 2, maximumOrder)
    
    # b-splines derivatives
    
    b_deriv = zeros(Num,numberNodes, numberNodes,maximumOrder,maximumOrder)
   
    for ι in 0:1:maximumOrder-1
        
        local neighbour = Int((1 - (-1)^ι) / 2)
        floor = Int((ι - neighbour) / 2)
        ceiling = Int((ι + neighbour) / 2)
        if ι === 0
            for ν in nodeIndices # this will run for all the ν related to nodes
                tmpν = ν - νₗ + 1
            
                b[tmpν, tmpν, ι+1] = 1
            end
        else
    
            for ν in nodeIndices # this will run for all the ν related to nodes
                tmpν = ν - νₗ + 1
                # the denominator for the ι>0
    
                denominator = BigInt(ι) * Δx
    
                # for the upgoing part
    
                if ν - neighbour <= νᵣ && ν - neighbour >= νₗ
                    rightlimit = minimum((numberNodes, tmpν + floor))
                    leftlimit = maximum((1, tmpν - ceiling))
                    for νSegment in leftlimit:1:rightlimit
                        tmpνSegment = νSegment #- νₗ + 1
                        numerator = x - BigInt(ν - ceiling) * Δx
                        b[tmpνSegment, tmpν, ι+1] += mySimplify(numerator / denominator * b[tmpνSegment, tmpν-neighbour, ι])
                    end
                end
    
                #for the downgoing part
    
                if ν - neighbour + 1 <= νᵣ && ν - neighbour + 1 >= νₗ
                    rightlimit = minimum((numberNodes, tmpν + floor + 1))
                    leftlimit = maximum((1, tmpν - ceiling + 1))
                    for νSegment in leftlimit:1:rightlimit
                        tmpνSegment = νSegment #- νₗ + 1
                        numerator = BigInt(ν + floor + 1) * Δx - x
                        b[tmpνSegment, tmpν, ι+1] += mySimplify(numerator / denominator * b[tmpνSegment, tmpν-neighbour+1, ι])
                    end
    
                end
    
            end
        end
    
        # unfortunately the νₗ and νᵣ related functions should recover the 'lost' amplitudes due to the trunctation 
        # we need to re-build based on the residual
    
        for ν in nodeIndices
            tmpν = ν - νₗ + 1
            
            #b[:,tmpν,ι+1] .= 0
            if ν===νₗ  # this loop runs only for the extremeties!!!!!
                # the nodes involved at most are: νₗ, νₗ + 1, ⋯, νₗ + ceiling 
                bX[1:ceiling+1,1,ι+1] .=1
                for νSegment in 1:1:ceiling+1
                    for νFunction in 2:1:numberNodes
                        bX[νSegment,1,ι+1] -= b[νSegment,νFunction,ι+1]
                    end
                end
            elseif ν===νᵣ
                bX[numberNodes-ceiling-1:numberNodes,2,ι+1] .=1
                for νSegment in numberNodes-ceiling-1:1:numberNodes
                    for νFunction in 1:1:numberNodes-1
                        bX[νSegment,2,ι+1] -= b[νSegment,νFunction,ι+1]
                    end
                end
            end
    
        end
    
        b[:,1,ι+1]=bX[:,1,ι+1]
        b[:,end,ι+1]=bX[:,end,ι+1]
    
        
        # computing the derivatives 
    
        for i in 0:1:maximumOrder-1 
            if i === 0
                b_deriv[:,:,i+1,ι+1] = b[:,:,ι+1]
            else
                b_deriv[:,:,i+1,ι+1] = mySimplify.(∂x.(b_deriv[:,:,i,ι+1]))
            
            end
    
        end
    
    end
        
    
    
    b_deriv_ξ = similar(b_deriv)
    for ι in 0:1:maximumOrder-1
        for i in 0:1:maximumOrder-1
            for ν in nodeIndices # this will run for all the ν related to nodes
                tmpν = ν - νₗ + 1
                for νSegment in nodeIndices
                    tmpνSegment = νSegment - νₗ + 1
                    b_deriv_ξ[tmpνSegment,tmpν,i+1,ι+1]=substitute(b_deriv[tmpνSegment,tmpν,i+1,ι+1],Dict(x=>ξ+Δx*(tmpν-1)))
                end
            end
        end
    end
    
    
    integral_b= zeros(Num,numberNodes)
    gvec = (@variables (g(x))[1:maximumOrder])[1]
    
    for ι in 0:1:maximumOrder-1
        for i in 0:1:maximumOrder-1
            for ν in nodeIndices # this will run for all the ν related to nodes
                tmpν = ν - νₗ + 1
                for νSegment in nodeIndices
                    tmpνSegment = νSegment - νₗ + 1
                    
                    diff = Δx * (tmpνSegment - tmpν)
                    tmpDic = Dict(x=>diff)
                    integral_b[tmpν] -= (-1)^(i+1)*substitute(gvec[i+1],tmpDic)*substitute(b_deriv[tmpνSegment,tmpν,i+1,ι+1],Dict(x=>nodesSymbolic[tmpνSegment]))
    
                    diff = Δx * (tmpνSegment - tmpν+1)
                    tmpDic = Dict(x=>diff)
                    integral_b[tmpν] += (-1)^(i+1)*substitute(gvec[i+1],tmpDic)*substitute(b_deriv[tmpνSegment,tmpν,i+1,ι+1],Dict(x=>nodesSymbolic[tmpνSegment+1]*Δx))
                    
                end
            end
        end
    end
    
    #display.(b_deriv)
    b_deriv_ξ=mySimplify.(b_deriv_ξ)
    #b_deriv_extremes=mySimplify.(b_deriv_extremes)
    integral_b=mySimplify.(integral_b)
    
    # for special g(x) = C ξ^N 
    @variables C N
    
    
    dictionaryForSubstitute =Dict()
    
    
    taylorNum = 1
    
    for i in 1:1:maximumOrder
        taylorNum *= N+i
        dictionaryForSubstitute[gvec[i]]=x^(N+i)/taylorNum
        
    end
    #@show dictionaryForSubstitute
    integral_b_polys= zeros(Num,numberNodes,maximumOrder)
    for ι in 0:1:maximumOrder-1
        for i in 0:1:maximumOrder-1
            for ν in nodeIndices # this will run for all the ν related to nodes
                tmpν = ν - νₗ + 1
                for νSegment in nodeIndices
                    tmpνSegment = νSegment - νₗ + 1
                    
                    tmp_b_deriv = b_deriv[tmpνSegment,tmpν,i+1,ι+1]
                    tmpG = substitute(gvec[i+1],dictionaryForSubstitute)
                        

                    if 1 <= tmpνSegment <= numberNodes-1
                        
                        diff = nodesSymbolic[tmpνSegment]-nodesSymbolic[tmpν]
                        tmpDic = Dict(x=>diff)

                        integral_b_polys[tmpν,ι+1] -= (-1)^(i)*substitute(tmpG,tmpDic)*substitute(tmp_b_deriv,Dict(x=>nodesSymbolic[tmpνSegment]))

                        diff = nodesSymbolic[tmpνSegment+1]-nodesSymbolic[tmpν]
                        tmpDic = Dict(x=>diff)

                        integral_b_polys[tmpν,ι+1] += (-1)^(i)*substitute(tmpG,tmpDic)*substitute(tmp_b_deriv,Dict(x=>nodesSymbolic[tmpνSegment+1]))

                    end



#=
                    diff = nodesSymbolic[tmpνSegment]-nodesSymbolic[tmpν]
                    tmpDic = Dict(x=>diff)
                    if tmpνSegment !== tmpν
                        integral_b_polys[tmpν,ι+1] -= (-1)^(i)*substitute(tmpG,tmpDic)*substitute(tmp_b_deriv,Dict(x=>nodesSymbolic[tmpνSegment]))
                    end
    
                    diff = nodesSymbolic[tmpνSegment+1]-nodesSymbolic[tmpν]
                    tmpDic = Dict(x=>diff)
                    if tmpνSegment+1 !== tmpν
                        integral_b_polys[tmpν,ι+1] += (-1)^(i)*substitute(tmpG,tmpDic)*substitute(tmp_b_deriv,Dict(x=>nodesSymbolic[tmpνSegment+1]))
                    end
=#


                    
                end
            end
        end
    end

    #=
    individualIntegral_b_polys = copy(integral_b_polys)
    integral_b_polys= zeros(Num,numberNodes,maximumOrder)

    for ι in 0:1:maximumOrder-1
        for ν in nodeIndices # this will run for all the ν related to nodes
            tmpν = ν - νₗ + 1
            for i in 0:ι
                integral_b_polys[tmpν,ι+1]+=individualIntegral_b_polys[tmpν,i+1]
            end
        end
    end
    =#

    integral_b_polys=mySimplify.(integral_b_polys)
    BsplineIntegraters=(numberNodes=numberNodes,integral_b_polys=integral_b_polys,N=N,Δx =Δx)
    return @strdict(BsplineIntegraters)
end
    