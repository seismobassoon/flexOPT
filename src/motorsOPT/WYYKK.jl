# these functions prepare analytically the 1D integral of WYYKK (cf. eq 52 of FD202X)


function getIntegralWYYKK(params::Dict)
    # this will give the integral of WYYKK of the equation 54 in 1 dimension
    @unpack oB, oWB, νCoord, LCoord, ΔCoord, l_n_max = params
    kernels = Array{Float64,4}(undef,LCoord,LCoord,l_n_max+1,l_n_max+1)
    modμ = nothing
    for l_n_field in 0:1:l_n_max
        for l_n_variable in 0:1:l_n_max
            for μ in 1:1:LCoord
                for μᶜ in 1:1:LCoord
                    #kernels[μᶜ,μ,l_n_variable+1,l_n_field+1],modμ=integralBsplineTaylorKernels1DWithWindow1D(oB,oWB,μᶜ,μ,νCoord,LCoord,ΔCoord,l_n_variable,l_n_field)
                    paramsBsplineTaylorIntegral1D=@strdict BsplineOrder = oB WBsplineOrder = oWB μᶜ = μᶜ μ = μ ν = νCoord L = LCoord Δ = ΔCoord l_n_variable = l_n_variable l_n_field = l_n_field
                    output=myProduceOrLoad(integralBsplineTaylorKernels1DWithWindow1D,paramsBsplineTaylorIntegral1D,"BsplineInt","BsplineTaylorIntegral1D")
                    kernels[μᶜ,μ,l_n_variable+1,l_n_field+1]=Num2Float64(output["kernels"])
                    modμ=output["modμ"]
                end
            end
        end
    end

    
    return @strdict(intKernelforνLΔ=kernels,modμ=modμ)
    # the target
    #integral1DWYYKK[iCoord][pointsIndices[linearμᶜ][iCoord],pointsIndices[linearμ][iCoord],l_n_variable,l_n_field]
end


function integralBsplineTaylorKernels1DWithWindow1D(params::Dict)
    @unpack BsplineOrder,WBsplineOrder,μᶜ,μ,ν,L,Δ,l_n_variable,l_n_field = params
    kernels,modμ=integralBsplineTaylorKernels1DWithWindow1D(BsplineOrder,WBsplineOrder,μᶜ,μ,ν,L,Δ,l_n_variable,l_n_field)
    return @strdict(kernels=kernels,modμ=modμ)
end


function integralBsplineTaylorKernels1DWithWindow1D(BsplineOrder,WBsplineOrder,μᶜ,μ,ν,L,Δ::Float64,l_n_variable,l_n_field)

    # Δ should be strictly Float64

    # this computes the analytical value of the 1D integral between B-spline fns and weighted Taylor kernels
    # \int dx Bspline Y_μᶜ Y_μ  K_{lᶜ-nᶜ}(y-y_μᶜ) K_{l-n}(y-y_μ)

    # unlike the previous integralBsplineTaylorKernels1D, it computes for a specific ν
    # Cˡη;μ are computed for a specific geometry, so even though the boundary condition reduce
    # the number of available points, each Taylor expansion for K_{l-n}(y-y_μ) should be Ok
    
    # however, the 'forgotten' μ (due to the whole) should be treated carefully 
    # (which I have not yet implemented here): I think Y_ignoredμ should be added to the Y_availableneighbouringμ

    # or maybe the 'forgotten' μ is anyways not available (and thus very probably not continuous)
    # so we just let this be forgotten 

        

    #@show μᶜ,μ,ν,L,Δ,l_n_variable,l_n_field

    kernelValue=0.0
   
    maximumOrder = maximum((BsplineOrder,WBsplineOrder,0)) # even if we only use the indicator functions (classical FD), we will call it for modμ (which we do not use for classical FD)
    params=@strdict maximumOrder numberNodes = L

    #output,_=@produce_or_load(BsplineTimesPolynomialsIntegrated,params,datadir("BsplineInt");filename = config -> savename("Bspline",params))
    
    output=myProduceOrLoad(BsplineTimesPolynomialsIntegrated,params,"BsplineInt","Bspline")
    #@show output
    nodeIndices,nodesSymbolic,b_deriv,integral_b,Δx,extFns,x,modμ =output["BsplineIntegraters"]


    if BsplineOrder=== -1
        # this is for an indicator function
        if l_n_variable === 0 && l_n_field === 0
            kernelValue=1.0
        else
            kernelValue=0.0
        end

    else
        
        
        # here we make a function Y_μ' Y_μ K_μ' K_μ (details ommitted)
        # note that ν is somewhere middle or at extremeties and 'ν+' expression is ommitted 

        Y_μᶜ = zeros(Num,L)
        Y_μ = zeros(Num,L)

        if WBsplineOrder === -1

            if μᶜ === ν
                Y_μᶜ = ones(Num,L) 
            end

            if μ === ν
                Y_μ = ones(Num,L)
            end

        else
            Y_μᶜ=b_deriv[:,μᶜ,1,WBsplineOrder+1]
            Y_μ =b_deriv[:,μ ,1,WBsplineOrder+1]
        end

        # modμ[3,:,ι+1] is the symbolic expression of the centre to compute Taylor kernels, which can be staggered!!!

        K_μᶜ=(x-Δx*modμ[3,μᶜ,WBsplineOrder+1])^l_n_variable
        K_μ =(x-Δx*modμ[3,μ,WBsplineOrder+1])^l_n_field

        # the convoluted function of all above
        F = mySimplify.(Y_μᶜ .* Y_μ .* K_μᶜ .* K_μ)


        # the target kernel integral

        targetKernel = integral_b[ν,BsplineOrder+1]

        dictionaryForSubstitute = Dict()
    
        for i in 0:1:BsplineOrder
            F = integrateTaylorPolynomials.(F,x) # integrate already for the 1st partial of W
            
            # mathematically I need to understand why but F cannot be disturbed by supplementary complexities due to constants 
            # that are arbitrarily put during the integral
            
            #F .-= substitute(F[ν],Dict(x=>nodesSymbolic[ν]))

            # F should be continuous, in general but since all the integral by parts
            # is done piecewisely we might not need this
            #=
            for iSegment in 2:1:L # this should be sequential
                lastValue = substitute(F[iSegment-1],Dict(x=>nodesSymbolic[iSegment]))
                startValue = substitute(F[iSegment],Dict(x=>nodesSymbolic[iSegment]))
                shiftValue = lastValue - startValue
                F[iSegment]=F[iSegment]+shiftValue
            end
            =#

            F .= mySimplify(F)

            for iSegment in nodeIndices
                dictionaryForSubstitute[extFns[1,iSegment,i+1]]=substitute(F[iSegment],Dict(x=>nodesSymbolic[iSegment]))
                dictionaryForSubstitute[extFns[2,iSegment,i+1]]=substitute(F[iSegment],Dict(x=>nodesSymbolic[iSegment+1]))
            end
        end

        #@show dictionaryForSubstitute,targetKernel
        

        kernelValue = substitute(targetKernel,dictionaryForSubstitute)  
        
        kernelValue = substitute(kernelValue,Dict(Δx=>Δ))/(BigInt(factorial(l_n_field))*BigInt(factorial(l_n_variable)))
        kernelValue = Num2Float64(kernelValue)
        
        #a= (Δ^(l_n_variable+l_n_field+1)-(-Δ)^(l_n_variable+l_n_field+1))/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        #@show a
    end
    return kernelValue,modμ
    
end



function integralBsplineTaylorKernels1D(BsplineOrder,Δ,l_n_variable,l_n_field)

    #call_new_fn(fn, args...) = Base.invokelatest(fn, args...)
    
    # this will compute \int dx Bspline K_{l-n} K_{lᶜ-nᶜ}
    middle_value = 0
    extreme_value = 0
    midPoint = BsplineOrder+2
    maxPoint = (BsplineOrder+1)*2 + 3
    nearboundaries_values=Array{Any,1}(undef,maxPoint)

    if BsplineOrder=== -1
        # this is for a delta function
        if l_n_variable === 0 && l_n_field === 0
            middle_value=1
            nearboundaries_values=ones(3)
        else
            middle_value=0
            nearboundaries_values=zeros(3)
        end

    elseif BsplineOrder >= 0
        maximumOrder = BsplineOrder
        params=@strdict maximumOrder numberNodes=L
        #output,_=@produce_or_load(BsplineTimesPolynomialsIntegrated,params,datadir("BsplineInt");filename = config -> savename("Bspline",params))
        output=myProduceOrLoad(BsplineTimesPolynomialsIntegrated,params,"BsplineInt","Bspline")
        numberNodes,integral_b_polys,N,Δx=output["BsplineIntegraters"]
        #fns=eval.(build_function.(integral_b_polys,N,Δx))
        middleNode = numberNodes ÷ 2
        #@show call_new_fn(fns[middleNode, BsplineOrder+1], l_n_variable + l_n_field + 1, Δ)
        #@show middle_value = call_new_fn(fns[middleNode,BsplineOrder+1](l_n_variable+l_n_field+1,Δ))/(factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        #@show integral_b_polys[middleNode, BsplineOrder+1]
     
        middle_value = Symbolics.substitute(integral_b_polys[middleNode, BsplineOrder+1],Dict(N=>l_n_variable+l_n_field,Δx => Δ))/(factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        #@show (Δ^(l_n_variable+l_n_field+1)-(-Δ)^(l_n_variable+l_n_field+1))/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        for iNode in 1:midPoint+2
            nearboundaries_values[iNode] = Symbolics.substitute(integral_b_polys[iNode, BsplineOrder+1],Dict(N=>l_n_variable+l_n_field,Δx => Δ))/(factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
            nearboundaries_values[maxPoint-iNode+1] = Symbolics.substitute(integral_b_polys[numberNodes-iNode+1, BsplineOrder+1],Dict(N=>l_n_variable+l_n_field,Δx => Δ))/(factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        end
    end

#=
    elseif BsplineOrder === 0
        #
    elseif BsplineOrder=== 1
        middle_value = (Δ^(l_n_variable+l_n_field+1)-(-Δ)^(l_n_variable+l_n_field+1))/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(BigInt(l_n_variable))*factorial(BigInt(l_n_field)))
        #extreme_value = (Δ^{l_n_variable+l_n_field+1})/((l_n_variable+l_n_field+2)*(l_n_variable+l_n_field+1)*factorial(l_n_variable)*factorial(l_n_field))
    end
    =#


    return middle_value,nearboundaries_values
end
