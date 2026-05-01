using KernelAbstractions


function TaylorOptions(itplParams,supplementaryOrder)
    options=(YorderBspace=itplParams.YorderBspace,YorderBtime=itplParams.YorderBtime,supplementaryOrder=supplementaryOrder,pointsμInSpace=itplParams.ptsSpace,pointsμInTime=itplParams.ptsTime,offsetμInΔyInSpace=itplParams.offsetSpace,offsetμInΔyInTime=itplParams.offsetTime)
    return options
end

function makeOPTsemiSymbolic(params::Dict)
    @unpack famousEquationType, Δ, orderBtime, orderBspace, pointsInSpace, pointsInTime, supplementaryOrder, fieldItpl, materItpl = params
    Δnum = SVector(Δ)
    # construction of NamedTuples
    trialFunctionsCharacteristics=(orderBtime=orderBtime,orderBspace=orderBspace,pointsInSpace=pointsInSpace,pointsInTime=pointsInTime)

    # here we can compute the different interpolated Taylor expansion options
    TaylorOptionsμ=TaylorOptions(fieldItpl,supplementaryOrder)
    TaylorOptionsμᶜ=TaylorOptions(materItpl,supplementaryOrder)

    equationCharacteristics,equationCharacteristicsForce=famousEquations(famousEquationType)

    # compact coefficients for l.h.s. of the equation
    equationCharacteristics=equationCharacteristics
    numbersOfTheSystem=numbersOfTheExpression(equationCharacteristics,trialFunctionsCharacteristics,TaylorOptionsμ,TaylorOptionsμᶜ)
    _,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ=investigateDependencies(equationCharacteristics,numbersOfTheSystem,trialFunctionsCharacteristics,TaylorOptionsμ,TaylorOptionsμᶜ)
    bigα,varM=bigαFinder(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ)
    Ajiννᶜ,AjiννᶜU,Ulocal=constructAmatrix(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ,Δnum,bigα,varM)
    lhs=(Ajiννᶜ=Ajiννᶜ,AjiννᶜU=AjiννᶜU,Ulocal=Ulocal,varM=varM)

    # compact coefficients for r.h.s. of the equation
    equationCharacteristics=equationCharacteristicsForce
    numbersOfTheSystem=numbersOfTheExpression(equationCharacteristics,trialFunctionsCharacteristics,TaylorOptionsμ,TaylorOptionsμᶜ)
    _,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ=investigateDependencies(equationCharacteristics,numbersOfTheSystem,trialFunctionsCharacteristics,TaylorOptionsμ,TaylorOptionsμᶜ)
    bigα,varM=bigαFinder(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ)
    Γjiννᶜ,ΓjiννᶜF,Flocal =constructAmatrix(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ,Δnum,bigα,varM)
    rhs=(Γjiννᶜ=Γjiννᶜ,ΓjiννᶜF=ΓjiννᶜF,Flocal=Flocal,varF=varM)

    #
    nodes=configsTaylorμ.availablePointsConfigurations
    centresIndices=configsTaylorμ.centrePointConfigurations
    timeMarching=numbersOfTheSystem.timeMarching
    recette=(lhs=lhs,rhs=rhs,nodes=nodes,centresIndices=centresIndices,numbersOfTheSystem=numbersOfTheSystem,timeMarching=timeMarching)
    return @strdict(recette)

end


function constructAmatrix(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,Δnum,bigα,varM)
    return constructAmatrix(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμ,configsTaylorμ,Δnum,bigα,varM)
end

function constructAmatrix(equationCharacteristics,numbersOfTheSystem,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ,Δnum,bigα,varM;ImakeReport=true)
    
    # for the future develpments: ν can move but it's already more or less coded! look at pointν and nGeometry

    # This function is for one iConfigGeometry

    #region preparation
    @unpack fields,vars = equationCharacteristics
    @unpack nCoordinates,NtypeofExpr,NtypeofExpr,NtypeofFields = numbersOfTheSystem
    @unpack multiOrdersIndices,availablePointsConfigurations,centrePointConfigurations,availableμPoints,availableμaxes, numberGeometries = configsTaylorμ
    availableμᶜPoints = configsTaylorμᶜ.availableμPoints
    availableμᶜaxes = configsTaylorμᶜ.availableμaxes
    orderBspline = ordersForSplinesμ.orderBspline
    YorderBsplineμ = ordersForSplinesμ.YorderBspline
    YorderBsplineμᶜ= ordersForSplinesμᶜ.YorderBspline
    #nCoordinates = nCoordinates isa Val ? nCoordinates.val : nCoordinates

    #for iConfigGeometry ∈ 1:numberGeometries
    nGeometry = numberGeometries
    nGeometry =1
    iConfigGeometry = 1

    @show pointsIndices=availablePointsConfigurations[iConfigGeometry] # CartesianIndices
    @show middleLinearν=centrePointConfigurations[iConfigGeometry] # scalar
    @show μPoints = availableμPoints[iConfigGeometry] # Array(SVector)
    @show μᶜPoints = availableμᶜPoints[iConfigGeometry]
    @show μaxes = availableμaxes[iConfigGeometry]
    @show μᶜaxes = availableμᶜaxes[iConfigGeometry]
    @show size(μPoints)
    @show pointν = pointsIndices[middleLinearν] # SVector
    
    # this is fully GPU optimised version of ASymbolic 
    
    nPoints = length(pointsIndices)
    nLs = length(multiOrdersIndices)

    # 

    L_MINUS_N = multiOrdersIndices
    L_MINUS_N = L_MINUS_N .-L_MINUS_N[1] 
    # here L_MINUS_N is truly \mathbf{l}-\mahtbf{n} ∈ \mathbb{Z}_{≥0}

    #endregion

    #region we compute the integral of WYYKK in 1D domain 
    
    # look at the debug1DKernelIntegral.ipynb    
    
    coefWYYKK = Array{Any,1}(undef,nCoordinates) # Array of (l_n_max+1,lᶜ_nᶜ_max+1,length(μs),length(μᶜs),length(ν) ) times nCoordinates
    for iCoord ∈ 1:nCoordinates
        maxNodes = pointsIndices[end][iCoord]
        nodesFromOne = [1,2,3] # ∈ Z like [1,2,3], an array of integers collect(1:1:Npoints) (nothing else!!)
        ν = (pointν[iCoord]) # this should be one point (for the moment)
        lᶜ_nᶜ_max = L_MINUS_N[end][iCoord] # variable
        l_n_max = L_MINUS_N[end][iCoord] # field
        params = (orderBspline1D=orderBspline[iCoord], YorderBspline1Dμᶜ=YorderBsplineμᶜ[iCoord], YorderBspline1Dμ=YorderBsplineμ[iCoord], μᶜs=μᶜaxes[iCoord], μs=μaxes[iCoord], maxNode = pointsIndices[end][iCoord], ν =(pointν[iCoord]), lᶜ_nᶜ_max=lᶜ_nᶜ_max, l_n_max=l_n_max,  Δ=Δnum[iCoord],ImakeReport=ImakeReport)
        coefWYYKK[iCoord] = WYYKKIntegralNumerical(params) 
    end

    #endregion

    #region get CˡηGlobal (for ν)

    @show typeof(μPoints), μPoints[1], typeof(pointsIndices)

    coefInversionDict = @strdict multiOrdersIndices pointsIndices μpointsIndices=μPoints Δ=Δnum
    output=myProduceOrLoad(TaylorCoefInversion,coefInversionDict,"taylorCoefInv")
    Cˡη=output["CˡηGlobal"]

    coefInversionDict = @strdict multiOrdersIndices pointsIndices μpointsIndices=μᶜPoints Δ=Δnum
    output=myProduceOrLoad(TaylorCoefInversion,coefInversionDict,"taylorCoefInv")
    Cˡηᶜ=output["CˡηGlobal"]

    #endregion 

    #region 
    @show nTotalSmallα = sum(length(bigα[iExpr, iField]) for iField ∈ 1:NtypeofFields, iExpr ∈ 1:NtypeofExpr)

    # but this will already include bigα so the coefficients for each α_{nn'ji} should be given here
    #endregion


    #region useful LinearIndices conversion functions
    
    LI_points = LinearIndices(pointsIndices)
    LI_L_MINUS_N_plus_1 = LinearIndices(L_MINUS_N.+vec2car(ones(Int,nCoordinates)))
  
    #endregion


    #region make the table for each (x',x,n',n) (x=η+μ)
    
    tableForLoop = Array{Int32,3}(undef,2+nCoordinates*2,nLs*nLs,nTotalSmallα)

    fill!(tableForLoop, 0)
    indexLinearα = 1
    for iExpr ∈ 1:NtypeofExpr,iField ∈ 1:NtypeofFields
        α = bigα[iExpr,iField]
        for eachα ∈ α
            nᶜ = eachα.nᶜ - vec2car(ones(Int,nCoordinates))
            n = eachα.n - vec2car(ones(Int,nCoordinates))
            #nodeValue = eachα.node # not important at this point
            # Available indices
            Lᶜ_avail = (nᶜ .+ L_MINUS_N) ∩ L_MINUS_N
            L_avail = (n .+ L_MINUS_N) ∩ L_MINUS_N
            Lᶜ_Nᶜ_avail = Lᶜ_avail .- nᶜ
            L_N_avail = L_avail .- n
            
            iL=1
            
            for l ∈ L_avail, lᶜ ∈ Lᶜ_avail
                
                tableForLoop[1,iL,indexLinearα]= LI_L_MINUS_N_plus_1[lᶜ+vec2car(ones(Int,nCoordinates))]
                tableForLoop[2,iL,indexLinearα]= LI_L_MINUS_N_plus_1[l+vec2car(ones(Int,nCoordinates))]
                #tableForLoop[3,iL,indexLinearα]= LI_L_MINUS_N_plus_1[lᶜ-nᶜ+vec2car(ones(Int,nCoordinates))]
                #tableForLoop[4,iL,indexLinearα]= LI_L_MINUS_N_plus_1[l-n+vec2car(ones(Int,nCoordinates))]
                tmplᶜ_nᶜ = lᶜ-nᶜ+vec2car(ones(Int,nCoordinates))
                tmpl_n = l-n+vec2car(ones(Int,nCoordinates))
                for iCoord ∈ 1:nCoordinates
                    tableForLoop[2+iCoord,iL,indexLinearα] = tmplᶜ_nᶜ[iCoord]
                    tableForLoop[2+iCoord+nCoordinates,iL,indexLinearα] = tmpl_n[iCoord]
                end

                iL += 1
            end
       
            indexLinearα += 1
        end
    end

    #endregion


    #region make a dictionary for μ ∈ μPoints and its linearised version

    tableForμPoints = Array{Int32,2}(undef,nCoordinates,length(μPoints))
    for iμ ∈ CartesianIndices(μPoints), iCoord ∈ 1:nCoordinates
            tableForμPoints[iCoord,iμ]=iμ[iCoord]
    end

    tableForμᶜPoints = Array{Int32,2}(undef,nCoordinates,length(μᶜPoints))
    for iμᶜ ∈ CartesianIndices(μᶜPoints), iCoord ∈ 1:nCoordinates        
        tableForμᶜPoints[iCoord,iμᶜ]=iμᶜ[iCoord]
    end
    


    #endregion


    
    #region adapt the arrays to the GPU backend
    tableForμPoints_gpu = Adapt.adapt(backend,tableForμPoints)
    tableForμᶜPoints_gpu = Adapt.adapt(backend,tableForμᶜPoints)
    tableForLoop_gpu = Adapt.adapt(backend,tableForLoop)
    C_gpu = Adapt.adapt(backend,Float32.(Cˡη))
    Cᶜ_gpu = Adapt.adapt(backend,Float32.(Cˡηᶜ))

    #endregion

    #region preparation for GPU computation
    # Collect the size of each array
    @show all_sizes = collect.(size.(coefWYYKK))  # vector of tuples

    # Element-wise maximum across all dimensions
    max_size = map((xs...) -> maximum(xs), all_sizes...)

    int_total_float32 = Array{Float32,6}(undef, nCoordinates, max_size...)
    int_total_float32 .= 0.f0

    for iCoord ∈ 1:nCoordinates
        small_size = size(coefWYYKK[iCoord])
        tmpMatrix = Float32.(coefWYYKK[iCoord])
        tmpRange = CartesianIndices(tmpMatrix)
        int_total_float32[iCoord,tmpRange] = tmpMatrix
    end
    int_gpu = Adapt.adapt(backend,int_total_float32)

    output_gpu = Adapt.adapt(backend, zeros(Float32, nPoints, nPoints, nTotalSmallα, nGeometry))

    #endregion



    #region launch GPU computation

    # scalars in Int32
    P      = Int32(nPoints)
    Pμᶜ    = Int32(length(μᶜPoints))
    Pμ     = Int32(length(μPoints))
    L      = Int32(nLs)
    nDim   = Int32(nCoordinates)
    nα     = Int32(nTotalSmallα)
    nGeometry_gpu=Int32(nGeometry)
    int1max = Int32(max_size[1])
    int2max = Int32(max_size[2])

    @show typeof(output_gpu)         # must be MtlArray{Float32,3}
    #@show typeof(C_gpu)              # must be MtlArray{Float32,3}
    #@show typeof(int_gpu)            # must be MtlArray{Float32,5}
    #@show typeof(tableForLoop_gpu)   # must be MtlArray{Int32,3}
    #@show typeof(tableForPoints_gpu) # must be MtlMatrix{Int32}
    #@show typeof(P) typeof(L) typeof(nDim) typeof(nα) typeof(int1max) typeof(int2max)

    kernel! = windowContraction!(backend,(8,8,8))#,128,size(output_gpu))
    kernel!(output_gpu, C_gpu, Cᶜ_gpu, int_gpu, tableForLoop_gpu,tableForμPoints_gpu,tableForμᶜPoints_gpu,
       P,Pμᶜ,Pμ,L,nDim,nα,int1max,int2max,nGeometry_gpu,ndrange=(P,P,nα*nGeometry))
    KernelAbstractions.synchronize(backend)

    # Here output_gpu[x',x,eachα] ∑μ' ∑μ ∑l' ∑ l C[x',l',μ'] C[x,l,μ] ∏_iCoord K[iCoord][μ',μ,l'-n',l-n]
    # with (n',n) depends on eachα and x=η+μ0

    @show "GPU computation of Ajiννᶜ: done"
    newCoef = Adapt.adapt(CPU(), output_gpu)
    #endregion


    iConfigGeometry = 1

    #region contruct Ulocal

    # the order is: (νᶜ,) ν, i, j  here

    Ulocal = Array{Num,2}(undef,length(pointsIndices),NtypeofFields)
    for iField in eachindex(fields)
        newstring=split(string(fields[iField]),"(")[1]
        Ulocal[:,iField]=Symbolics.variables(Symbol(newstring),1:length(pointsIndices))
    end

    #endregion

    #region make Ajiννᶜ and AjiννᶜU symbolically (which we will soon remove!)
    Ajiννᶜ = Array{Num,4}(undef,length(pointsIndices),NtypeofFields,NtypeofExpr,nGeometry)
    AjiννᶜU = Array{Num,2}(undef,NtypeofExpr,nGeometry)
    
    # this is the cost function for ν point so the number of elements is just the number of expressions (governing equations)
    Ajiννᶜ .= 0
    AjiννᶜU .= 0
    indexLinearα = 1
   
    for iExpr ∈ 1:NtypeofExpr,iField ∈ 1:NtypeofFields
        α = bigα[iExpr,iField]
        for eachα ∈ α
            @show nodeValue=eachα.node
            for x ∈ pointsIndices
                xLinear = LI_points[svec2car(x)]
                
                localmapηᶜ=Dict()
                for iVar ∈ eachindex(vars)
                    localmapηᶜ[vars[iVar]]=varM[iVar,xLinear][]
                end
                for xᶜ ∈ pointsIndices
                    xᶜLinear = LI_points[svec2car(xᶜ)]
                    U_HERE = Ulocal[xᶜLinear,iField]                    
                    substitutedValue = substitute(nodeValue, localmapηᶜ)
                    Ajiννᶜ[xᶜLinear,iField,iExpr] += newCoef[xLinear,xᶜLinear,indexLinearα,iConfigGeometry] *substitutedValue
                    AjiννᶜU[iExpr]+= Ajiννᶜ[xᶜLinear,iField,iExpr] * U_HERE
                end

            end
            indexLinearα += 1
        end

    end

    #return coefWYYKK,Cˡη,tableForμPoints, newCoef
    return Ajiννᶜ,AjiννᶜU,Ulocal
end






@kernel function windowContraction!(
    output::AbstractArray{Float32,4},      # (P, P, nα, nGeometry)
    C::AbstractArray{Float32,3},           # (P, L, Pμ)
    Cc::AbstractArray{Float32,3},          # (P, L, Pμᶜ)
    int::AbstractArray{Float32,6},         # (nDim, int1max, int1max, int2max, int2max, nGeometry)
    table::AbstractArray{Int32,3},         # (2 + 2*nDim, nLoop, nα)
    tableμPoints::AbstractArray{Int32,2},  # (nDim, Pμ)
    tableμᶜPoints::AbstractArray{Int32,2}, # (nDim, Pμᶜ)
    P::Int32, Pμᶜ::Int32, Pμ::Int32, L::Int32,
    nDim::Int32, nα::Int32, int1max::Int32, int2max::Int32, nGeometry::Int32
)
    (xᶜ, x, ag) = @index(Global, NTuple)

    if xᶜ <= P && x <= P && ag <= nα * nGeometry
        α = Int32(((ag - 1) % nα) + 1)
        iGeometry = Int32(((ag - 1) ÷ nα) + 1)

        acc = 0.0f0

        @inbounds for idx in 1:size(table, 2)
            lᶜ = table[1, idx, α]
            l  = table[2, idx, α]

            if lᶜ > 0 && l > 0
                for μᶜ in 1:Pμᶜ
                    for μ in 1:Pμ
                        prod_int = 1.0f0

                        for iDim in 1:nDim
                            k = table[2 + iDim, idx, α]
                            lp = table[2 + nDim + iDim, idx, α]

                            if k > 0 && lp > 0
                                mᶜ = tableμᶜPoints[iDim, μᶜ]
                                m  = tableμPoints[iDim, μ]

                                prod_int *= int[iDim, k, lp, m, mᶜ, iGeometry]
                            else
                                prod_int = 0.0f0
                                break
                            end
                        end

                        acc += Cc[xᶜ, lᶜ, μᶜ] * C[x, l, μ] * prod_int
                    end
                end
            end
        end

        output[xᶜ, x, α, iGeometry] = acc
    end
end
