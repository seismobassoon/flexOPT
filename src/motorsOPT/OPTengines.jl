# this is the GPU version of OPTEngines/OPTnewEngines.jl 
# December 2025 Nobuaki Fuji (IPGP/UPC/IUF)



function OPTobj_temporarily_obsolete(operatorConfigurations::Dict)
    # this is just a wrapper for the OPTobj function below, for DrWatson package
    @unpack myEquationInside, Δnum, orderBtime, orderBspace, WorderBtime,WorderBspace,supplementaryOrder,pointsInSpace, pointsInTime,IneedExternalSources, iExperiment= operatorConfigurations

    TaylorOptions=(WorderBtime=WorderBtime,WorderBspace=WorderBspace,supplementaryOrder=supplementaryOrder)
    trialFunctionsCharacteristics=(orderBtime=orderBtime,orderBspace=orderBspace,pointsInSpace=pointsInSpace,pointsInTime=pointsInTime)
    @time operatorData=OPTobj(myEquationInside,Δnum; trialFunctionsCharacteristics=trialFunctionsCharacteristics,TaylorOptions=TaylorOptions,iExperiment=iExperiment)
    #AjiννᶜU=operatorData[1]
    #utilities=operatorData[2]   

    operatorForceData=nothing
    # if you do not want to apply external forces, it is possible to skip below
    if IneedExternalSources 
        @time operatorForceData=OPTobj(myEquationInside,Δnum; trialFunctionsCharacteristics=trialFunctionsCharacteristics,TaylorOptions=TaylorOptions,iExperiment=iExperiment)  
        #@show Γg = operatorForceData[1]
        #utilitiesForce = operatorForceData[2]
    end
    eqInfo=(exprs=exprs,fields=fields,vars=vars,extexprs=extexprs,extfields=extfields,extvars=extvars,coordinates=coordinates)
    operators=(operatorPDE=operatorData, operatorForce=operatorForceData,eqInfo=eqInfo)
    return @strdict(operators)
end



function OPTobj(OPTconfig::Dict)
    @unpack myEquationInside, Δnum, TaylorOptions, trialFunctionsCharacteristics = OPTconfig

    Ajiννᶜs,AjiννᶜUs,Ulocals,utilities=OPTobj(myEquationInside.exprs,myEquationInside.fields,myEquationInside.vars,myEquationInside.coordinates,myEquationInside.∂,Δnum;TaylorOptions=TaylorOptions, trialFunctionsCharacteristics=trialFunctionsCharacteristics)
    
    Γjkννᶜs,ΓjkννᶜFs,Flocals= nothing,nothing,nothing
    if myEquationInside.extvars !== nothing
        Γjkννᶜs,ΓjkννᶜFs,Flocals,_ = OPTobj(myEquationInside.extexprs,myEquationInside.extfields,myEquationInside.extvars,myEquationInside.coordinates,myEquationInside.∂,Δnum;TaylorOptions=TaylorOptions, trialFunctionsCharacteristics=trialFunctionsCharacteristics)
    end
    output = (Ajiννᶜs=Ajiννᶜs,AjiννᶜUs=AjiννᶜUs,Ulocals=Ulocals,Γjkννᶜs=Γjkννᶜs,ΓjkννᶜFs=ΓjkννᶜFs,Flocals=Flocals,utilities=utilities)
    return @strdict(output)
end


function OPTobj(exprs,fields,vars,coordinates,∂,Δnum;TaylorOptions=(WorderBtime=1,WorderBspace=1,supplementaryOrder=2), trialFunctionsCharacteristics=(orderBtime=1,orderBspace=1, pointsInSpace=2,pointsInTime=2))

    #region General introduction, some cautions

    # Nobuaki Fuji @ IPGP/UPC/IUF since 2024
    #
    # thanks to CLEEDIs and hackathons in black forest, the function is getting matured, with GPU vectorisation 
    #
    #
    # encouraged by Thibault Duretz @ U. Frankfurt Goethe, Kurama Okubo @ NIED
    #
    # with some inputs from Giacomo Aloisi @ ETH during Julia hackathon 
    #                                  in the black forest October 2024
    #
    # intermediate presentations: IPGP-CIA workshop October 2024; IPGP-ERI workshop November 2024; lighttalk @ systemI December 2024
    #

    # This function SHOULD offer an OPT operator for a given expression 
    # and the boundary conditions on the 8 sides of a local 4D Cartesian box 

     # B-spline order (-1: Dirac's delta function (finite difference); 0: step function (2 points); 1: triangle (3 points); 2+: beautiful curve) 
    
    # this should be corrected according to the orderPDEs that should be given by the gvn. eqn.
    # orderBtime and orderBspace should be given


    # This function should expand the compact PDE given by a textbook or by an imagination
    # in order to organise an OPT operators for 4D Cartesian coordinates. 

    # orderBtime and orderBspace are the order of B-spline (if -1: delta function) to be used as test functions
    # as the paper in prep. (as of August 2024) shows, the numbers of points used (pointsUsed) are not necessarily related to orderBtime and orderBspace

    # even for spherical harmonic version with only z-dependency, this should work (but the expansion should be done a posteriori)

    # Note that variableDependency defines the (maximum) dimension of physics to be solved. 



    # some notes on March 2025

    # highestOrderPartial should be somehow determined more naturally!!
    # (l-n) max is the same as pointsUsedForFields here!

    # if timeMarching === true then we consider that the last coordinate is time
    
    # CˡηSymbolicInversion is highly recommended to be false since it takes really a big effort for nothing
    # pointsUsed(ForFields) should be the highest order of PDE + 1 at least




    #endregion

    #region initialising, unpacking etc. 


    myEquationOneSide = (exprs=exprs,fields=fields,vars,coordinates=coordinates)

    timeMarching = any(a -> a === timeDimensionString, string.(coordinates))


    @unpack orderBtime, orderBspace, pointsInSpace, pointsInTime = trialFunctionsCharacteristics
    @unpack WorderBtime, WorderBspace,supplementaryOrder = TaylorOptions

    NtypeofExpr=length(exprs)   # number of governing equations
    NtypeofMaterialVariables=length(vars) # number of material coefficients
    NtypeofFields=length(fields) # number of unknown fields
    
    Ndimension = length(coordinates) # we do not change this for the moment, especially for the time-marching scheme
    pointsUsed = ones(Int, Ndimension).*(pointsInSpace+1)
    if timeMarching
        pointsUsed[end]=pointsInTime+1
    end


    if length(Δnum) !== Ndimension 
        @error "the numerical delta increment has not the same dimension!"
    end

 


    #endregion

    #region investigation of all the fields and vars dependencies in terms of x-y-z-t

    variableDependency=ones(Int,Ndimension)
    fieldDependency=ones(Int,Ndimension)
    eachVariableDependency=ones(Int,Ndimension,NtypeofMaterialVariables) 
    eachFieldDependency=ones(Int,Ndimension,NtypeofFields)
  
    for iFields in 1:NtypeofFields
        eachFieldDependency[:,iFields]=findCartesianDependency(fields[iFields],Ndimension,∂)
        fieldDependency = fieldDependency .* (ones(Int,Ndimension).-eachFieldDependency[:,iFields])
    end


    for iVars in 1:NtypeofMaterialVariables
        eachVariableDependency[:,iVars]=findCartesianDependency(vars[iVars],Ndimension,∂)
        variableDependency = variableDependency .* (ones(Int,Ndimension).-eachVariableDependency[:,iVars])
    end

    

    fieldDependency = ones(Int,Ndimension).-fieldDependency
    variableDependency = ones(Int,Ndimension).-variableDependency

    # here we correct variableDependency with fieldDependency: if fieldDependency is zero then we do not take care of that dimension for the variables
    variableDependency = variableDependency .* fieldDependency

    #endregion

    #region definition of points in time and space to be used

    # heaviside(x) = x > 0 ? 1 : x == 0 ? 0 : -1

    # the orders of B-spline functions, depending on fields 

    orderBspline=zeros(Int,Ndimension)
    WorderBspline=zeros(Int,Ndimension)

    if timeMarching
        orderBspline[Ndimension]=orderBtime*fieldDependency[Ndimension]
        orderBspline[1:Ndimension-1]=orderBspace*fieldDependency[1:Ndimension-1]
        WorderBspline[Ndimension]=WorderBtime*fieldDependency[Ndimension]
        WorderBspline[1:Ndimension-1]=WorderBspace*fieldDependency[1:Ndimension-1]
    else
        orderBspline[1:Ndimension]=orderBspace*fieldDependency[1:Ndimension]
        WorderBspline[1:Ndimension]=WorderBspace*fieldDependency[1:Ndimension]
    end
    
    # the maximum number of points used in the vicinity of the node, which is independent of the order of B-spline functions (see our paper)
    pointsUsedForFields=(pointsUsed.-1).*fieldDependency.+1

    # orderExpressions is the maximal orders of partials that we could expect in the expressions
    orderExpressions=pointsUsedForFields
    
    # numbers of points to evaluate the integral for the governing equation filtered by the test functions
    
    # orderU is the maximum orders for the fields that we will use for OPT coefficients' exploration
    orderU = (orderExpressions .-1) .+ (supplementaryOrder .*fieldDependency).+1 
    # we restore this orderU since we need to control this 

    #endregion

    #region analysis of expressions to obtain the α_{n'nji}

    bigα=Array{Any,2}(missing,NtypeofFields,NtypeofExpr)
    varM=nothing
    for iExpr in eachindex(exprs)
        for iField in eachindex(fields)
            
            tmpNonZeroAlphas=PDECoefFinder(orderExpressions,coordinates,exprs[iExpr],fields[iField],vars) 
            # we assume that the pointsUsedForFields represent the highest order of partials
            bigα[iField,iExpr]=unique(tmpNonZeroAlphas)
        end
    end
    varM=varMmaker(pointsUsedForFields,coordinates,vars,∂)
    @show bigα,varM
    #endregion

    #region Preparation for Taylor expansion
    
    orderTaylors=Array{Any,Ndimension}(undef,Tuple(orderU))
    pointsInSpaceTime=Array{Any,Ndimension}(undef,Tuple(pointsUsedForFields))
    
    
    multiOrdersIndices=CartesianIndices(orderTaylors)

    #availablePointsConfigurations = []
    #centrePointConfigurations=[]

    availablePointsConfigurations = Array{Array{Vector{Int64},Ndimension},1}()
    centrePointConfigurations=Array{Int64,1}()

    #endregion

    #region Cartesian indices that can be available to use (normally: iGeometry=1)

    multiPointsIndices=CartesianIndices(pointsInSpaceTime)
    # this is the whole local Cartesian grids (without any lacking points)
    
    tmpVecForMiddlePoint = ((car2vec(multiPointsIndices[end]).-1 ).÷2 ).+1 # only valid for testOnlyCentre
    midTimeCoord = nothing
    if timeMarching
        midTimeCoord=car2vec(multiPointsIndices[end])[end]-1
        tmpVecForMiddlePoint[end]=midTimeCoord
        #AjiννᶜU = Array{Num,2}(undef,length(multiPointsIndices)÷(midTimeCoord+1),NtypeofExpr)
    end
    #@show tmpVecForMiddlePoint 
    middleν=vec2car(tmpVecForMiddlePoint)


    availablePointsConfigurations=push!(availablePointsConfigurations,car2vec.(multiPointsIndices))
    centrePointConfigurations=push!(centrePointConfigurations,LinearIndices(multiPointsIndices)[middleν])

    #endregion

    #region obtaining the semi-symbolic expression of cost function based on eqns. 52 and 53.

    # before calling AuSymbolic we can manipulate pointsIndices for various boundary configurations
    Ajiννᶜs=[]
    AjiννᶜUs=[]
    Ulocals=[]


    #AmatrixSemiSymbolicGPU(coordinates,multiOrdersIndices,pointsIndices,multiPointsIndices,middleLinearν,Δnum,varM,bigα,orderBspline,WorderBspline,NtypeofExpr,NtypeofFields)


    for iConfigGeometry in eachindex(availablePointsConfigurations) 
        pointsIndices=availablePointsConfigurations[iConfigGeometry]
        middleLinearν=centrePointConfigurations[iConfigGeometry]
        #varM is given above for the max number of points used 
        #tmpAjiννᶜU,tmpUlocal=AuSymbolic(coordinates,multiOrdersIndices,pointsIndices,multiPointsIndices,middleLinearν,Δ,varM,bigα,orderBspline,WorderBspline,NtypeofExpr,NtypeofFields)
        coefsSemiSymbolic=AmatrixSemiSymbolicGPU(myEquationOneSide,multiOrdersIndices,pointsIndices,multiPointsIndices,middleLinearν,Δnum,varM,bigα,orderBspline,WorderBspline,NtypeofExpr,NtypeofFields)
        Ajiννᶜs=push!(Ajiννᶜs,coefsSemiSymbolic.Ajiννᶜ)
        AjiννᶜUs=push!(AjiννᶜUs,coefsSemiSymbolic.AjiννᶜU)
        Ulocals=push!(Ulocals,coefsSemiSymbolic.Ulocal)
    end


    #endregion

    #region outputs
    
    utilities=(middlepoint=middleν,middlepointLinear=centrePointConfigurations[1],localPointsIndices=multiPointsIndices,localMaterials=varM,localFields=Ulocals[1],timeMarching=timeMarching)
   
    return Ajiννᶜs,AjiννᶜUs,Ulocals,utilities
    

    #endregion
    
end


function AmatrixSemiSymbolicGPU(myEquationInside,multiOrdersIndices,pointsIndices,multiPointsIndices,middleLinearν,Δ,varM,bigα,orderBspline,WorderBspline,NtypeofExpr,NtypeofFields;threads = 256,backend=backend)

    #region preparation

    exprs=myEquationInside.exprs
    fields=myEquationInside.fields
    vars=myEquationInside.vars
    coordinates=myEquationInside.coordinates

    # this is fully GPU optimised version of ASymbolic 
    nCoordinates = length(coordinates)
    nPoints = length(pointsIndices)
    nLs = length(multiOrdersIndices)
    # normally we put only the centre points (in the present time)
    # but the boundary condition requires some points a bit more challenging
    νIndices=[pointsIndices[middleLinearν]]
    # νIndices can be pointsIndices if one needs to look for
    # boundary operators

    L_MINUS_N = multiOrdersIndices
    L_MINUS_N = L_MINUS_N .-L_MINUS_N[1]
    # here L_MINUS_N is truly \mathbf{l}-\mahtbf{n} ∈ \mathbb{Z}_{≥0}

    #endregion

    #region we compute the integral of WYYKK in 1D domain (μ ∈ \mathbb{Z}/2 has not yet been implemented)

    integral1DWYYKK = Array{Any,1}(undef,nCoordinates)
    modifiedμ=Array{Any,1}(undef,nCoordinates)
    for iCoord in eachindex(coordinates) # for each 
        integralParams = @strdict oB =orderBspline[iCoord] oWB = WorderBspline[iCoord] νCoord=pointsIndices[middleLinearν][iCoord] LCoord = multiPointsIndices[end][iCoord] ΔCoord=Δ[iCoord] l_n_max=L_MINUS_N[end][iCoord]
        output = myProduceOrLoad(getIntegralWYYKK,integralParams,"intKernel")
        integral1DWYYKK[iCoord] = output["intKernelforνLΔ"]
        modifiedμ[iCoord] = output["modμ"] # this can be still 'nothing'
    end



    #endregion

    #region get CˡηGlobal (for ν)

    coefInversionDict = @strdict coordinates multiOrdersIndices pointsIndices Δ WorderBspline modifiedμ

    #@show coordinates multiOrdersIndices pointsIndices Δ WorderBspline modifiedμ

    output=myProduceOrLoad(TaylorCoefInversion,coefInversionDict,"taylorCoefInv")
    Cˡη=output["CˡηGlobal"]

    #endregion 

    #region objectives (Ajiννᶜ)

    # the order is: (νᶜ,) ν, i, j  here
    Ajiννᶜ = Array{Num,4}(undef,nPoints,nPoints,NtypeofFields,NtypeofExpr)

    # but this will already include bigα so the coefficients for each α_{nn'ji} should be given here

    #endregion

    #region useful LinearIndices conversion functions
    
    LI_points = LinearIndices(pointsIndices)
    LI_L_MINUS_N_plus_1 = LinearIndices(L_MINUS_N.+vec2car(ones(Int,nCoordinates)))
  
    #endregion

    #region make the table for each (x',x,n',n) (x=η+μ)
    nTotalSmallα = sum(length(bigα[iExpr, iField]) for iExpr in 1:NtypeofExpr, iField in 1:NtypeofFields)

    tableForLoop = Array{Int32,3}(undef,2+nCoordinates*2,nLs*nLs,nTotalSmallα) # for every α it will give the l and l-n
    
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
            for lᶜ ∈ Lᶜ_avail, l ∈ L_avail
                
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

    #region make a dictionary for μ ∈ pointsIndices and its linearised version

    tableForPoints = Array{Int32,2}(undef,nCoordinates,nPoints)

    for μ ∈ pointsIndices, iCoord ∈ 1:nCoordinates
        linearisedμ = LI_points[vec2car(μ)]
        tableForPoints[iCoord,linearisedμ]=vec2car(μ)[iCoord]
    end

    #endregion

    #region adapt the arrays to the GPU backend
    tableForPoints_gpu = makeGPUarray(backend,tableForPoints)
    tableForLoop_gpu = makeGPUarray(backend,tableForLoop)
    C_gpu = makeGPUarray(backend,Float32.(Cˡη))


    # Collect the size of each array
    all_sizes = collect.(size.(integral1DWYYKK))  # vector of tuples

    # Element-wise maximum across all dimensions
    max_size = map((xs...) -> maximum(xs), all_sizes...)

    int_total_float32 = Array{Float32,5}(undef, nCoordinates, max_size...)
    int_total_float32 .= 0.f0
    for iCoord in 1:nCoordinates
        small_size = size(integral1DWYYKK[iCoord])
        tmpMatrix = Float32.(integral1DWYYKK[iCoord])
        tmpRange = CartesianIndices(tmpMatrix)
        int_total_float32[iCoord,tmpRange] = tmpMatrix
    end
    int_gpu = makeGPUarray(backend,int_total_float32)
    


    output_gpu = makeGPUarray(backend, zeros(Float32, nPoints, nPoints, nTotalSmallα))
    
    #


    #endregion
    
    #region launch GPU computation

    # scalars in Int32
    P      = Int32(nPoints)
    L      = Int32(nLs)
    nDim   = Int32(nCoordinates)
    nα     = Int32(nTotalSmallα)
    int1max = Int32(max_size[1])
    int2max = Int32(max_size[2])

    #@show typeof(output_gpu)         # must be MtlArray{Float32,3}
    #@show typeof(C_gpu)              # must be MtlArray{Float32,3}
    #@show typeof(int_gpu)            # must be MtlArray{Float32,5}
    #@show typeof(tableForLoop_gpu)   # must be MtlArray{Int32,3}
    #@show typeof(tableForPoints_gpu) # must be MtlMatrix{Int32}
    #@show typeof(P) typeof(L) typeof(nDim) typeof(nα) typeof(int1max) typeof(int2max)

    kernel! = windowContraction!(backend,(8,8,8))#,128,size(output_gpu))
    kernel!(output_gpu, C_gpu, int_gpu, tableForLoop_gpu,tableForPoints_gpu,
       P,L,nDim,nα,int1max,int2max,ndrange=size(output_gpu))
    KernelAbstractions.synchronize(backend)

    # Here output_gpu[x',x,eachα] ∑μ' ∑μ ∑l' ∑ l C[x',l',μ'] C[x,l,μ] ∏_iCoord K[iCoord][μ',μ,l'-n',l-n]
    # with (n',n) depends on eachα and x=η+μ

    @show "GPU computation of Ajiννᶜ: done"
    output = makeGPUarray(CPU(),output_gpu)
    #endregion

    #region contruct Ulocal

    # the order is: (νᶜ,) ν, i, j  here

    Ulocal = Array{Num,2}(undef,length(pointsIndices),NtypeofFields)
    for iField in eachindex(fields)
        newstring=split(string(fields[iField]),"(")[1]
        Ulocal[:,iField]=Symbolics.variables(Symbol(newstring),1:length(pointsIndices))
    end

    #endregion

    #region make Ajiννᶜ and AjiννᶜU symbolically (which we will soon remove!)
    Ajiννᶜ = Array{Num,3}(undef,length(pointsIndices),NtypeofFields,NtypeofExpr)
    AjiννᶜU = Array{Num,1}(undef,NtypeofExpr)
    
    # this is the cost function for ν point so the number of elements is just the number of expressions (governing equations)
    Ajiννᶜ .= 0
    AjiννᶜU .= 0
    indexLinearα = 1
   
    for iExpr ∈ 1:NtypeofExpr,iField ∈ 1:NtypeofFields
        α = bigα[iExpr,iField]
        for eachα ∈ α
            @show nodeValue=eachα.node
            for x ∈ pointsIndices
                xLinear = LI_points[vec2car(x)]
                
                localmapηᶜ=Dict()
                for iVar ∈ eachindex(vars)
                    localmapηᶜ[vars[iVar]]=varM[iVar,xLinear][]
                end
                for xᶜ ∈ pointsIndices
                    xᶜLinear = LI_points[vec2car(xᶜ)]
                    U_HERE = Ulocal[xᶜLinear,iField]                    
                    substitutedValue = substitute(nodeValue, localmapηᶜ)
                    Ajiννᶜ[xᶜLinear,iField,iExpr] += Float64(output[xᶜLinear,xLinear,indexLinearα])*substitutedValue
                    AjiννᶜU[iExpr]+= Ajiννᶜ[xᶜLinear,iField,iExpr] * U_HERE
                end

            end
            indexLinearα += 1
        end

    end


    #endregion

    coefs=(Ajiννᶜ=output,Ulocal=Ulocal,AjiννᶜU=AjiννᶜU)
    return coefs
end 



@kernel function windowContraction!(
    output::AbstractArray{Float32,3},   # Array{Float32,3} on GPU: (P,P,nalpha)
    C::AbstractArray{Float32,3},        # Array{Float32,3} on GPU: (P,L,P)
    int::AbstractArray{Float32,5},    # Array{Float32,5} on GPU: (nDim, int1max,int1max,int2max,int2max)
    table::AbstractArray{Int32,3}, # Array{Int32,3} on GPU: (2+nDim*2,L*L, nalpha)
    tablePoints::AbstractArray{Int32,2}, # Array{Int32,2} on GPU: (nDim,P)
    #sizeOf1DIntegrals::AbstractArray{Int32,1}, # max size of 1D WYYKK integrals
    P::Int32, L::Int32, nDim::Int32, nalpha::Int32, int1max::Int32,int2max::Int32
)

    (xᶜ, x, α) = @index(Global, NTuple)   # ← 3-D indices
    
    @inbounds begin
    
        if xᶜ ≤ size(output,1) && x ≤ size(output,2) && α ≤ size(output,3)

            acc = zero(eltype(output))

            for idx in 1:size(table, 2)  # iterate over all entries in table for this alpha
                i = table[1, idx, α]
                j = table[2, idx, α]
                

                if i >0 && j >0
                        
                    for μᶜ in 1:P
                        for μ in 1:P
                            prod_int = 1.0f0
                            for iDim in 1:nDim
                                
                                k = table[2+iDim, idx, α]
                                l = table[2+iDim+nDim, idx, α]
                                if k >0 && l >0
                    
                                    mᶜ=tablePoints[iDim,μᶜ]
                                    m=tablePoints[iDim,μ]

                                    prod_int *= int[iDim, mᶜ, m, k, l]
                        
                                end
                            end
                            acc += C[xᶜ, i, μᶜ] * C[x, j, μ] * prod_int  # <- corrected x' in first C
                        end
                    end
                
                end
            end

            output[xᶜ, x, α] = acc
        end
    end
end
