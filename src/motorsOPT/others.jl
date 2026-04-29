
timeDimensionString="t" 
# if the user does not want to use "t" for the time marching scheme, it should be changed
# and this "t" should be the last element in coordinates

function findCartesianDependency(expression,Ndimension,∂)
    expressionDependency=ones(Int,Ndimension)
    for iDimension in 1:Ndimension
        if typeof(expand_derivatives(∂[iDimension](expression))==0) <: Bool 
            if expand_derivatives(∂[iDimension](expression))==0
                expressionDependency[iDimension] = 0
            end
        end
    end
    return expressionDependency
end

function makeMixPartials(orders,coordinates;field=identity)
    # this function will give a matrix of mixed partial derivative operators 
    # coordinates should be an array of Symbolics variables 
    # orders is a matrix 


    Ndimension=length(coordinates)
    if length(orders)!== Ndimension
        @error "the highest orders array has not the same dimension as that of the coordinates"
    end

    ∂ = []
    for iDim in 1:Ndimension
        ∂ = push!(∂,Differential(coordinates[iDim]))
    end

    ∇ = Array{Any,Ndimension}(undef, Tuple(orders))
    R=CartesianIndices(∇)
        
    ∇ .= field
    for I in R
        for iDim in 1:Ndimension
            ∇[I] = (∂[iDim]^(I[iDim]-1))(∇[I])
        end
    end
   
    return ∇
end

function varMmaker(maxPointsUsed,coordinates,vars,∂) 
    # this will make an array of material coeffs for with a local Cartesian grid (max points used for a node)
    Ndimension = length(coordinates)

    R = CartesianIndices(Tuple(maxPointsUsed))

    varM=Array{Any,2}(undef,length(vars),length(R))
   
    for iVar in eachindex(vars)


        newstring=split(string(vars[iVar]),"(")[1]
     
        
        CartesianDependency=findCartesianDependency(vars[iVar],Ndimension,∂)
       
        smallVarM=Symbolics.variables(Symbol(newstring),1:length(R))
        for j in R
            linearJ=LinearIndices(R)[j]
            realJ=(car2svec(j).-1).*CartesianDependency .+1 # if there is no dependence on a direction, it should get the same name
            linearRealJ=LinearIndices(R)[CartesianIndex(realJ...)]
            smallVarM[linearJ]=smallVarM[linearRealJ]
        end
        varM[iVar,:]=smallVarM
    end
    
    return varM
end

function PDECoefFinder(orders,coordinates,expr,field,vars)
    # PDECoefFinder cannot detect the material partials × material partials for the moment!! 
    # I know how to do it, but eq. 40 should be then more generalised (kind of the product of partials of different materials)

    # maxPolynomialOrderMaterial is also a chelou thing, that I need to work on more systematically
    # like the powers of partials should also be included but here search for Rm[1], yeah, that's what I am doing


    Ndimension = length(coordinates)
    alpha=[]
    
    maxPolynomialOrderMaterial = 2*(maximum(orders)-1)
    ∇=makeMixPartials(orders,coordinates;field=field)
    R=CartesianIndices(∇)
    expr=mySimplify(expr)


    for i in R
        term_searched = ∇[i]

        tmpCoeff = myCoeff(expr,term_searched)
        if tmpCoeff !== 0
            isTmpCoeffAConstant=true
            for iVar in eachindex(vars)
                
                ∇m=makeMixPartials(orders,coordinates;field=vars[iVar]) # material partials
                Rm=CartesianIndices(∇m)
                for j in Rm
                    term_material_searched = ∇m[j]
                    tmpCoeffMaterial = myCoeff(tmpCoeff,term_material_searched)
                    
                    if tmpCoeffMaterial !==0
                        isTmpCoeffAConstant=false
                        isOKtoinclude =true
                        # This is to avoid partials of other material 
                        for jVar in eachindex(vars)
                            if jVar !== iVar
                                ∇n = makeMixPartials(orders,coordinates;field=vars[jVar]) 
                                for jj in Rm
                                    if jj !== Rm[1]
                                        term_material_searched_plus = ∇n[jj]
                                        differentialCoeff = myCoeff(tmpCoeffMaterial,term_material_searched_plus)
                                        if differentialCoeff !==0
                                            isOKtoinclude = false
                                        end
                                    end
                                end
                            end
                        end
                        if isOKtoinclude
                            #tmpCoeffMaterial=substitute(tmpCoeffMaterial,mapping)
                            specificMaterialTerm=tmpCoeffMaterial*vars[iVar]
                            tmpAlphaIJ = (node=specificMaterialTerm,nᶜ=j,n=i) # the famous n prime and n in the equation 56 or 40
                            alpha = push!(alpha,tmpAlphaIJ)
                        end
                    end
                end
                for matPower in 2:maxPolynomialOrderMaterial
                    tmpCoeffMaterial = myCoeff(tmpCoeff,vars[iVar]^matPower)
                    #tmpCoeffMaterial=substitute(tmpCoeffMaterial,mapping)
                    if tmpCoeffMaterial !==0
                        specificMaterialTerm=tmpCoeffMaterial*vars[iVar]^matPower
                        tmpAlphaIJ = (node=specificMaterialTerm,nᶜ=Rm[1],n=i) # the famous n prime and n in the equation 56 or 40
                        alpha = push!(alpha,tmpAlphaIJ)
                    end
                end
            end
            if isTmpCoeffAConstant
                specificMaterialTerm = tmpCoeff
                tmpAlphaIJ=(node=specificMaterialTerm,nᶜ=R[1],n=i)
                alpha = push!(alpha,tmpAlphaIJ)
            end
        end
    end
    alpha=unique(alpha)

    return alpha # varM: iVar and linearised cartesian indices
end 

function TaylorCoefInversion(coefInversionDict::Dict)

    # this function will give Float64 array (not Any)

    # the user might want to have a look on illposedTaylorCoefficientsInversion_legend, which is deprecated as of 10/06/2025.

    # as of 26/04/2026, it can compute for an arbitrary μPoints

    # based on the equation 27 (of the version 10/06/2025 FD2025 : \psi_{;\mu,\nu}^{(l)}[\nu+\mu]=\sum_\eta C_{\mu+\eta;\mu,\nu}^{(l)} \psi[\nu+\mu+\eta]), we need to perform this inversion anyways for all the point \mu inside L(\nu) (the concerned points for \nu)


    # be careful that pointsIndices is now a 1D array of integer vectors!!

    @unpack multiOrdersIndices, pointsIndices, μpointsIndices, Δ  = coefInversionDict


    @show typeof(pointsIndices)
    
    @show multiOrdersIndices, pointsIndices, μpointsIndices 
    
    numberOfEtas = length(pointsIndices)
    numberOfLs   = length(multiOrdersIndices)
    numberOfMus = length(μpointsIndices) 

    CˡηGlobal = Array{Float64,3}(undef,numberOfEtas,numberOfLs,numberOfMus)

    # this is the C^{(l)}_{\mu+\eta; μ, \nu}

    for μ_oneD in axes(CˡηGlobal,3)
        #@show typeof(multiOrdersIndices),typeof(pointsIndices)
        CˡηGlobal[:,:,μ_oneD]=TaylorCoefInversion(numberOfLs,numberOfEtas,multiOrdersIndices,pointsIndices,Δ,μpointsIndices[μ_oneD])
    end 

    return @strdict(CˡηGlobal)

end

function TaylorCoefInversion(numberOfLs,numberOfEtas,multiOrdersIndices,pointsIndices,Δ,μPoint)


    # the old version is : illposedTaylorCoefficientsInversionSingleCentre

    # in fact, available points depend on the position of μ (=k here), we need to 'mute' some points
    # with Y_μ
    
    # for this pointsIndices are filtered for every μ

    tmpTaylorExpansionCoeffs = Array{Any,2}(undef,numberOfLs,numberOfEtas)

    for iAvailablePoint in eachindex(pointsIndices)
        η_μ = pointsIndices[iAvailablePoint]
        η = Float64.(η_μ) .- μPoint
        distances= η .* Δ
        for j in multiOrdersIndices
            linearJ = LinearIndices(multiOrdersIndices)[j]
            orders = car2svec(j).-1
            numerator = prod(distances .^orders)
            denominator=prod(factorial.(orders))
            tmpTaylorCoeffs = numerator/denominator
            tmpTaylorExpansionCoeffs[linearJ,iAvailablePoint]=tmpTaylorCoeffs 
        end
    end

    # here we do the famous ill-posed inversion (ttttttt) 
    
    aa=transpose(tmpTaylorExpansionCoeffs)*tmpTaylorExpansionCoeffs
    aa=Num2Float64.(aa)
    typeof(aa),aa,size(aa)
    invaa= myInv(aa)
    tmpCˡηlocal=invaa*transpose(tmpTaylorExpansionCoeffs)

    return tmpCˡηlocal
end

function numbersOfTheExpression(equationCharacteristics,
                               trialFunctionsCharacteristics,
                               TaylorOptions1, TaylorOptions2)

    numbersOfTheSystem1=numbersOfTheExpression(equationCharacteristics,trialFunctionsCharacteristics,TaylorOptions1)
    numbersOfTheSystem2=numbersOfTheExpression(equationCharacteristics,trialFunctionsCharacteristics,TaylorOptions2)
    @unpack timeMarching,NtypeofExpr,NtypeofMaterialVariables,NtypeofFields,nCoordinates,Ndimension,pointsUsed,pointsμUsed,offsetsμUsed = numbersOfTheSystem1

    return (
        timeMarching = timeMarching,
        NtypeofExpr = NtypeofExpr,
        NtypeofMaterialVariables = NtypeofMaterialVariables,
        NtypeofFields = NtypeofFields,
        nCoordinates = nCoordinates,
        Ndimension = Ndimension,   # 🔥 key change
        pointsUsed = pointsUsed,
        pointsμUsed = pointsμUsed,
        offsetsμUsed = offsetsμUsed,
        pointsμᶜUsed = numbersOfTheSystem2.pointsμUsed,
        offsetsμᶜUsed = numbersOfTheSystem2.offsetsμUsed,
    )
end

function numbersOfTheExpression(equationCharacteristics,
                               trialFunctionsCharacteristics,
                               TaylorOptions)

    @unpack exprs, fields, vars, coordinates = equationCharacteristics
    @unpack pointsInSpace, pointsInTime = trialFunctionsCharacteristics
    @unpack pointsμInSpace, pointsμInTime,
            offsetμInΔyInSpace, offsetμInΔyInTime = TaylorOptions

    timeMarching = any(a -> a === timeDimensionString, string.(coordinates))

    NtypeofExpr = length(exprs)
    NtypeofMaterialVariables = length(vars)
    NtypeofFields = length(fields)

    Ndimension = length(coordinates)

    # 🔥 cleaner + no broadcast
    pointsUsed   = fill(pointsInSpace, Ndimension)
    pointsμUsed  = fill(pointsμInSpace, Ndimension)
    offsetsμUsed = fill(Float64(offsetμInΔyInSpace), Ndimension)

    if timeMarching
        pointsUsed[end]   = pointsInTime
        pointsμUsed[end]  = pointsμInTime
        offsetsμUsed[end] = Float64(offsetμInΔyInTime)
    end


    return (
        timeMarching = timeMarching,
        NtypeofExpr = NtypeofExpr,
        NtypeofMaterialVariables = NtypeofMaterialVariables,
        NtypeofFields = NtypeofFields,
        nCoordinates = Ndimension,
        Ndimension = Val(Ndimension),   # 🔥 key change
        pointsUsed = pointsUsed,
        pointsμUsed = pointsμUsed,
        offsetsμUsed = offsetsμUsed
    )
end

function investigateDependencies(equationCharacteristics,
                                 numbersOfTheSystem,
                                 trialFunctionsCharacteristics,
                                 TaylorOptionsμ,TaylorOptionsμᶜ)
    dependencies,ordersForSplinesμ,configsTaylorμ=investigateDependencies(equationCharacteristics,numbersOfTheSystem,trialFunctionsCharacteristics,TaylorOptionsμ)
    _,ordersForSplinesμᶜ,configsTaylorμᶜ= investigateDependencies(equationCharacteristics,numbersOfTheSystem,trialFunctionsCharacteristics,TaylorOptionsμᶜ)
    return dependencies,ordersForSplinesμ,configsTaylorμ,ordersForSplinesμᶜ,configsTaylorμᶜ
end


function investigateDependencies(equationCharacteristics,
                                 numbersOfTheSystem,
                                 trialFunctionsCharacteristics,
                                 TaylorOptions)

        @unpack Ndimension = numbersOfTheSystem

    return _investigateDependencies(Ndimension,
                                   equationCharacteristics,
                                   numbersOfTheSystem,
                                   trialFunctionsCharacteristics,
                                   TaylorOptions)
end

function _investigateDependencies(::Val{N},
                                 equationCharacteristics,
                                 numbersOfTheSystem,
                                 trialFunctionsCharacteristics,
                                 TaylorOptions) where N

    @unpack exprs,fields,vars,∂ = equationCharacteristics
    @unpack timeMarching,NtypeofExpr,NtypeofMaterialVariables,NtypeofFields,
            pointsUsed,pointsμUsed,offsetsμUsed = numbersOfTheSystem
    @unpack orderBtime,orderBspace = trialFunctionsCharacteristics
    @unpack YorderBtime,YorderBspace,supplementaryOrder = TaylorOptions

    # ---------------- Dependencies ----------------

    variableDependency = ones(Int, N)
    fieldDependency    = ones(Int, N)

    eachVariableDependency = zeros(Int, N, NtypeofMaterialVariables)
    eachFieldDependency    = zeros(Int, N, NtypeofFields)

    for iFields in 1:NtypeofFields
        dep = findCartesianDependency(fields[iFields], N, ∂)
        eachFieldDependency[:,iFields] = dep
        fieldDependency .*= (1 .- dep)
    end

    for iVars in 1:NtypeofMaterialVariables
        dep = findCartesianDependency(vars[iVars], N, ∂)
        eachVariableDependency[:,iVars] = dep
        variableDependency .*= (1 .- dep)
    end

    fieldDependency .= 1 .- fieldDependency
    variableDependency .= (1 .- variableDependency) .* fieldDependency

    # ---------------- Orders ----------------

    orderBspline  = zeros(Int, N)
    YorderBspline = zeros(Int, N)

    if timeMarching
        orderBspline[end] = orderBtime * fieldDependency[end]
        orderBspline[1:end-1] .= orderBspace .* fieldDependency[1:end-1]

        YorderBspline[end] = YorderBtime * fieldDependency[end]
        YorderBspline[1:end-1] .= YorderBspace .* fieldDependency[1:end-1]
    else
        orderBspline  .= orderBspace .* fieldDependency
        YorderBspline .= YorderBspace .* fieldDependency
    end

    pointsUsedForFields = (pointsUsed .- 1) .* fieldDependency .+ 1
    orderExpressions = pointsUsedForFields
    orderU = (orderExpressions .- 1) .+ (supplementaryOrder .* fieldDependency) .+ 1

    # ---------------- Taylor grids ----------------

    orderTaylors = Array{Any,N}(undef, Tuple(orderU))
    pointsInSpaceTime = Array{Any,N}(undef, Tuple(pointsUsedForFields))

    multiOrdersIndices = CartesianIndices(orderTaylors)
    multiPointsIndices = CartesianIndices(pointsInSpaceTime)

    availablePointsConfigurations = Vector{Any}()
    availableμPoints = Vector{Any}()
    availableμaxes = Vector{Any}()
    centrePointConfigurations = Int[]


    # in order to propose other middle points, we need to make a loop below

    numberGeometries = 1

    for i in 1:numberGeometries

        # ---------------- Middle point ----------------
    

        tmpVec = ((car2vec(multiPointsIndices[end]) .- 1) .÷ 2) .+ 1

        if timeMarching
            tmpVec[end] = car2vec(multiPointsIndices[end])[end] - 1
        end

        middleν = vec2car(tmpVec)

        push!(availablePointsConfigurations, car2svec.(multiPointsIndices))
        push!(centrePointConfigurations,
            LinearIndices(multiPointsIndices)[middleν])

        # ---------------- μ coordinates ----------------

        tmpμCoordinates = Array{SVector{N,Float64}}(undef, pointsμUsed...)

        tmpDistances = Float64.(availablePointsConfigurations[1][end] .- 1)
        tmpΔμ = tmpDistances .- 2 .* offsetsμUsed

        for i in 1:N
            if pointsμUsed[i] > 1
                tmpΔμ[i] /= (pointsμUsed[i] - 1)
            end
        end

        μaxes = ntuple(d -> begin
            [1.0 + offsetsμUsed[d] + (j - 1) * tmpΔμ[d] for j in 1:size(tmpμCoordinates, d)]
        end, N)

        for I in CartesianIndices(tmpμCoordinates)
            tmpμCoordinates[I] = SVector{N}(ntuple(d -> μaxes[d][I[d]], N))
        end


        push!(availableμPoints, tmpμCoordinates)
        push!(availableμaxes,μaxes)
    end
    
    # here we use the same bases for fields and materials
    availableμᶜPoints=availableμPoints
    availableμᶜaxes=availableμaxes

    # ---------------- Outputs ----------------

    dependencies = (
        variableDependency = variableDependency,
        fieldDependency = fieldDependency,
        eachVariableDependency = eachVariableDependency,
        eachFieldDependency = eachFieldDependency
    )

    ordersForSplines = (
        orderBspline = orderBspline,
        YorderBspline = YorderBspline,
        orderExpressions = orderExpressions,
        orderU = orderU
    )

    configsTaylor = (
        numberGeometries = numberGeometries,
        multiOrdersIndices = multiOrdersIndices, # we still need this since available points can differ
        availablePointsConfigurations = availablePointsConfigurations,
        centrePointConfigurations = centrePointConfigurations,
        availableμPoints = availableμPoints,
        availableμaxes = availableμaxes,
    )

    return dependencies, ordersForSplines, configsTaylor
end

function bigαFinder(equationCharacteristics,numbersOfTheSystem,ordersForSplines)

    @unpack NtypeofExpr, NtypeofFields = numbersOfTheSystem
    @unpack exprs,fields,vars,coordinates,∂ = equationCharacteristics
    @unpack orderExpressions = ordersForSplines

    bigα=Array{Any,2}(missing,NtypeofFields,NtypeofExpr)
    varM=nothing
    pointsUsedForFields=orderExpressions

    for iExpr in eachindex(exprs)
        for iField in eachindex(fields)
            
            tmpNonZeroAlphas=PDECoefFinder(orderExpressions,coordinates,exprs[iExpr],fields[iField],vars) 
            # we assume that the pointsUsedForFields represent the highest order of partials
            bigα[iField,iExpr]=unique(tmpNonZeroAlphas)
        end
    end
    varM=varMmaker(pointsUsedForFields,coordinates,vars,∂)
    return bigα,varM
end
#
#
#
#
#
# below are the obsolete functions


function numbersOfTheExpression_(equationCharacteristics,trialFunctionsCharacteristics,TaylorOptions)

    @unpack exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂² = equationCharacteristics
    @unpack orderBtime,orderBspace,pointsInSpace,pointsInTime = trialFunctionsCharacteristics
    @unpack YorderBtime,YorderBspace,supplementaryOrder,pointsμInSpace,pointsμInTime,offsetμInΔyInSpace,offsetμInΔyInTime = TaylorOptions

    timeMarching = any(a -> a === timeDimensionString, string.(coordinates))
   
    NtypeofExpr=length(exprs)   # number of governing equations
    NtypeofMaterialVariables=length(vars) # number of material coefficients
    NtypeofFields=length(fields) # number of unknown fields
    
    Ndimension = length(coordinates) # we do not change this for the moment, especially for the time-marching scheme
    pointsUsed = ones(Int, Ndimension).*(pointsInSpace)
    pointsμUsed = ones(Int, Ndimension).*(pointsμInSpace)
    offsetsμUsed = ones(Float64, Ndimension).*offsetμInΔyInSpace 
    if timeMarching
        pointsUsed[end]=pointsInTime
        pointsμUsed[end]=pointsμInTime
        offsetsμUsed[end]=offsetμInΔyInTime
    end


    if length(Δnum) !== Ndimension 
        @error "the numerical delta increment has not the same dimension!"
    end
    
    numbersOfTheSystem=(timeMarching=timeMarching,NtypeofExpr=NtypeofExpr,NtypeofMaterialVariables=NtypeofMaterialVariables,NtypeofFields=NtypeofFields,Ndimension=Ndimension,pointsUsed=pointsUsed,pointsμUsed=pointsμUsed,offsetsμUsed=offsetsμUsed)

    return numbersOfTheSystem

end

function investigateDependencies_(equationCharacteristics,numbersOfTheSystem,trialFunctionsCharacteristics,TaylorOptions)

    @unpack exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂² = equationCharacteristics
    @unpack timeMarching,NtypeofExpr,NtypeofMaterialVariables,NtypeofFields,Ndimension,pointsUsed,pointsμUsed,offsetsμUsed=numbersOfTheSystem
    @unpack orderBtime,orderBspace,pointsInSpace,pointsInTime = trialFunctionsCharacteristics
    @unpack YorderBtime,YorderBspace,supplementaryOrder,pointsμInSpace,pointsμInTime,offsetμInΔyInSpace,offsetμInΔyInTime = TaylorOptions


    #region overall coordinate dependency check

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
    YorderBspline=zeros(Int,Ndimension)

    if timeMarching
        orderBspline[Ndimension]=orderBtime*fieldDependency[Ndimension]
        orderBspline[1:Ndimension-1]=orderBspace*fieldDependency[1:Ndimension-1]
        YorderBspline[Ndimension]=YorderBtime*fieldDependency[Ndimension]
        YorderBspline[1:Ndimension-1]=YorderBspace*fieldDependency[1:Ndimension-1]
    else
        orderBspline[1:Ndimension]=orderBspace*fieldDependency[1:Ndimension]
        YorderBspline[1:Ndimension]=YorderBspace*fieldDependency[1:Ndimension]
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

    #region Preparation for Taylor expansion
    
        orderTaylors=Array{Any,Ndimension}(undef,Tuple(orderU))
        pointsInSpaceTime=Array{Any,Ndimension}(undef,Tuple(pointsUsedForFields))
        
        multiOrdersIndices=CartesianIndices(orderTaylors)

        availablePointsConfigurations = Array{Array{Vector{Int64},Ndimension},1}()
        availableμPoints = Array{Array{Vector{Float64},Ndimension},1}()
        centrePointConfigurations=Array{Int64,1}()

    #endregion


    #region Cartesian indices that can be available to use (normally: iGeometry=1)
    # here for the moment we just use one mid point in the space and time so iGeometry = 1 only

        multiPointsIndices=CartesianIndices(pointsInSpaceTime)
        # this is the whole local Cartesian grids (without any lacking points)
        
        tmpVecForMiddlePoint = ((car2svec(multiPointsIndices[end]).-1 ).÷2 ).+1 # only valid for testOnlyCentre
        midTimeCoord = nothing
        if timeMarching
            midTimeCoord=car2svec(multiPointsIndices[end])[end]-1
            tmpVecForMiddlePoint[end]=midTimeCoord
            #AjiννᶜU = Array{Num,2}(undef,length(multiPointsIndices)÷(midTimeCoord+1),NtypeofExpr)
        end
        #@show tmpVecForMiddlePoint 
        middleν=svec2car(tmpVecForMiddlePoint)

        availablePointsConfigurations=push!(availablePointsConfigurations,car2svec.(multiPointsIndices))
        centrePointConfigurations=push!(centrePointConfigurations,LinearIndices(multiPointsIndices)[middleν])

        # μ points for interpolated Taylor expansion
        # pointsμUsed = the number of μ ; offsetsμUsed = offsets from the extremeties (in Δy)

        tmpμCoordinates = Array{SVector{Ndimension,Float64}}(undef, pointsμUsed...)

        tmpDistancesInΔy = Float64.(availablePointsConfigurations[1][end].-ones(Float64,Ndimension))
        tmpΔμ=(tmpDistancesInΔy.-2.0*offsetsμUsed)
        for iCoord in 1:Ndimension
            if pointsμUsed[iCoord]>1
                tmpΔμ[iCoord]=tmpΔμ[iCoord]/(pointsμUsed[iCoord]-1)
            end
        end
        @show tmpμCoordinates
        for iμ in CartesianIndices(tmpμCoordinates)
            #tmpμCoordinates[iμ]=ones(Float64,Ndimension).+offsetsμUsed+(Float64.(car2svec(iμ)).-ones(Float64,Ndimension)).*tmpΔμ
            @show iμ
            #tmpμCoordinates[iμ]=offsetsμUsed
        end
        availableμPoints=push!(availableμPoints,tmpμCoordinates)
        
        #@show size(availablePointsConfigurations)
    #endregion





    

    dependencies = (variableDependency=variableDependency,fieldDependency=fieldDependency,eachVariableDependency=eachVariableDependency,eachFieldDependency=eachFieldDependency)
    ordersForSplines = (orderBspline=orderBspline,YorderBspline=YorderBspline,orderExpressions=orderExpressions,orderU=orderU)
    configsTaylor = (multiOrdersIndices=multiOrdersIndices,availablePointsConfigurations=availablePointsConfigurations,centrePointConfigurations=centrePointConfigurations,availableμPoints=availableμPoints)

    return dependencies,ordersForSplines,configsTaylor
end


function illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices;testOnlyCentre=true,Δ=nothing,timeMarching=false)

    # this function is deprecated, no more testOnlyCentre/timeMarching options are allowed (nor Symbolics Δ)

    # here we propose a big Taylor expansion matrix with Δcoordinates, symbolically (when Δ=nothing or Δ as a symbolic array)
    #   or numerically otherwise


    Ndimension=length(coordinates)

    if Δ === nothing
       Δ = Symbolics.variables(:Δ,1:Ndimension)
    else
        if length(Δ) !== Ndimension
            @error "the numerical delta increment has not the same dimension!"
        end
    end

    numberOfEtas = length(multiPointsIndices)
    numberOfLs   = length(multiOrdersIndices)

    #TaylorExpansionCoeffs=Array{Any,Ndimension*2}(undef, multiOrdersIndices[end],multiPointsIndices[end])
    #TaylorExpansionCoeffs=Array{Any,Ndimension*2}(undef, Tuple(vcat(collect(Tuple(multiOrdersIndices[end])),collect(Tuple(multiPointsIndices[end])))))
    
    CˡηGlobal = Array{Any,3}(undef,numberOfEtas,numberOfLs,numberOfEtas)
   

    tmpVecForMiddlePoint = (car2svec(multiPointsIndices[end]).-1 ).÷2 .+1 # only valid for testOnlyCentre
    midTimeCoord=nothing
    if timeMarching
        midTimeCoord=car2svec(multiPointsIndices[end])[end]-1
        tmpVecForMiddlePoint[end]=midTimeCoord
    end
    midK=svec2car(tmpVecForMiddlePoint)

    for k in multiPointsIndices
        linearK = LinearIndices(multiPointsIndices)[k]
        
        if !testOnlyCentre || k === midK || (timeMarching && car2svec(k)[end] === midTimeCoord && !testOnlyCentre) # because we cannot predict more than one futures
            CˡηGlobal[:,:,linearK]=illposedTaylorCoefficientsInversionSingleCentre(numberOfLs,numberOfEtas,multiOrdersIndices,multiPointsIndices,Δ,k)
        end
    end 

    if testOnlyCentre
        CˡηCentre = CˡηGlobal[:,:,LinearIndices(multiPointsIndices)[midK]]
        CˡηGlobal = nothing
        return CˡηCentre,Δ,multiOrdersIndices
    else
        #@show CˡηGlobal
        return CˡηGlobal,Δ,multiOrdersIndices
    end
end

function illposedTaylorCoefficientsInversionSingleCentre(numberOfLs,numberOfEtas,multiOrdersIndices,multiPointsIndices,Δ,k)

    # this function is deprecated as of 10/06/2025

    # this should be completely numerical
    TaylorExpansionCoeffs = Array{Num,2}(undef,numberOfLs,numberOfEtas)
    for i in multiPointsIndices
        linearI = LinearIndices(multiPointsIndices)[i]
        η = car2svec(i-k)
        distances= η .* Δ
        for j in multiOrdersIndices
            linearJ = LinearIndices(multiOrdersIndices)[j]
            orders = car2svec(j).-1
            numerator = prod(distances .^orders)
            denominator=prod(factorial.(orders))
            tmpTaylorCoeffs = numerator/denominator
            TaylorExpansionCoeffs[linearJ,linearI]=tmpTaylorCoeffs 

        end
    end
    # here we do the famous inversion (ttttttt) even though this code is essentially a forward problem
    
    aa=transpose(TaylorExpansionCoeffs)*TaylorExpansionCoeffs
    invaa= myInv(aa)
    Cˡηlocal=invaa*transpose(TaylorExpansionCoeffs)
    return Cˡηlocal
end



function spaceCoordinatesConversionfunctions(absorbingBoundaries, NdimensionMinusTime)
    offset_model = svec2car(absorbingBoundaries[1, 1:NdimensionMinusTime])
    #offset_empty = svec2car(spacePointsUsed)

    model2whole(a::CartesianIndex) = a + offset_model
    whole2model(a::CartesianIndex) = a - offset_model
    #whole2empty(a::CartesianIndex) = a + offset_empty
    #empty2whole(a::CartesianIndex) = a - offset_empty
    #model2empty(a::CartesianIndex) = whole2empty(model2whole(a))
    #empty2model(a::CartesianIndex) = whole2model(empty2whole(a))
    return(; model2whole, whole2model)
    #return (; model2whole, whole2model, whole2empty, empty2whole, model2empty, empty2model)
end

function BouncingCoordinates(a::CartesianIndex,PointsUsed)
    #
    # this will bounce the boundary inside the PointsUsed vector
    #
    # i.e. get the nearby coordinates inside the domain to fake the continuity

    if length(a) !== length(PointsUsed)
        @error "cannot bound this CartesianIndex due to the dimension mismatch"
    end
    avector=car2svec(a)
    for iCoord in eachindex(avector)
        if avector[iCoord] < 1
            avector[iCoord] = 1
        elseif avector[iCoord] > PointsUsed[iCoord]
            avector[iCoord] = PointsUsed[iCoord]
        end
    end
    a=svec2car(avector)
    return a
end