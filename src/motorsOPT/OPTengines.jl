# this is the GPU version of OPTEngines/OPTnewEngines.jl 
# December 2025 Nobuaki Fuji (IPGP/UPC/IUF)



function OPTobj(operatorConfigurations::Dict)
    # this is just a wrapper for the OPTobj function below, for DrWatson package
    @unpack famousEquationType, Δnum, orderBtime, orderBspace, WorderBtime,WorderBspace,supplementaryOrder,pointsInSpace, pointsInTime,IneedExternalSources, iExperiment= operatorConfigurations

    exprs,fields,vars,extexprs,extfields,extvars,coordinates,∂,∂² = famousEquations(famousEquationType)

    TaylorOptions=(WorderBtime=WorderBtime,WorderBspace=WorderBspace,supplementaryOrder=supplementaryOrder)
    trialFunctionsCharacteristics=(orderBtime=orderBtime,orderBspace=orderBspace,pointsInSpace=pointsInSpace,pointsInTime=pointsInTime)
    @time operatorData=OPTobj(exprs,fields,vars; coordinates=coordinates,trialFunctionsCharacteristics=trialFunctionsCharacteristics,TaylorOptions=TaylorOptions,Δnum = Δnum,iExperiment=iExperiment)
    #AjiννᶜU=operatorData[1]
    #utilities=operatorData[2]   

    operatorForceData=nothing
    # if you do not want to apply external forces, it is possible to skip below
    if IneedExternalSources 
        @time operatorForceData=OPTobj(extexprs,extfields,extvars; coordinates=coordinates,trialFunctionsCharacteristics=trialFunctionsCharacteristics,TaylorOptions=TaylorOptions,Δnum = Δnum,iExperiment=iExperiment)  
        #@show Γg = operatorForceData[1]
        #utilitiesForce = operatorForceData[2]
    end
    eqInfo=(exprs=exprs,fields=fields,vars=vars,extexprs=extexprs,extfields=extfields,extvars=extvars,coordinates=coordinates)
    operators=(operatorPDE=operatorData, operatorForce=operatorForceData,eqInfo=eqInfo)
    return @strdict(operators)
end

function OPTobj(exprs,fields,vars; coordinates=(x,y,z,t), TaylorOptions=(WorderBtime=1,WorderBspace=1,supplementaryOrder=2), trialFunctionsCharacteristics=(orderBtime=1,orderBspace=1, pointsInSpace=2,pointsInTime=2),iExperiment =nothing)

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
        eachFieldDependency[:,iFields]=findCartesianDependency(fields[iFields],Ndimension)
        fieldDependency = fieldDependency .* (ones(Int,Ndimension).-eachFieldDependency[:,iFields])
    end


    for iVars in 1:NtypeofMaterialVariables
        eachVariableDependency[:,iVars]=findCartesianDependency(vars[iVars],Ndimension)
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
    varM=varMmaker(pointsUsedForFields,coordinates,vars)
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


    Δ = Δnum
    

    AjiννᶜU=[]
    Ulocal=[]

    for iConfigGeometry in eachindex(availablePointsConfigurations) 
        pointsIndices=availablePointsConfigurations[iConfigGeometry]
        middleLinearν=centrePointConfigurations[iConfigGeometry]
        #varM is given above for the max number of points used 
        tmpAjiννᶜU,tmpUlocal=AuSymbolic(coordinates,multiOrdersIndices,pointsIndices,multiPointsIndices,middleLinearν,Δ,varM,bigα,orderBspline,WorderBspline,NtypeofExpr,NtypeofFields)
        AjiννᶜU=push!(AjiννᶜU,tmpAjiννᶜU)
        Ulocal=push!(Ulocal,tmpUlocal)
    end


    #endregion



    #region outputs
    
    utilities=(middlepoint=middleν,middlepointLinear=centrePointConfigurations[1],localPointsIndices=multiPointsIndices,localMaterials=varM,localFields=Ulocal[1])
   
    return AjiννᶜU,utilities
    

    #endregion
    
end
