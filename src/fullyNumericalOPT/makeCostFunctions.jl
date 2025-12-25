
function quasiNumericalOperatorConstruction(operators,modelName,models,forceModels,famousEquationType,modelPoints,IneedExternalSources;maskedRegionForFieldInSpace = nothing,maskedRegionForSourcesInSpace=nothing,iExperiment=nothing)

    NpointsUsed=iExperiment.iPointsUsed
    # this is a big wrapper that reads the semi symbolic expressions to give a set of numerical operators (with symbolic expression in time)
    # which will call wrappers of onstructingNumericalDiscretisedEquations(Masked)

    operators=operators["operators"]
    operatorPDE,operatorForce,eqInfo = operators
    exprs,fields,vars,extexprs,extfields,extvars,coordinates=eqInfo

    AjiννᶜU,utilities=operatorPDE
    if IneedExternalSources 
        Γg,utilitiesForce=operatorForce
    end

    lhsConfigurations = @strdict semiSymbolicOpt=AjiννᶜU coordinates modelName models fields vars famousEquationType modelPoints utilities maskedRegion=maskedRegionForFieldInSpace NpointsUsed

    #numOperators,file = @produce_or_load(constructingNumericalDiscretisedEquations,lhsConfigurations,datadir("numOperators",savename(lhsConfigurations));filename = config -> savename(lhsConfigurations; ignores=["vars", "fields"]))
    numOperators = myProduceOrLoad(constructingNumericalDiscretisedEquations,lhsConfigurations,"numOperators","lhs")


    # left-hand side, which is far more recyclable than r.h.s.
    
    costfunctionsLHS,fieldLHS,_=numOperators["numOperators"]
    

    costfunctionsRHS = similar(costfunctionsLHS)
    costfunctionsRHS .= 0.
    fieldRHS = similar(fieldLHS)
    fieldRHS .= 0.
    
    champsLimité = nothing

    if IneedExternalSources 
        rhsConfigurations = @strdict semiSymbolicOpt=Γg coordinates modelName models=forceModels fields=extfields vars=extvars famousEquationType modelPoints utilities=utilitiesForce maskedRegion=maskedRegionForSourcesInSpace NpointsUsed
        numOperators,file=produce_or_load(constructingNumericalDiscretisedEquations,rhsConfigurations,datadir("numOperators",savename(rhsConfigurations));filename = config -> savename("source",rhsConfigurations; ignores=["vars", "fields"]))

        numOperators = myProduceOrLoad(constructingNumericalDiscretisedEquations,rhsConfigurations,"numOperators","rhs")
       
        costfunctionsRHS,fieldRHS,champsLimité=numOperators["numOperators"]

    end
  
    #@show size(costfunctionsLHS),size(costfunctionsRHS)
    #@show costfunctionsRHS[1,430],costfunctionsRHS[1,431],costfunctionsRHS[1,434]
    #costfunctions=0.#costfunctionsLHS[1,1]-costfunctionsLHS[1,1]
    costfunctions = costfunctionsLHS .- costfunctionsRHS
    #numOperators=(costfunctions=costfunctions)
    return costfunctions,fieldLHS,fieldRHS,champsLimité
end


function getModelPoints(models,pointsInTime,timeMarching)
    
    fakeNt = 1
    if timeMarching
        fakeNt = pointsInTime+1
        modelPoints = (size(models)...,fakeNt) # Nx, Ny etc thing. Nt is also mentioned and it should be the last element!
    else
        modelPoints = (size(models)...,1)
    end
    return modelPoints 
end

function constructingNumericalDiscretisedEquations(config::Dict)
    # just a wrapper
    @unpack semiSymbolicOpt,coordinates,modelName,models,fields,vars,famousEquationType,modelPoints,utilities, maskedRegion, NpointsUsed = config
   
    costfunctions,場,champsLimité=constructingNumericalDiscretisedEquations(semiSymbolicCoefs,myEquationInside,models;initialCondition=0.0)
    numOperators=(costfunctions=costfunctions,場=場,champsLimité=champsLimité)

    #@show costfunctions
    return @strdict(numOperators)
end

function constructingNumericalDiscretisedEquations(semiSymbolicCoefs,myEquationInside,models,modelPoints;absorbingBoundaries=nothing,initialCondition=0.0)

    # if modelPoints = nothing -> models[1] will be the reference for the modelPoints (in space)
    # if timeMarching modelPoints will get an additional dimension (>1)
    # if !timeMarching modelPoints will get an additional dimension (=1)

    #region todo list
    #todo list
    # 
    # this function is tooooooo complicated! I think I can simplify very much this!
    #
    #
    #  need to work on the bc, same like the masked thing (limited region of source)
    #
    # absorbing boundaries : I think we can already put the bc inside the numerical operators but be careful with the time marching: search for weightingCerjan
    # 
    # need extend to 4 points with the same test functions (3 points) -> staggered grid
    #  
    # I have to include some complex initial condition for 場
    #
    # have to write:
    #  function illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices,midPoint,Δ)
    #endregion
    
    #region general introduction
    #
    #
    # after the construction of local (semi-)symbolic expressions with linearised operators#
    # here we will read the model parameters and construct the numerical operators
    #
    # Nobuaki Fuji @ IPGP/UPC/IUF since 2024
    #
    # 
    # encouraged by Thibault Duretz @ U. Frankfurt Goethe, Kurama Okubo @ NIED
    #
    # Julia hackathon October 2024, March 2025
    #
    #

    # coordinates: Model: the real model domain; Whole: computation domain with absorbing boundaries; 
    #              Empty: Whole + some more points to avoid missing reference to the field and material (they should be just zeros)
    #

    # intermediate presentations: IPGP-CIA workshop October 2024; IPGP-ERI workshop November 2024; lighttalk @ systemI December 2024
    #             EGU @ Vienna May 2025
    #    
    #     Fuji & Duretz in preparation
    #
    #
    #
    #
    #endregion

    #region unpacking, N-dimensionalising all the models 

    # ce dont j'ai besoin
    semiSymbolicsOperators,maskedRegionInSpace
    
    # not yet resolved:


    testOnlyCentre = false
 
    if size(semiSymbolicsOperators)[1] === 1
        testOnlyCentre = true
    elseif ndims(semiSymbolicsOperators) !==2
        @error "the semi symbolic operators are not computed correctly!"
    end



    coordinates=myEquationInside.coordinates
    fields=myEquationInside.fields
    vars=myEquationInside.vars

    utilities=semiSymbolicCoefs["output"].utilities

    @unpack middlepoint,middlepointLinear,localPointsIndices,localMaterials,localFields,timeMarching = utilities
    timeMarching = utilities.timeMarching


    # not like allsame = all(s -> s == size.(models)[1], size.(models))

   
    # the last coordinate should be cosidered as time

    if !timeMarching
        localPointsIndices=CartesianIndices(Tuple([car2vec(localPointsIndices[end]);1]))
        middlepoint=CartesianIndex([car2vec(middlepoint);1]...)
        tmpT=Symbolics.variable(timeDimensionString)
        coordinates = (coordinates...,tmpT)
        #modelPoints = (modelPoints...,1)
        tmp_del = Symbolics.variable("∂"*timeDimensionString)
        tmp_del = Differential(tmpT)
        ∂ .= push!(∂,tmp_del)
        tmp_del_2 = Symbolics.variable("∂"*timeDimensionString*"²")
        tmp_del_2 = Differential(tmpT)*Differential(tmpT)
        ∂² .= push!(∂²,tmp_del_2)
    end

    #@show coordinates,∂,∂²

    localPointsSpaceIndices=CartesianIndices(Tuple(car2vec(localPointsIndices[end])[1:end-1]))
    Ndimension=length(coordinates)
    
   
    modelPoints=collect(modelPoints)

    Models=[]

    NtypeofMaterialVariables = length(vars)
    NtypeofFields = length(fields)
    NtypeofExpr = size(semiSymbolicsOperators)[end]

    Models=Array{Any,1}(undef,NtypeofMaterialVariables)
    ModelPoints=Array{Int,2}(undef,Ndimension,NtypeofMaterialVariables)

    if length(models) !== NtypeofMaterialVariables 
        @error "Each material has to have its own model"
    end
    
    
    for iVar in eachindex(vars)
        CartesianDependency=findCartesianDependency(vars[iVar],length(coordinates))
        if ndims(models[iVar]) === CartesianDependency
            @error "model parameter dimension is not what you declared in the equation!"
        end
        if sum(CartesianDependency) === 0 # when it is a constant
            tmpModel=Array{Any,Ndimension}(undef,(ones(Int, Ndimension)...)...)
            ModelPoints[:,iVar] = ones(Int, Ndimension)
            tmpModel[vec2car(ones(Int, Ndimension))] = models[iVar]
            Models[iVar]=tmpModel
        else
            #@show models[iVar],iVar,CartesianDependency, vars[iVar]
            newCoords=expandVectors(size(models[iVar]),CartesianDependency)
            ModelPoints[:,iVar] = newCoords
 
            tmpModel=reshape(models[iVar],newCoords...)
            Models[iVar]=tmpModel

            for iCoord in eachindex(newCoords)
                if newCoords[iCoord]!== modelPoints[iCoord] && newCoords[iCoord] !== 1
                    @error "the model should have the same dimension! (or constant)"
                end
            end
        end
     
    end

    #endregion

    #region construction of the fields

    wholeRegionPoints = nothing

    if absorbingBoundaries === nothing
        wholeRegionPoints=modelPoints
        absorbingBoundaries = zeros(Int,2, Ndimension)
    elseif absorbingBoundaries === "CerjanBoundary"
        wholeRegionPoints=modelPoints
        absorbingBoundaries = ones(Int,2, Ndimension-1)*CerjanGridPoints
        absorbingBoundaries=[absorbingBoundaries; 0 0]
        wholeRegionPoints=modelPoints.+sum(absorbingBoundaries,1) 
    else
        # absorbingBoundaries should be two column array 
        if size(absorbingBoundaries)[1] !== 2
            @error "you have to give us the left and right values for absorbing boundaries"
        elseif size(absorbingBoundaries)[2] !== size(modelPoints)[1] && !timeMarching
            @error "you have to give us the values for each direction for absorbing boundaries"
        elseif size(absorbingBoundaries)[2] === size(modelPoints)[1]-1 && timeMarching
            absorbingBoundaries=[absorbingBoundaries; 0 0]
        end
        wholeRegionPoints=modelPoints.+sum(absorbingBoundaries,1) 
    end


    # some useful stuff


    spacePointsUsed=car2vec(localPointsIndices[end])[1:end-1]
    timePointsUsedForOneStep=car2vec(localPointsIndices[end])[end]
    wholeRegionPointsSpace=wholeRegionPoints[1:end-1]
   
 
    # we need to put the left and right regions in order that centre ν configuration can pass

    #emptyRegionPointsSpace=wholeRegionPointsSpace.+ 2 .* spacePointsUsed

    #場dummy=Array{Any,2}(undef,NtypeofFields,timePointsUsedForOneStep)
    場 = Array{Any,2}(undef,NtypeofFields,timePointsUsedForOneStep)
    

    for it in 1:timePointsUsedForOneStep
        for iField in eachindex(fields)
            newstring=split(string(fields[iField]),"(")[1]*"_mod"*"_t="*string(it)
            場[iField,it]=Symbolics.variables(Symbol(newstring),Base.OneTo.(Tuple(wholeRegionPointsSpace))...)
        end
    end

    #since everything is super clumsy, here we make several useful functions to change one coordinate to another
    
    conv=spaceCoordinatesConversionfunctions(absorbingBoundaries[:,1:end-1], Ndimension-1)

    #endregion 

    #region making a maskingField (for limited source areas, boundary conditions, etc.)

    maskingField=Array{Any,Ndimension-1}(undef,Tuple(wholeRegionPointsSpace)) # maskingField is defined only for whole domain
    champsLimité = nothing
    if maskedRegionInSpace === nothing
        maskingField .= 1.0
    elseif typeof(maskedRegionInSpace) === Array{CartesianIndex,1}
        champsLimité = Array{Any,2}(undef,NtypeofFields,timePointsUsedForOneStep)
        for it in 1:timePointsUsedForOneStep
            for iField in eachindex(fields)
                newstring=split(string(fields[iField]),"(")[1]*"_mod_limited"*"_t="*string(it)
                champsLimité[iField,it] = Array{Any,1}(undef,length(maskedRegionInSpace))
            end
        end
        maskingField .= 0.0
        tmpIndex=1
        for iSpace in maskedRegionInSpace
            jSpace = conv.model2whole(iSpace)
            maskingField[jSpace] =1.0
            for it in 1:timePointsUsedForOneStep
                for iField in eachindex(fields)
                    
                    #tmpChampsLimitéContents= (jSpace,場[iField,it][jSpace])
                    champsLimité[iField,it][tmpIndex]=場[iField,it][jSpace]
                end
            end
            tmpIndex += 1
        end
    else
        @error "maskedRegionInSpace should be a 1D array of CartesianIndex (if it is CartesianIndices, you need to collect(Tuple()))"
    end


    #endregion

    #region relative ν to be considered, especially around the boundaries, useful for the following sections

    PointsSpace=CartesianIndices(Tuple(wholeRegionPointsSpace))
    NpointsSpace=length(PointsSpace) # number of points in space

    NtestfunctionsInSpace=NpointsSpace # this assumption is valid only for test functions related to grid points

    νWhole=Array{Any,1}(undef,NtestfunctionsInSpace) # the coordinate in wholeRegionPointsSpace: for the moment mapping from testfunction to ν is bijective

    # below is only for the bijective projection between test functions and ν

    for iPoint in eachindex(νWhole)

        νWhole[iPoint] = PointsSpace[iPoint] # this should be not true for higher B-spline test functions

    end


    νRelative=Array{Any,1}(undef,NtestfunctionsInSpace) # the relative coordinate to take (the coordinate used for the semi-symbolic operator derivation)
    νRelative.=middlepoint


    
    #endregion

    #region we construct the numerical operators for each test function that is related to its corresponding point

    # first we compute the νRelative more seriously if testOnlyCentre we might do nothing at all

    if testOnlyCentre
        # If we compute only the operators without boundaries, we use kind of 'truncated' crazy operators 
        # derived at the centre point and we do not talk about it, just believe the absorbing boundaries
        # like, tant pis, il n'y a pas de points donc j'ignore juste !

        # maybe we do not do anything!

    else

        # we use this clause only if we are interested in a serious boundary conditions, 
        # i.e. not the artificial cartesian box boundaries (we can just apply some stupid absorbing boundaries in that case)
        # Hence this clause should be more generalised even it kills the performance
        #  like, we give an array of free surface or discontinuities in CartesianIndex arrays 
        # and we look all the points concerned

        # we need to explore everywhere in the wholeRegionPoints! Free surface etc. should be very much affected 
        #
      
        boundaryPointsSpace=[]
        for iDimSpace in 1:Ndimension-1 # we take care of the boundaries of the Cartesian box (it should be the same for internal/external topography)
            # points concerned
            leftstart=1
            leftend=spacePointsUsed[iDimSpace]÷2
            rightstart=NpointsSpace[iDimSpace]-spacePointsUsed[iDimSpace]÷2+1
            rightend=NpointsSpace[iDimSpace]
            # suppose that the domain is sufficiently big (maybe it is not the case for some crazy topography ...)
            for iCoord in range(leftstart,leftend)
                


            end

            for iCoord in range(rightstart,rightend)

            end
        end
    end

    # here the number of test functions should not be necessarily the number of points but I will work later

    costFunctions=Array{Any,2}(undef,NtypeofExpr,NtestfunctionsInSpace)

    #@show semiSymbolicsOperators
    #@show localMaterials[1,15],localFields,size(localFields)
    #@show Models[1][10,15,1]

    for iTestFunctions in eachindex(νWhole)
        # here each test function is connected to one ν point 
        # We need to be careful that this can be no more true for different basis functions other than linear B-spline
        iPoint = iTestFunctions 
        νtmpWhole=νWhole[iPoint]
        # be careful with the two lines above, they are based on the assumption that each test function is linked to one collocated point


        #νtmpModel=conv.whole2model(νtmpWhole)
        νᶜtmpWhole = localPointsSpaceIndices .+ (νtmpWhole - carDropDim(νRelative[iPoint])) # this is the shift vector
        νᶜtmpModel = conv.whole2model.(νᶜtmpWhole)

        # examine νᶜtmpWhole if it is out of the range 



        for iExpr in eachindex(semiSymbolicsOperators[1,:])

            tmpMapping=Dict()


            for iT in 1:timePointsUsedForOneStep
                
                for iVar in eachindex(vars)
                    
                    spaceModelBouncedPoints=ModelPoints[1:end-1,iVar]

                    if ModelPoints[end,iVar] > 1
                        iiT=iT
                    else
                        iiT = 1
                    end


                    # model parameters should be bounced at the whole region limits
                    νᶜtmpModelTruncated = BouncingCoordinates.(νᶜtmpModel, Ref(spaceModelBouncedPoints))

                    for jPoint in νᶜtmpWhole
                        jPointLocal = jPoint - νtmpWhole + carDropDim(νRelative[iPoint])
                        jPointTLocal = carAddDim(jPointLocal,iT)
                        linearjPointTLocal=LinearIndices(localPointsIndices)[jPointTLocal]
              
                        tmpMapping[localMaterials[iVar,linearjPointTLocal]] = Models[iVar][carAddDim(νᶜtmpModelTruncated[jPointLocal],iiT)]
                        
                    end

                end

                for jPoint in νᶜtmpWhole
                    #@show iPoint, jPoint
                    jPointLocal = jPoint - νtmpWhole + carDropDim(νRelative[iPoint])
                    jPointTLocal = carAddDim(jPointLocal,iT)
                    linearjPointTLocal=LinearIndices(localPointsIndices)[jPointTLocal]
                    #jPointT=carAddDim(jPoint,iT)
                    #linearjPointT=LinearIndices(localPointsIndices)[jPointT]
                    for iField in eachindex(fields)
                        if is_all_less_than_or_equal(CartesianIndex(ones(Int, Ndimension-1)...),conv.whole2model(jPoint)) && is_all_less_than_or_equal(conv.whole2model(jPoint),vec2car(ModelPoints[1:end-1]))
                            # when it is inside the model domain box
                            tmpMapping[localFields[linearjPointTLocal,iField]] = 場[iField,iT][jPoint]*maskingField[jPoint]

                        elseif is_all_less_than_or_equal(νWhole[1],jPoint) && is_all_less_than_or_equal(jPoint,νWhole[end])
                            # if it is in the absorbing boundary zones we apply a simple Cerjan
                            if iT === timePointsUsedForOneStep # the last one (the future) will be using un-weighted operators
                                tmpMapping[localFields[linearjPointTLocal,iField]] = 場[iField,iT][jPoint]*maskingField[jPoint]
                            else
                                distance2 = distance2_point_to_box(conv.whole2model(jPoint),CartesianIndex(ones(Int, Ndimension-1)...), vec2car(ModelPoints[1:end-1]))
                                tmpMapping[localFields[linearjPointTLocal,iField]] = 場[iField,iT][jPoint]*maskingField[jPoint]*CerjanBoundaryCondition(distance2)
                            end
                        else
                            tmpMapping[localFields[linearjPointTLocal,iField]]=0.
                            #jPoint, νWhole[1],νWhole[end]
                        end
                    end
                end

                        

                        #場[iField,it]=string_as_varname(newstring, Array{Any,Ndimension-1}(undef,Tuple(wholeRegionPointsSpace)))
                        #
                        #νᶜtmpWholeMissing = ReplacerHorsLimiteParMissing(νᶜtmpWhole,PointsSpace[end])
                        # field values are defined only at the whole region and not at the Empty
                        #replace!(x -> x>0.2 ? missing : x, Array{Union{Float64, Missing}}(A) )
                        #replace!(x -> )
                        
                        
                
               
            end
            tmpAddress = nothing
            if testOnlyCentre
                tmpAddress = 1
            else
                tmpAddress=carDropDim(νRelative[iPoint])
            end
            costFunctions[iExpr,iTestFunctions]=substitute(semiSymbolicsOperators[tmpAddress,iExpr],tmpMapping)
            # be careful that semiSymbolicsOperators could be 2D
        end
    end

    #@show costFunctions

    #endregion

    return costFunctions,場,champsLimité

end

