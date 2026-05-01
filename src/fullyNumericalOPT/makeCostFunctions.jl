function numericalOperatorConstruction(params::Dict)
    @unpack optRec,modelFam,absorbingBoundaries,maskedRegionInSpace=params
    costLHS=quasiNumericalOperatorConstruction(optRec,modelFam,"left";absorbingBoundaries=absorbingBoundaries,maskedRegionInSpace=maskedRegionInSpace)
    costRHS=quasiNumericalOperatorConstruction(optRec,modelFam,"right";absorbingBoundaries=absorbingBoundaries,maskedRegionInSpace=maskedRegionInSpace)
    costfunctions=costLHS.costFunctions-costRHS.costFunctions
    fieldLHS=costLHS.場
    fieldRHS=costRHS.場
    champsLimité=costRHS.場
    numOperators=(costfunctions=costfunctions,fieldLHS=fieldLHS,fieldRHS=fieldRHS,champsLimité=champsLimité)
    return @strdict(numOperators)
end



function numericalOperatorConstruction(optRec,modelFam,side;absorbingBoundaries=nothing,maskedRegionInSpace=nothing)

    #region general introduction
    #
    # 01/05/2026 Nobuaki Fuji
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

    @unpack models, modelName, modelPoints = modelFam


    @unpack lhs, rhs, nodes, centresIndices, numbersOfTheSystem,fieldNames=optRec["recette"]

    @unpack fields, extfields = fieldNames

    fields = nothing
    if side=="left"
        symA=lhs.Ajiννᶜ
        varM=lhs.varM
        Ulocal=lhs.Ulocal
        fields=fields
    elseif side=="right"
        symA=rhs.Γjiννᶜ
        varM=rhs.varF
        Ulocal=rhs.Flocal
        fields=extfields
    end

    

    @unpack timeMarching,nCoordinates,NtypeofExpr,NtypeofFields,NtypeofMaterialVariables=numbersOfTheSystem.numbersOfTheSystemL
    NtypeofExtFields=numbersOfTheSystem.numbersOfTheSystemL.NtypeofFields
    NtypeofExtMaterialVariables=numbersOfTheSystem.numbersOfTheSystemR.NtypeofMaterialVariables
    nConfigurations=numbersOfTheSystem.nConfigurations
    # normally the geometry configurations should be proposed in the preferred order
    nGeometry=nConfigurations
    Ndimension=nCoordinates

    

    # the last coordinate should be cosidered as time

    localPointsIndices = Vector{Any}(undef, nGeometry)
    middlepoints = Vector{Any}(undef, nGeometry)
    lhs_CartesianDependencies = lhs.CartesianDependencies
    rhs_CartesianDependencies = rhs.CartesianDependencies

    newD=Ndimension

    if !timeMarching
        newD = Ndimension + 1
        modelPoints = (modelPoints..., 1)

        lhs_CartesianDependencies = vcat(lhs.CartesianDependencies,zeros(Int, 1, NtypeofMaterialVariables),)
        rhs_CartesianDependencies = vcat(rhs.CartesianDependencies,zeros(Int, 1, NtypeofFields),)

        for iGeometry in 1:nGeometry
            localPointsIndices[iGeometry] = [SVector{newD,Int}(p..., 1) for p in nodes[iGeometry]]
            selected = nodes[iGeometry][centresIndices[iGeometry]]
            middlepoints[iGeometry] = SVector{newD,Int}(selected..., 1)
        end
    else
        
        for iGeometry in 1:nGeometry
            localPointsIndices[iGeometry] = nodes[iGeometry]
            selected = nodes[iGeometry][centresIndices[iGeometry]]
            middlepoints[iGeometry] = selected
        end
    end



    if length(models) !== NtypeofMaterialVariables 
        @error "Each material has to have its own model"
    end
    
    Models=Array{Any,1}(undef,NtypeofMaterialVariables)
    ModelPoints = Array{Int,2}(undef, newD, NtypeofMaterialVariables)
    
    for iVar ∈ 1:NtypeofMaterialVariables
        CartesianDependency=lhs_CartesianDependencies[:,iVar]
        if ndims(models[iVar]) !== sum(CartesianDependency)
            @error "model parameter dimension is not what you declared in the equation!"
        end
        if sum(CartesianDependency) == 0 # when it is a constant
            tmpModel = Array{Any,newD}(undef, (ones(Int, newD)...)...)
            ModelPoints[:, iVar] = ones(Int, newD)
            tmpModel[vec2car(ones(Int, newD))] = models[iVar]
            Models[iVar]=tmpModel
        else
            @show size(models[iVar])
            @show newCoords=expandVectors(size(models[iVar]),CartesianDependency)
            @show ModelPoints[:,iVar] = newCoords
            @show tmpModel=reshape(models[iVar],newCoords...)
            @show Models[iVar]=tmpModel

            for iCoord in eachindex(newCoords)
                if newCoords[iCoord]!== modelPoints[iCoord] && newCoords[iCoord] !== 1
                    @error "the model should have the same dimension! (or constant)"
                end
            end
        end
    end
    
    #endregion
    @show Models, ModelPoints


    
    #region construction of the fields

    wholeRegionPoints = nothing

    if absorbingBoundaries === nothing
        wholeRegionPoints=modelPoints
        absorbingBoundaries = zeros(Int,2, newD)
    elseif absorbingBoundaries === "CerjanBoundary"
        wholeRegionPoints=modelPoints
        absorbingBoundaries = ones(Int,2, newD-1)*CerjanGridPoints
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
        wholeRegionPoints=modelPoints.+ sum(absorbingBoundaries, dims=1)[:]
    end
    wholeRegionPointsSpace=wholeRegionPoints[1:end-1]

    #endregion


    #region 

    # Preferred geometry at each spatial point.
    # Later this can be filled from a table or a classifier.
    geometryPreference = fill(1, Tuple(wholeRegionPointsSpace))

    spacePointsUsed=Vector{Any}(undef, nGeometry)
    timePointsUsedForOneStep=Vector{Any}(undef, nGeometry)
    for iGeometry ∈ 1:nGeometry
        localPointVecs = localPointsIndices[iGeometry]
        spacePointsUsed[iGeometry] = localPointVecs[end][1:end-1]
        timePointsUsedForOneStep[iGeometry] = localPointVecs[end][end]
    end

    # For now, all points use geometry 1, so this is the active time depth.
    activeTimePoints = timePointsUsedForOneStep[1]

    # Fields
 
    場 = Array{Any,2}(undef, NtypeofFields, activeTimePoints)

    for it ∈ 1:activeTimePoints
        for iField ∈ 1:NtypeofFields
            newstring = split(string(fields[iField]), "(")[1] * "_mod_t=" * string(it)
            場[iField, it] = Symbolics.variables(
                Symbol(newstring),
                Base.OneTo.(Tuple(wholeRegionPointsSpace))...
            )
        end
    end


    #since everything is super clumsy, here we make several useful functions to change one coordinate to another
    
    conv=spaceCoordinatesConversionfunctions(absorbingBoundaries[:,1:end-1], newD-1)
    #endregion 

    # Useful point lists
    PointsSpace = CartesianIndices(Tuple(wholeRegionPointsSpace))
    NpointsSpace = length(PointsSpace)

    # For now, test functions are still identified with spatial points
    νWhole = collect(PointsSpace)
    wholeMin = νWhole[1]
    wholeMax = νWhole[end]
    modelMin = CartesianIndex(ones(Int, newD-1)...)



    # Pointwise preferred geometry and corresponding relative centre
    νGeometry = Vector{Int}(undef, NpointsSpace)
    νRelative = Vector{Any}(undef, NpointsSpace)

    for iPoint in eachindex(νWhole)
        νtmpWhole = νWhole[iPoint]
        iGeometry = geometryPreference[νtmpWhole]
        νGeometry[iPoint] = iGeometry

        # For now take the first centre of that geometry.
        # Later this can be refined if one point wants another centre within the same geometry.
        νRelative[iPoint] = middlepoints[iGeometry]
    end

 

    #endregion


    #region making a maskingField (for limited source areas, boundary conditions, etc.)

    maskingField=Array{Any,newD-1}(undef,Tuple(wholeRegionPointsSpace)) # maskingField is defined only for whole domain
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

    #region
    # one operator per test point in space, for now
    NtestfunctionsInSpace = length(νWhole)
    costFunctions = Array{Any,2}(undef, NtypeofExpr, NtestfunctionsInSpace)
    costFunctions .= 0
    for iTestFunctions in eachindex(νWhole)
        iPoint = iTestFunctions
        νtmpWhole = νWhole[iPoint]

        # preferred geometry for this test point
        iGeometry = νGeometry[iPoint]
        localPointsHere = localPointsIndices[iGeometry]
        middlepointHere = νRelative[iPoint]              # for now = middlepoints[iGeometry][1]
        localPointsSpaceHere = carDropDim.(localPointsHere)
        localPointsSpaceIndicesHere = CartesianIndices(Tuple(localPointsSpaceHere[end]))
        middlepointSpaceHere = svec2car(carDropDim(middlepointHere))
        
        # if time is appended, the active time depth depends on the geometry
        iTimeMax = timePointsUsedForOneStep[iGeometry]

        # shift local stencil to the current whole-space test point
     
        νᶜtmpWhole = localPointsSpaceIndicesHere .+ (νtmpWhole .- svec2car(carDropDim(middlepointHere)))
        νᶜtmpModel = conv.whole2model.(νᶜtmpWhole)

        # linear indexing on the chosen local stencil
        localLinearIndices = LinearIndices(Tuple(localPointsHere[end]))


        for iExpr in 1:NtypeofExpr
            tmpMapping = Dict()

            for iT in 1:iTimeMax
                # material variables
                for iVar in 1:NtypeofMaterialVariables
                    spaceModelBouncedPoints = ModelPoints[1:end-1, iVar]
                    iiT = ModelPoints[end, iVar] > 1 ? iT : 1

                    
                    νᶜtmpModelTruncated = BouncingCoordinates.(νᶜtmpModel, Ref(spaceModelBouncedPoints))

                    for jPoint in νᶜtmpWhole
                        jPointLocal = jPoint - νtmpWhole + svec2car(carDropDim(middlepointHere))
                        jPointTLocal = carAddDim(jPointLocal, iT)
                        linearjPointTLocal = localLinearIndices[jPointTLocal]

                        tmpMapping[varM[iVar, linearjPointTLocal]] =
                            Models[iVar][carAddDim(νᶜtmpModelTruncated[jPointLocal], iiT)]
                    end


                end

                # fields
                for jPoint in νᶜtmpWhole
                    jPointLocal = jPoint - νtmpWhole + middlepointSpaceHere
                    jPointTLocal = carAddDim(jPointLocal, iT)
                    linearjPointTLocal = localLinearIndices[jPointTLocal]

                    jPointInWhole = is_all_less_than_or_equal(wholeMin, jPoint) &&
                                    is_all_less_than_or_equal(jPoint, wholeMax)

                    if !jPointInWhole
                        for iField in 1:NtypeofFields
                            tmpMapping[Ulocal[linearjPointTLocal, iField]] = 0.0
                        end
                        continue
                    end

                    jPointModel = conv.whole2model(jPoint)

                    for iField in 1:NtypeofFields
                        if is_all_less_than_or_equal(modelMin, jPointModel) &&
                        is_all_less_than_or_equal(jPointModel, vec2car(ModelPoints[1:end-1, 1]))

                            #tmpMapping[Ulocal[linearjPointTLocal, iField]] =
                            #    場[iField, iT][jPoint] * maskingField[jPoint]
                            costFunctions[iExpr, iTestFunctions] += 
                                場[iField, iT][jPoint] * maskingField[jPoint] *
                                substitute(symA[linearjPointTLocal,iField,iExpr,iGeometry],tmpMapping)
                        else
                            if iT == iTimeMax
                                #tmpMapping[Ulocal[linearjPointTLocal, iField]] =
                                #    場[iField, iT][jPoint] * maskingField[jPoint]
                                costFunctions[iExpr, iTestFunctions] += 
                                    場[iField, iT][jPoint] * maskingField[jPoint] *
                                    substitute(symA[linearjPointTLocal,iField,iExpr,iGeometry],tmpMapping)
                            else
                                distance2 = distance2_point_to_box(
                                    jPointModel,
                                    modelMin,
                                    vec2car(ModelPoints[1:end-1, 1]),
                                )
                                #tmpMapping[Ulocal[linearjPointTLocal, iField]] =
                                #    場[iField, iT][jPoint] *
                                #    maskingField[jPoint] *
                                #    CerjanBoundaryCondition(distance2)
                                costFunctions[iExpr, iTestFunctions] += 
                                    場[iField, iT][jPoint] * maskingField[jPoint] *
                                    CerjanBoundaryCondition(distance2) *
                                    substitute(symA[linearjPointTLocal,iField,iExpr,Geometry],tmpMapping)
                            end
                        end
                    end
                end

            end

            # centre-only address: always use the geometry centre now
            #costFunctions[iExpr, iTestFunctions] +=
            #    substitute(Av[iExpr,iGeometry], tmpMapping)
        end
    end



    #endregion



    return (costFunctions=costFunctions,場=場,champsLimité=champsLimité)

end

carDropDim(v::SVector{N,T}) where {N,T} = SVector{N-1,T}(ntuple(i -> v[i], N-1))

function spaceCoordinatesConversionfunctions(absorbingBoundaries, NdimensionMinusTime)
    offset_model = vec2car(absorbingBoundaries[1, 1:NdimensionMinusTime])
    #offset_empty = vec2car(spacePointsUsed)

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
    avector=car2vec(a)
    for iCoord in eachindex(avector)
        if avector[iCoord] < 1
            avector[iCoord] = 1
        elseif avector[iCoord] > PointsUsed[iCoord]
            avector[iCoord] = PointsUsed[iCoord]
        end
    end
    a=vec2car(avector)
    return a
end