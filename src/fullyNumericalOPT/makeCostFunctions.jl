struct CouplingTable{T}
    rows::Vector{Int32}
    cols::Vector{Int32}
    vals::Vector{T}
end

struct MatrixFreeNumericalOperator{T}
    table::CouplingTable{T}
    size::Tuple{Int,Int}
    side::String
    backend::Symbol
    representation::Symbol
    geometry
end

struct SparseNumericalOperator{T,S}
    matrix::S
    table::CouplingTable{T}
    side::String
    backend::Symbol
    representation::Symbol
    geometry
end

function applyNumericalOperator!(y, op::MatrixFreeNumericalOperator, x)
    fill!(y, zero(eltype(y)))
    rows = op.table.rows
    cols = op.table.cols
    vals = op.table.vals
    @inbounds for k in eachindex(vals)
        y[rows[k]] += vals[k] * x[cols[k]]
    end
    return y
end

function applyNumericalOperator!(y, op::SparseNumericalOperator, x)
    mul!(y, op.matrix, x)
    return y
end

_as_symbol(x::Symbol) = x
_as_symbol(x::AbstractString) = Symbol(x)

_field_symbol(fieldSyms::AbstractArray, i) = fieldSyms[i]
_field_symbol(fieldSyms, i) = fieldSyms

function _numerical_value(x)
    v = Symbolics.value(x)
    return v isa Number ? v : error("non-numerical coefficient remains after material substitution: $v")
end

function _push_coupling!(rows, cols, vals, row, col, val)
    v = _numerical_value(val)
    iszero(v) && return nothing
    push!(rows, Int32(row))
    push!(cols, Int32(col))
    push!(vals, v)
    return nothing
end

function _finalize_numerical_operator(table::CouplingTable{T}, nrow::Int, ncol::Int, side, backend, representation, geometry) where T
    backend_sym = _as_symbol(backend)
    repr_sym = representation === :auto ? :sparse : _as_symbol(representation)

    if backend_sym == :mpi
        repr_sym = repr_sym == :auto ? :matrixfree : repr_sym
        return MatrixFreeNumericalOperator(table, (nrow, ncol), side, backend_sym, repr_sym, geometry)
    elseif repr_sym == :matrixfree
        return MatrixFreeNumericalOperator(table, (nrow, ncol), side, backend_sym, repr_sym, geometry)
    elseif repr_sym == :sparse || repr_sym == :blocksparse
        A = sparse(Int.(table.rows), Int.(table.cols), table.vals, nrow, ncol)
        if backend_sym == :cuda
            if isdefined(Main, :CUDA)
                A = Main.CUDA.CUSPARSE.CuSparseMatrixCSR(A)
            else
                @warn "backend=:cuda requested, but CUDA is not loaded in Main; returning CPU sparse matrix"
                backend_sym = :cpu
            end
        elseif backend_sym != :cpu
            throw(ArgumentError("backend must be :cpu, :cuda, or :mpi; got $backend"))
        end
        return SparseNumericalOperator(A, table, side, backend_sym, repr_sym, geometry)
    else
        throw(ArgumentError("representation must be :auto, :matrixfree, :sparse, or :blocksparse; got $representation"))
    end
end



function numericalOperatorConstruction(params::Dict)
    @unpack optRec,modelFam,absorbingBoundaries,maskedRegionInSpace=params
    backend = _paramget(params, :backend, :cpu)
    representation = _paramget(params, :representation, :auto)
    if !hasproperty(modelFam, :modelName) || modelFam.modelName === nothing
        modelFam = merge(modelFam, (; modelName = "model_" * Dates.format(now(), "yyyymmdd_HHMMSS")))
    end
    costLHS=numericalOperatorConstruction(optRec,modelFam,"left";absorbingBoundaries=absorbingBoundaries,maskedRegionInSpace=maskedRegionInSpace,backend=backend,representation=representation)
    costRHS=numericalOperatorConstruction(optRec,modelFam,"right";absorbingBoundaries=absorbingBoundaries,maskedRegionInSpace=maskedRegionInSpace,backend=backend,representation=representation)
    costfunctions=costLHS.costFunctions-costRHS.costFunctions
    fieldLHS=costLHS.場
    fieldRHS=costRHS.場
    champsLimité=costRHS.champsLimité
    numericalOperators=(left=costLHS.operator,right=costRHS.operator)
    numOperators=(costfunctions=costfunctions,fieldLHS=fieldLHS,fieldRHS=fieldRHS,champsLimité=champsLimité,numericalOperators=numericalOperators)
    return @strdict(numOperators)
end

function numericalOperatorConstruction(optRec,modelFam,side;absorbingBoundaries=nothing,maskedRegionInSpace=nothing,
    backend=:cpu,        # :cpu, :cuda, :mpi
    representation=:auto # :matrixfree, :sparse, :blocksparse
)

    numericalGeometry = prepareNumericalOperatorGeometry(optRec,modelFam,side; absorbingBoundaries, maskedRegionInSpace)

    @unpack νWhole, νGeometry, νRelative, localPointsIndices,
        ModelPoints, Models, maskingField, conv,
        timePointsUsedForOneStep, activeTimePoints, wholeMin, wholeMax,
        modelMin, wholeRegionPointsSpace = numericalGeometry

    @unpack lhs, rhs, numbersOfTheSystem, fieldNames = optRec["recette"]
    @unpack fields, extfields = fieldNames
    @unpack NtypeofExpr, NtypeofFields, NtypeofMaterialVariables = numbersOfTheSystem.numbersOfTheSystemL

    if side == "left"
        symA = lhs.Ajiννᶜ
        varM = lhs.varM
        fieldSyms = fields
    elseif side == "right"
        symA = rhs.Γjiννᶜ
        varM = rhs.varF
        fieldSyms = extfields
    else
        throw(ArgumentError("side must be \"left\" or \"right\", got $side"))
    end

    fieldLinearIndices = LinearIndices(
        (NtypeofFields, activeTimePoints, Tuple(wholeRegionPointsSpace)...)
    )
    residualLinearIndices = LinearIndices((NtypeofExpr, length(νWhole)))

    rows = Int32[]
    cols = Int32[]
    vals = Float64[]

    for iTestFunctions in eachindex(νWhole)
        νtmpWhole = νWhole[iTestFunctions]
        iGeometry = νGeometry[iTestFunctions]
        localPointsHere = localPointsIndices[iGeometry]
        middlepointHere = νRelative[iTestFunctions]
        localPointsSpaceHere = carDropDim.(localPointsHere)
        localPointsSpaceIndicesHere = CartesianIndices(Tuple(localPointsSpaceHere[end]))
        middlepointSpaceHere = svec2car(carDropDim(middlepointHere))
        iTimeMax = timePointsUsedForOneStep[iGeometry]

        νᶜtmpWhole = localPointsSpaceIndicesHere .+ (νtmpWhole .- middlepointSpaceHere)
        νᶜtmpModel = conv.whole2model.(νᶜtmpWhole)
        localLinearIndices = LinearIndices(Tuple(localPointsHere[end]))

        materialMapping = Dict()
        for iT in 1:iTimeMax
            for iVar in 1:NtypeofMaterialVariables
                spaceModelBouncedPoints = ModelPoints[1:end-1, iVar]
                iiT = ModelPoints[end, iVar] > 1 ? iT : 1
                νᶜtmpModelTruncated = BouncingCoordinates.(νᶜtmpModel, Ref(spaceModelBouncedPoints))

                for jPoint in νᶜtmpWhole
                    jPointLocal = jPoint - νtmpWhole + middlepointSpaceHere
                    jPointTLocal = carAddDim(jPointLocal, iT)
                    linearjPointTLocal = localLinearIndices[jPointTLocal]
                    materialMapping[varM[iVar, linearjPointTLocal]] =
                        Models[iVar][carAddDim(νᶜtmpModelTruncated[jPointLocal], iiT)]
                end
            end
        end

        for iExpr in 1:NtypeofExpr
            row = residualLinearIndices[iExpr, iTestFunctions]

            for iT in 1:iTimeMax
                for jPoint in νᶜtmpWhole
                    jPointLocal = jPoint - νtmpWhole + middlepointSpaceHere
                    jPointTLocal = carAddDim(jPointLocal, iT)
                    linearjPointTLocal = localLinearIndices[jPointTLocal]

                    jPointInWhole = is_all_less_than_or_equal(wholeMin, jPoint) &&
                                    is_all_less_than_or_equal(jPoint, wholeMax)
                    boundaryWeight = 0.0

                    if jPointInWhole
                        jPointModel = conv.whole2model(jPoint)
                        jPointInModel =
                            is_all_less_than_or_equal(modelMin, jPointModel) &&
                            is_all_less_than_or_equal(jPointModel, vec2car(ModelPoints[1:end-1, 1]))

                        if jPointInModel || iT == iTimeMax
                            boundaryWeight = 1.0
                        else
                            distance2 = distance2_point_to_box(
                                jPointModel,
                                modelMin,
                                vec2car(ModelPoints[1:end-1, 1]),
                            )
                            boundaryWeight = CerjanBoundaryCondition(distance2)
                        end
                    end

                    weight = jPointInWhole ? maskingField[jPoint] * boundaryWeight : 0.0
                    iszero(weight) && continue

                    for iField in 1:NtypeofFields
                        col = fieldLinearIndices[iField, iT, Tuple(jPoint)...]
                        coef = weight * substitute(symA[linearjPointTLocal, iField, iExpr, iGeometry], materialMapping)
                        _push_coupling!(rows, cols, vals, row, col, coef)
                    end
                end
            end
        end
    end

    table = CouplingTable(rows, cols, vals)
    nrow = NtypeofExpr * length(νWhole)
    ncol = NtypeofFields * activeTimePoints * length(νWhole)
    operator = _finalize_numerical_operator(table, nrow, ncol, side, backend, representation, numericalGeometry)

    fieldArrays = Array{Any,2}(undef, NtypeofFields, activeTimePoints)
    for iT in 1:activeTimePoints
        for iField in 1:NtypeofFields
            fieldArrays[iField, iT] = reshape(
                [_field_symbol(fieldSyms, iField) for _ in 1:length(νWhole)],
                Tuple(wholeRegionPointsSpace),
            )
        end
    end

    costFunctions = operator isa SparseNumericalOperator ? operator.matrix : operator
    return (
        costFunctions = costFunctions,
        operator = operator,
        couplingTable = table,
        場 = fieldArrays,
        champsLimité = nothing,
        geometry = numericalGeometry,
    )
end


function numericalOperatorConstruction_too_heavy(optRec,modelFam,side;absorbingBoundaries=nothing,maskedRegionInSpace=nothing)

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

    #region

    numericalGeometry = prepareNumericalOperatorGeometry(optRec,modelFam,side; absorbingBoundaries, maskedRegionInSpace)

    @unpack νWhole, νGeometry, νRelative, localPointsIndices, middlepoints,
        ModelPoints, Models, maskingField, conv,
        timePointsUsedForOneStep, wholeMin, wholeMax, modelMin,
        wholeRegionPoints, wholeRegionPointsSpace = numericalGeometry

    @unpack models, modelName, modelPoints = modelFam

    @unpack lhs, rhs, nodes, centresIndices, numbersOfTheSystem,fieldNames=optRec["recette"]

    @unpack fields, extfields = fieldNames
    
    @unpack timeMarching,nCoordinates,NtypeofExpr,NtypeofFields,NtypeofMaterialVariables=numbersOfTheSystem.numbersOfTheSystemL
    
    nConfigurations=numbersOfTheSystem.nConfigurations

    if side == "left"
        symA = lhs.Ajiννᶜ
        varM = lhs.varM
        Ulocal = lhs.Ulocal
        fieldSyms = fields
    elseif side == "right"
        symA = rhs.Γjiννᶜ
        varM = rhs.varF
        Ulocal = rhs.Flocal
        fieldSyms = extfields
    else
        throw(ArgumentError("side must be \"left\" or \"right\", got $side"))
    end




    #region
    # one operator per test point in space, for now
    NtestfunctionsInSpace = length(νWhole)
    costFunctions = Array{Any,2}(undef, NtypeofExpr, NtestfunctionsInSpace)
    #costFunctionsCoefs = Array{Number,6}(undef, NtestfunctionsInSpace, NtypeofFields, NtypeofExpr, iTimeMax, NtestfunctionsInSpace)
    costFunctionsCoefs .= 0
    Threads.@threads for iTestFunctions in eachindex(νWhole)

        localCost = Array{Any,1}(undef,NtypeofExpr)
        localCost .= 0.0
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
                            localCost[iExpr] += 
                                場[iField, iT][jPoint] * maskingField[jPoint] *
                                substitute(symA[linearjPointTLocal,iField,iExpr,iGeometry],tmpMapping)

                            #costFunctionsCoefs[jPoint] += NtypeofFields,  NtypeofExpr,iTimeMax, NtestfunctionsInSpace)
                        else
                            if iT == iTimeMax
                                #tmpMapping[Ulocal[linearjPointTLocal, iField]] =
                                #    場[iField, iT][jPoint] * maskingField[jPoint]
                                localCost[iExpr] += 
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
                                localCost[iExpr] += 
                                    場[iField, iT][jPoint] * maskingField[jPoint] *
                                    CerjanBoundaryCondition(distance2) *
                                    substitute(symA[linearjPointTLocal,iField,iExpr,iGeometry],tmpMapping)
                            end
                        end
                    end
                end

            end

            # centre-only address: always use the geometry centre now
            #costFunctions[iExpr, iTestFunctions] +=
            #    substitute(Av[iExpr,iGeometry], tmpMapping)
        end
        costFunctions[:, iTestFunctions] .= localCost
    end



    #endregion

    return (costFunctions=costFunctions,場=場,champsLimité=champsLimité)

end

function prepareNumericalOperatorGeometry(optRec, modelFam,side; absorbingBoundaries, maskedRegionInSpace)


    @unpack models, modelName, modelPoints = modelFam


    @unpack lhs, rhs, nodes, centresIndices, numbersOfTheSystem,fieldNames=optRec["recette"]

    @unpack fields, extfields = fieldNames

    
    if side=="left"
        symA=lhs.Ajiννᶜ
        varM=lhs.varM
        Ulocal=lhs.Ulocal
        fields=fields
    elseif side=="right"
        fields = nothing
        symA=rhs.Γjiννᶜ
        varM=rhs.varF
        Ulocal=rhs.Flocal
        fields=extfields
    end

    

    @unpack timeMarching,nCoordinates,NtypeofExpr,NtypeofFields,NtypeofMaterialVariables=numbersOfTheSystem.numbersOfTheSystemL
    
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
            size(models[iVar])
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

   


    #since everything is super clumsy, here we make several useful functions to change one coordinate to another
    
    conv=spaceCoordinatesConversionfunctions(absorbingBoundaries[:,1:end-1], newD-1)
    #endregion 

    #region

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

    fieldLinearIndices = LinearIndices(
        (NtypeofFields, activeTimePoints, Tuple(wholeRegionPointsSpace)...)
    )

    residualLinearIndices = LinearIndices(
        (NtypeofExpr, length(νWhole))
    )



    return (modelPoints=modelPoints,
    wholeRegionPoints=wholeRegionPoints,
    absorbingBoundaries=absorbingBoundaries,
    νWhole=νWhole,
    νGeometry=νGeometry,
    νRelative=νRelative,
    localPointsIndices=localPointsIndices,
    middlepoints=middlepoints,
    ModelPoints=ModelPoints,
    geometryPreference=geometryPreference,
    maskingField=maskingField,
    conv=conv,
    Models = Models,
    timePointsUsedForOneStep = timePointsUsedForOneStep,
    activeTimePoints = activeTimePoints,
    wholeRegionPointsSpace = wholeRegionPointsSpace,
    wholeMin = wholeMin,
    wholeMax = wholeMax,
    modelMin = modelMin,
    )
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