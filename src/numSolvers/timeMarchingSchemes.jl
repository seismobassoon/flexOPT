
using JLD2,DrWatson,SparseArrays,SparseDiffTools,LinearAlgebra
function prepareNumericalOperators(numOperators; T=Float32)
    @unpack costfunctions, fieldLHS, fieldRHS, champsLimité = numOperators

    costvec = vec(costfunctions)

    # current layout: fieldLHS[iField, iT] is an Array{Num,Nspace}
    timePointsUsedForOneStep = size(fieldLHS, 2)
    NField = size(fieldLHS, 1)

    spaceShape = size(fieldLHS[1,1])
    Rcoord = CartesianIndices(spaceShape)
    linearR = LinearIndices(Rcoord)
    NpointsSpace = length(Rcoord)
    NknownTime = max(timePointsUsedForOneStep - 1, 0)

    # symbolic only for compilation
    symbUnknownField = Array{Num,2}(undef, NpointsSpace, NField)
    symbKnownField = Array{Num,3}(undef, NpointsSpace, NField, NknownTime)

    for iT in 1:NknownTime
        for iField in 1:NField
            A = fieldLHS[iField, iT]
            for j in Rcoord
                lj = linearR[j]
                symbKnownField[lj, iField, iT] = A[j]
            end
        end
    end

    for iField in 1:NField
        A = fieldLHS[iField, timePointsUsedForOneStep]
        for j in Rcoord
            lj = linearR[j]
            symbUnknownField[lj, iField] = A[j]
        end
    end

    if champsLimité === nothing
        NforcePoints = NpointsSpace
        symbKnownForce = Array{Num,3}(undef, NforcePoints, NField, timePointsUsedForOneStep)
        for iT in 1:timePointsUsedForOneStep
            for iField in 1:NField
                A = fieldRHS[iField, iT]
                for j in Rcoord
                    lj = linearR[j]
                    symbKnownForce[lj, iField, iT] = A[j]
                end
            end
        end
    else
        NforcePoints = length(champsLimité[1,1])
        symbKnownForce = Array{Num,3}(undef, NforcePoints, NField, timePointsUsedForOneStep)
        for iT in 1:timePointsUsedForOneStep
            for iField in 1:NField
                for j in 1:NforcePoints
                    symbKnownForce[j, iField, iT] = champsLimité[iField, iT][j]
                end
            end
        end
    end

    unknownInputs = vec(symbUnknownField)
    knownFieldInputs = NknownTime == 0 ? Num[] : vec(symbKnownField)
    knownForceInputs = vec(symbKnownForce)
    allInputs = vcat(unknownInputs, knownFieldInputs, knownForceInputs)

    # numeric arrays only from here on
    unknown_template = zeros(T, NpointsSpace, NField)
    known_lhs_template = zeros(T, NpointsSpace, NField, NknownTime)
    known_rhs_template = zeros(T, NforcePoints, NField, timePointsUsedForOneStep)
    residual_template = zeros(T, length(costvec))
    
    residual_expr = build_function(costvec, allInputs; expression=Val(false))
    residual_tuple = eval(residual_expr)

    residual_fun = residual_tuple[1]
    residual_fun! = residual_tuple[2]

    return (
        residual_fun = residual_fun,
        residual_fun! = residual_fun!,
        unknown_template = unknown_template,
        known_lhs_template = known_lhs_template,
        known_rhs_template = known_rhs_template,
        residual_template = residual_template,
        spaceShape = spaceShape,
        NpointsSpace = NpointsSpace,
        NforcePoints = NforcePoints,
        NField = NField,
        timePointsUsedForOneStep = timePointsUsedForOneStep,
    )

end

function makePreparedResidual(prepared)
    Nunknown = length(prepared.unknown_template)
    NknownLHS = length(prepared.known_lhs_template)
    NknownRHS = length(prepared.known_rhs_template)

    function residual!(F, unknownField, knownField, knownForce)
        T = eltype(unknownField)
        allInputs = Vector{T}(undef, Nunknown + NknownLHS + NknownRHS)

        allInputs[1:Nunknown] .= vec(unknownField)
        allInputs[Nunknown+1:Nunknown+NknownLHS] .= T.(vec(knownField))
        allInputs[Nunknown+NknownLHS+1:end] .= T.(vec(knownForce))

        prepared.residual_fun!(F, allInputs)
        return nothing
    end

    return residual!
end

function prepareSparseJacobian(prepared)
    knownField = copy(prepared.known_lhs_template)
    knownForce = copy(prepared.known_rhs_template)

    nU = length(prepared.unknown_template)
    nF = length(prepared.residual_template)

    function wrapped!(F, U)
        T = eltype(U)
        unknownField = reshape(U, size(prepared.unknown_template))
        knownFieldT = T.(knownField)
        knownForceT = T.(knownForce)

        residual! = makePreparedResidual(prepared)
        residual!(F, unknownField, knownFieldT, knownForceT)
        return nothing
    end

    U0 = zeros(Float64, nU)
    F0 = zeros(Float64, nF)

    Usym = Vector{Num}(undef, nU)
    Fsym = Vector{Num}(undef, nF)
    sparsity = Symbolics.jacobian_sparsity(wrapped!, Fsym, Usym)
    S = sparse(map(!iszero, sparsity))

    J = spzeros(Float64, size(S)...)
    cache = ForwardColorJacCache(wrapped!, U0; sparsity=S)

    return J, cache
end


function timeStepOptimisationPrepared!(
    residual!,
    unknownField,
    knownField,
    knownForce,
    J,
    cache,
    F;
    nIteration=10,
    smallNumber=1e-8,
    boundaryConditionForced=false,
)
    nEq = length(F)
    normalisation = 1.0 / nEq
    r1 = 1.0

    U = vec(copy(unknownField))

    wrapped! = (F, U) -> begin
        unknownField .= reshape(U, size(unknownField))
        residual!(F, unknownField, knownField, knownForce)
        if boundaryConditionForced
            F[1] = U[1] - 1.0
            F[end] = U[end] - 1.0
        end
        nothing
    end

    for iter in 1:nIteration
        wrapped!(F, U)

        r = norm(F) * normalisation
        if iter == 1
            r1 = r
        end
        if r == 0.0 || r / r1 < smallNumber
            break
        end

        forwarddiff_color_jacobian!(J, wrapped!, U, cache)
        δU = -(J \ F)
        U .+= δU
    end

    unknownField .= reshape(U, size(unknownField))
    return nothing
end

function timeMarchingSchemePrepared(
    prepared,
    Nt,
    Δnum,
    modelName;
    videoMode=true,
    sourceType="Ricker",
    t₀=50,
    f₀=0.04,
    initialCondition=0.0,
    sourceFull=nothing,
    iExperiment=nothing,
    boundaryConditionForced=false,
)
    if !isdir(datadir("fieldResults"))
        mkdir(datadir("fieldResults"))
    end

    if isnothing(iExperiment)
        iExperiment = (iH=1, iCase=1, iPointsUsed=1)
    end

    @unpack iH, iCase, iPointsUsed = iExperiment
    experiment_name = "$(iH)_$(iCase)_$(iPointsUsed)"

    sequentialFileName = datadir("fieldResults", savename((Nt, Δnum..., modelName, sourceType, experiment_name), "jld2"))
    compactFileName = datadir("fieldResults", savename("compact", (Nt, Δnum..., modelName, sourceType, experiment_name), "jld2"))

    residual! = makePreparedResidual(prepared)
    J, cache = prepareSparseJacobian(prepared)

    NField = prepared.NField
    NpointsSpace = prepared.NpointsSpace
    NforcePoints = prepared.NforcePoints
    timePointsUsedForOneStep = prepared.timePointsUsedForOneStep
    spaceShape = prepared.spaceShape
    NknownTime = max(timePointsUsedForOneStep - 1, 0)

    unknownField = fill(Float64(initialCondition), size(prepared.unknown_template))
    knownField = fill(Float64(initialCondition), size(prepared.known_lhs_template))
    knownForce = fill(Float64(initialCondition), size(prepared.known_rhs_template))
    F = copy(prepared.residual_template)

    itVec = collect(1:Nt)
    t = (itVec .- 1) .* Δnum[end]

    sourceTime = nothing
    if sourceType == "Ricker"
        myRicker(x) = Ricker(x, t₀, f₀)
        sourceTime = myRicker.(t)
        if timePointsUsedForOneStep != 1
            prepend!(sourceTime, zeros(timePointsUsedForOneStep))
        end
    elseif sourceType == "Explicit"
        expected = (NforcePoints, NField, Nt)
        size(sourceFull) == expected || error("sourceFull should have size $expected")
    end

    if !isfile(sequentialFileName)
        fieldFile = jldopen(sequentialFileName, "w")

        hm = nothing
        if videoMode
            fig = Figure()
            ax = Axis(fig[1, 1])
            firstField = reshape(unknownField[:, 1], spaceShape...)
            hm = heatmap!(ax, Float32.(firstField), colormap=:deep, colorrange=(-1e-5, 1e-5))
            display(fig)
        end

        for it in itVec
            if sourceType == "Ricker"
                knownForce .= 0
                knownForce[1, 1, :] .= sourceTime[it:it+timePointsUsedForOneStep-1]
            else
                knownForce .= sourceFull[:, :, it:it+timePointsUsedForOneStep-1]
            end

            timeStepOptimisationPrepared!(residual!, unknownField, knownField, knownForce, J, cache, F; boundaryConditionForced=boundaryConditionForced)

            if NknownTime > 0
                if NknownTime > 1
                    knownField[:, :, 1:end-1] .= knownField[:, :, 2:end]
                end
                knownField[:, :, end] .= unknownField
            end

            newField = reshape(unknownField, NpointsSpace, NField)
            fieldFile["timestep_$it"] = reshape(newField, spaceShape..., NField)

            if videoMode
                firstField = reshape(unknownField[:, 1], spaceShape...)
                hm[1] = Float32.(firstField)
            end
        end

        close(fieldFile)
    end

    file = jldopen(sequentialFileName, "r")
    a = [file["timestep_$it"] for it in itVec]
    close(file)

    if videoMode
        fig = Figure()
        ax = Axis(fig[1, 1])
        hm = heatmap!(ax, Float32.(selectdim(a[1], ndims(a[1]), 1)), colormap=:plasma, colorrange=(-1e-6, 1e-6))
        Colorbar(fig[1, 2], hm)
        display(fig)
    end

    acompact = reduce(hcat, a)
    jldsave(compactFileName; acompact=acompact)

    return acompact
end

function timeStepOptimisationPrepared!(
    residual!,
    unknownField,
    knownField,
    knownForce,
    J,
    cache,
    F;
    nIteration=10,
    smallNumber=1e-8,
    boundaryConditionForced=false,
)
    nEq = length(F)
    normalisation = 1.0 / nEq
    r1 = 1.0

    U = vec(copy(unknownField))

    wrapped! = (F, U) -> begin
        T = eltype(U)

        unknownFieldT = reshape(U, size(unknownField))
        knownFieldT = T.(knownField)
        knownForceT = T.(knownForce)

        residual!(F, unknownFieldT, knownFieldT, knownForceT)

        if boundaryConditionForced
            F[1] = U[1] - one(T)
            F[end] = U[end] - one(T)
        end
        nothing
    end

    for iter in 1:nIteration
        wrapped!(F, U)

        r = norm(F) * normalisation
        if iter == 1
            r1 = r
        end
        if r == 0.0 || r / r1 < smallNumber
            break
        end

        forwarddiff_color_jacobian!(J, wrapped!, U, cache)
        δU = -(J \ F)
        U .+= δU
    end

    unknownField .= reshape(U, size(unknownField))
    return nothing
end
