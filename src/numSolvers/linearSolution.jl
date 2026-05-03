using Symbolics
using SparseArrays
using LinearAlgebra
using JLD2
using DrWatson

# Optional plotting support inside timeMarchingSchemeLinear
# Uncomment if you want videoMode=true.
# using CairoMakie
# CairoMakie.activate!()

function prepareLinearSystem(numOperators; T=Float64)
    @unpack costfunctions, fieldLHS, fieldRHS, champsLimité = numOperators

    costvec = vec(costfunctions)

    timePointsUsedForOneStep = size(fieldLHS, 2)
    NField = size(fieldLHS, 1)

    spaceShape = size(fieldLHS[1, 1])
    Rcoord = CartesianIndices(spaceShape)
    linearR = LinearIndices(Rcoord)
    NpointsSpace = length(Rcoord)
    NknownTime = max(timePointsUsedForOneStep - 1, 0)

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
        NforcePoints = length(champsLimité[1, 1])
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
    knownInputs = vcat(knownFieldInputs, knownForceInputs)

    # residual = A*u + c  =>  A*u = -c
    A_sym = Symbolics.jacobian(costvec, unknownInputs)
    zero_map = Dict(u => 0 for u in unknownInputs)
    c_sym = Symbolics.substitute.(costvec, Ref(zero_map))
    b_sym = .-c_sym

    A_expr = build_function(vec(A_sym), knownInputs; expression=Val(false))
    b_expr = build_function(b_sym, knownInputs; expression=Val(false))

    A_tuple = eval(A_expr)
    b_tuple = eval(b_expr)

    A_fun = A_tuple[1]
    A_fun! = A_tuple[2]

    b_fun = b_tuple[1]
    b_fun! = b_tuple[2]

    nUnknown = length(unknownInputs)
    nEq = length(costvec)

    known_lhs_template = zeros(T, NpointsSpace, NField, NknownTime)
    known_rhs_template = zeros(T, NforcePoints, NField, timePointsUsedForOneStep)

    A_template = zeros(T, nEq, nUnknown)
    b_template = zeros(T, nEq)

    return (
        A_fun = A_fun,
        A_fun! = A_fun!,
        b_fun = b_fun,
        b_fun! = b_fun!,
        A_template = A_template,
        b_template = b_template,
        known_lhs_template = known_lhs_template,
        known_rhs_template = known_rhs_template,
        spaceShape = spaceShape,
        NpointsSpace = NpointsSpace,
        NforcePoints = NforcePoints,
        NField = NField,
        timePointsUsedForOneStep = timePointsUsedForOneStep,
    )
end

function makeLinearEvaluator(prepared)
    nKnownField = length(prepared.known_lhs_template)
    nKnownForce = length(prepared.known_rhs_template)
    knownInputs = zeros(eltype(prepared.b_template), nKnownField + nKnownForce)

    function eval!(A, b, knownField, knownForce)
        knownInputs[1:nKnownField] .= vec(knownField)
        knownInputs[nKnownField+1:nKnownField+nKnownForce] .= vec(knownForce)

        prepared.A_fun!(vec(A), knownInputs)
        prepared.b_fun!(b, knownInputs)
        return A, b
    end

    return eval!
end

function evaluateLinearSystem(prepared, knownField, knownForce; sparse_output=false)
    eval! = makeLinearEvaluator(prepared)
    A = copy(prepared.A_template)
    b = copy(prepared.b_template)
    eval!(A, b, knownField, knownForce)
    if sparse_output
        A = sparse(A)
    end
    return A, b
end

function applyBoundaryConditionForced!(A, b; leftValue=1.0, rightValue=1.0)
    A[1, :] .= 0
    A[1, 1] = 1
    b[1] = leftValue

    A[end, :] .= 0
    A[end, end] = 1
    b[end] = rightValue

    return A, b
end

function prepareConstantMatrix(preparedLin;
    boundaryConditionForced=false,
    sparse_output=true,
    leftValue=1.0,
    rightValue=1.0,
)
    eval! = makeLinearEvaluator(preparedLin)

    knownField0 = zero(preparedLin.known_lhs_template)
    knownForce0 = zero(preparedLin.known_rhs_template)

    A = copy(preparedLin.A_template)
    btmp = copy(preparedLin.b_template)

    eval!(A, btmp, knownField0, knownForce0)

    if sparse_output
        A = sparse(A)
    end

    if boundaryConditionForced
        applyBoundaryConditionForced!(A, btmp; leftValue=leftValue, rightValue=rightValue)
    end

    factor = lu(A)
    return A, factor
end

function timeMarchingSchemeLinear(
    preparedLin,
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
    sparse_output=true,
    boundaryConditionForced=false,
    assume_constant_matrix=true,
)
    mkpath(datadir("fieldResults"))


    if isnothing(iExperiment)
        iExperiment = (iH=1, iCase=1, iPointsUsed=1)
    end

    @unpack iH, iCase, iPointsUsed = iExperiment
    experiment_name = "$(iH)_$(iCase)_$(iPointsUsed)"

    sequentialFileName = datadir("fieldResults", savename((Nt, Δnum..., modelName, sourceType, experiment_name, "linear"), "jld2"))
    compactFileName = datadir("fieldResults", savename("compact", (Nt, Δnum..., modelName, sourceType, experiment_name, "linear"), "jld2"))

    NField = preparedLin.NField
    NpointsSpace = preparedLin.NpointsSpace
    NforcePoints = preparedLin.NforcePoints
    timePointsUsedForOneStep = preparedLin.timePointsUsedForOneStep
    spaceShape = preparedLin.spaceShape
    NknownTime = max(timePointsUsedForOneStep - 1, 0)

    knownField = fill(Float64(initialCondition), size(preparedLin.known_lhs_template))
    knownForce = fill(Float64(initialCondition), size(preparedLin.known_rhs_template))
    unknownField = fill(Float64(initialCondition), NpointsSpace, NField)

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
        fig = nothing
        if videoMode
            fig = Figure()
            ax = Axis(fig[1, 1])
            firstField = reshape(unknownField[:, 1], spaceShape...)
            hm = heatmap!(ax, Float32.(firstField), colormap=:deep, colorrange=(-1e-5, 1e-5))
            display(fig)
        end

        eval! = makeLinearEvaluator(preparedLin)
        Awork = copy(preparedLin.A_template)
        b = copy(preparedLin.b_template)

        factor = nothing
        if assume_constant_matrix
            _, factor = prepareConstantMatrix(
                preparedLin;
                boundaryConditionForced=boundaryConditionForced,
                sparse_output=sparse_output,
            )
        end

        for it in itVec
            if sourceType == "Ricker"
                knownForce .= 0.0
                knownForce[1, 1, :] .= sourceTime[it:it+timePointsUsedForOneStep-1]
            else
                knownForce .= sourceFull[:, :, it:it+timePointsUsedForOneStep-1]
            end

            if assume_constant_matrix
                preparedLin.b_fun!(b, vcat(vec(knownField), vec(knownForce)))
                if boundaryConditionForced
                    b[1] = 1.0
                    b[end] = 1.0
                end
            else
                eval!(Awork, b, knownField, knownForce)
                A = sparse_output ? sparse(Awork) : copy(Awork)
                if boundaryConditionForced
                    applyBoundaryConditionForced!(A, b; leftValue=1.0, rightValue=1.0)
                end
                factor = lu(A)
            end

            u = factor \ b
            unknownField .= reshape(u, NpointsSpace, NField)

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

    acompact = reduce(hcat, a)
    jldsave(compactFileName; acompact=acompact)

    return acompact
end
