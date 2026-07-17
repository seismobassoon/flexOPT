
using LinearAlgebra
using SparseArrays
using Statistics
using CairoMakie
CairoMakie.activate!()

function point_source_weights(preparedLin, point::CartesianIndex; normalise=:none)
    weights = zeros(Float64, preparedLin.NforcePoints)
    LI = LinearIndices(preparedLin.spaceShape)
    weights[LI[point]] = 1.0
    if normalise == :sum
        s = sum(weights); !iszero(s) && (weights ./= s)
    elseif normalise == :l1
        s = sum(abs, weights); !iszero(s) && (weights ./= s)
    elseif normalise == :max
        s = maximum(abs, weights); !iszero(s) && (weights ./= s)
    end
    return weights
end

function gaussian_field(spaceShape, center; sigma=12.0, amplitude=1.0)
    G = zeros(Float64, spaceShape)
    c = Tuple(center)
    for I in CartesianIndices(spaceShape)
        d2 = sum((Tuple(I) .- c).^2)
        G[I] = amplitude * exp(-0.5 * d2 / sigma^2)
    end
    return G
end

function prepare_fd2d_acoustic_baseline(velocity, delta; source_scale=:dt2)
    spaceShape = size(velocity)
    length(spaceShape) == 2 || error("baseline helper currently expects a 2D velocity model")
    dx, dz, dt = Float64.(delta)
    nPoints = prod(spaceShape)
    LI = LinearIndices(spaceShape)

    rowsA = Int[]; colsA = Int[]; valsA = Float64[]
    rowsL = Int[]; colsL = Int[]; valsL = Float64[]
    rowsR = Int[]; colsR = Int[]; valsR = Float64[]
    col_known(point, it) = point + (it - 1) * nPoints

    for I in CartesianIndices(spaceShape)
        p = LI[I]
        v = Float64(velocity[I])
        cx = (v * dt / dx)^2
        cz = (v * dt / dz)^2

        push!(rowsA, p); push!(colsA, p); push!(valsA, 1.0)
        push!(rowsL, p); push!(colsL, col_known(p, 1)); push!(valsL, 1.0)
        push!(rowsL, p); push!(colsL, col_known(p, 2)); push!(valsL, -2.0 + 2cx + 2cz)

        i, j = Tuple(I)
        if i > 1
            q = LI[CartesianIndex(i - 1, j)]
            push!(rowsL, p); push!(colsL, col_known(q, 2)); push!(valsL, -cx)
        end
        if i < spaceShape[1]
            q = LI[CartesianIndex(i + 1, j)]
            push!(rowsL, p); push!(colsL, col_known(q, 2)); push!(valsL, -cx)
        end
        if j > 1
            q = LI[CartesianIndex(i, j - 1)]
            push!(rowsL, p); push!(colsL, col_known(q, 2)); push!(valsL, -cz)
        end
        if j < spaceShape[2]
            q = LI[CartesianIndex(i, j + 1)]
            push!(rowsL, p); push!(colsL, col_known(q, 2)); push!(valsL, -cz)
        end

        rscale = source_scale == :dt2 ? dt^2 : Float64(source_scale)
        push!(rowsR, p); push!(colsR, p); push!(valsR, rscale)
    end

    A_unknown = sparse(rowsA, colsA, valsA, nPoints, nPoints)
    L_known = sparse(rowsL, colsL, valsL, nPoints, 2nPoints)
    R_force = sparse(rowsR, colsR, valsR, nPoints, nPoints * 3)
    b_template = zeros(Float64, nPoints)
    known_lhs_template = zeros(Float64, nPoints, 1, 2)
    known_rhs_template = zeros(Float64, nPoints, 1, 3)

    function b_fun!(b, knownInputs)
        nKnownField = length(known_lhs_template)
        knownFieldVec = @view knownInputs[1:nKnownField]
        knownForceVec = @view knownInputs[nKnownField+1:nKnownField+length(known_rhs_template)]
        b .= 0.0
        mul!(b, L_known, knownFieldVec, -1.0, 1.0)
        mul!(b, R_force, knownForceVec, 1.0, 1.0)
        return b
    end

    function A_fun!(avec, knownInputs)
        avec .= vec(A_unknown)
        return avec
    end

    return (
        is_numerical=true,
        is_fd2d_baseline=true,
        A_unknown=A_unknown,
        L_known=L_known,
        R_force=R_force,
        A_fun! = A_fun!,
        b_fun! = b_fun!,
        A_template=copy(A_unknown),
        b_template=b_template,
        known_lhs_template=known_lhs_template,
        known_rhs_template=known_rhs_template,
        spaceShape=spaceShape,
        NpointsSpace=nPoints,
        NforcePoints=nPoints,
        NField=1,
        NForceField=1,
        timePointsUsedForOneStep=3,
    )
end

function propagate_linear_frames_from_initial(preparedLin, initialPast, initialPresent, Nt; store_every=1, blowup_limit=Inf, diagonal_shift=0.0)
    NField = preparedLin.NField
    NpointsSpace = preparedLin.NpointsSpace
    timePointsUsedForOneStep = preparedLin.timePointsUsedForOneStep
    NknownTime = max(timePointsUsedForOneStep - 1, 0)
    NknownTime == 2 || error("initial helper expects exactly two known time levels; got $NknownTime")

    knownField = zeros(Float64, size(preparedLin.known_lhs_template))
    knownField[:, 1, 1] .= vec(initialPast)
    knownField[:, 1, 2] .= vec(initialPresent)
    knownForce = zero(preparedLin.known_rhs_template)
    unknownField = zeros(Float64, NpointsSpace, NField)

    A = sparse(preparedLin.A_template)
    if diagonal_shift != 0.0
        A = A + diagonal_shift * I
    end
    factor = try
        lu(A)
    catch err
        @error "A_unknown factorization failed" exception=(err, catch_backtrace()) matrix_report=implicit_matrix_report(preparedLin) diagonal_shift=diagonal_shift
        rethrow(err)
    end
    b = copy(preparedLin.b_template)

    frames = Vector{Array{Float64}}()
    push!(frames, copy(initialPresent))
    for it in 1:Nt
        knownInputs = vcat(vec(knownField), vec(knownForce))
        preparedLin.b_fun!(b, knownInputs)
        u = factor \ b
        unknownField .= reshape(real.(u), NpointsSpace, NField)

        umax = maximum(abs, unknownField)
        if !isfinite(umax) || umax > blowup_limit
            @warn "Stopping because wavefield blew up" it umax blowup_limit
            push!(frames, reshape(copy(unknownField[:, 1]), preparedLin.spaceShape...))
            break
        end

        if it % store_every == 0 || it == Nt
            push!(frames, reshape(copy(unknownField[:, 1]), preparedLin.spaceShape...))
        end

        knownField[:, :, 1] .= knownField[:, :, 2]
        knownField[:, :, 2] .= unknownField
    end
    return frames
end


function implicit_matrix_report(preparedLin)
    A = sparse(preparedLin.A_template)
    diagA = abs.(diag(A))
    row_abs = vec(sum(abs.(A); dims=2))
    row_nnz = vec(sum(abs.(A) .> 0; dims=2))
    return (
        size=size(A),
        nnz=nnz(A),
        zero_rows=count(iszero, row_nnz),
        diag_min=isempty(diagA) ? NaN : minimum(diagA),
        diag_median=isempty(diagA) ? NaN : median(diagA),
        diag_max=isempty(diagA) ? NaN : maximum(diagA),
        row_abs_min=minimum(row_abs),
        row_abs_max=maximum(row_abs),
    )
end

function wavefield_snapshot_report(frames)
    rows = NamedTuple[]
    for (i, frame) in pairs(frames)
        finite_vals = Float64.(frame[isfinite.(frame)])
        finite_max = isempty(finite_vals) ? NaN : maximum(abs, finite_vals)
        push!(rows, (
            frame=i,
            nbad=count(!isfinite, frame),
            finite_max=finite_max,
            minimum=isempty(finite_vals) ? NaN : minimum(finite_vals),
            maximum=isempty(finite_vals) ? NaN : maximum(finite_vals),
        ))
    end
    return rows
end

function center_of_mass(frame; power=2, floor_fraction=1e-12)
    W = abs.(Float64.(frame)).^power
    wmax = maximum(W)
    W[W .< floor_fraction * max(wmax, eps(Float64))] .= 0.0
    total = sum(W)
    total == 0.0 && return (x=NaN, z=NaN, total=0.0)
    sx = 0.0; sz = 0.0
    for I in CartesianIndices(frame)
        w = W[I]
        sx += Tuple(I)[1] * w
        sz += Tuple(I)[2] * w
    end
    return (x=sx/total, z=sz/total, total=total)
end

function drift_report(frames, center; power=2, floor_fraction=1e-12)
    c0 = Tuple(center)
    return [begin
        cm = center_of_mass(frame; power=power, floor_fraction=floor_fraction)
        (
            frame=i,
            cm_x=cm.x,
            cm_z=cm.z,
            drift_x=cm.x - c0[1],
            drift_z=cm.z - c0[2],
            maxabs=maximum(abs, frame),
        )
    end for (i, frame) in pairs(frames)]
end

function argmax_report(frames, center)
    c0 = Tuple(center)
    return [begin
        absf = abs.(frame)
        I = CartesianIndices(size(frame))[argmax(absf)]
        (
            frame=i,
            maxpoint=I,
            dx=Tuple(I)[1] - c0[1],
            dz=Tuple(I)[2] - c0[2],
            value=frame[I],
            maxabs=absf[I],
        )
    end for (i, frame) in pairs(frames)]
end

function symmetry_error(frame, center)
    cx, cz = Tuple(center)
    nx, nz = size(frame)
    err_lr = 0.0; den_lr = 0.0
    err_ud = 0.0; den_ud = 0.0
    err_diag = 0.0; den_diag = 0.0
    for I in CartesianIndices(frame)
        i, j = Tuple(I)
        im = 2cx - i
        jm = 2cz - j
        if 1 <= im <= nx
            a = frame[i, j]; b = frame[im, j]
            err_lr += abs(a - b); den_lr += abs(a) + abs(b)
        end
        if 1 <= jm <= nz
            a = frame[i, j]; b = frame[i, jm]
            err_ud += abs(a - b); den_ud += abs(a) + abs(b)
        end
        if 1 <= im <= nx && 1 <= jm <= nz
            a = frame[i, j]; b = frame[im, jm]
            err_diag += abs(a - b); den_diag += abs(a) + abs(b)
        end
    end
    return (
        lr=den_lr == 0 ? 0.0 : err_lr / den_lr,
        ud=den_ud == 0 ? 0.0 : err_ud / den_ud,
        diag=den_diag == 0 ? 0.0 : err_diag / den_diag,
    )
end

function symmetry_report(frames, center)
    return [(frame=i, symmetry_error(frame, center)...) for (i, frame) in pairs(frames)]
end

function plot_wave_snapshots(frames; indices=nothing, sourcePoint=nothing, clim=nothing, ncols=3, title="wave snapshots")
    isempty(frames) && error("frames is empty")
    if indices === nothing
        indices = unique(round.(Int, range(1, length(frames), length=min(9, length(frames)))))
    end
    maxima = [maximum(abs, frames[i]) for i in indices]
    if clim === nothing
        vmax = maximum(maxima)
        clim = vmax == 0 ? 1.0 : vmax
    end
    n = length(indices)
    nrows = cld(n, ncols)
    fig = Figure(size=(320*ncols, 300*nrows))
    for (k, iframe) in enumerate(indices)
        r = cld(k, ncols)
        c = k - (r - 1) * ncols
        ax = Axis(fig[r, c], aspect=DataAspect(), title="frame $iframe")
        heatmap!(ax, frames[iframe]; colormap=:balance, colorrange=(-clim, clim))
        if sourcePoint !== nothing
            xy = Tuple(sourcePoint)
            scatter!(ax, [xy[1]], [xy[2]], color=:red, markersize=8)
        end
    end
    Label(fig[0, :], title)
    return fig
end

function build_toy_opt_prepared(; velocity_value=2600.0, shape=(201,201), dx=100.0, cfl=0.45,
    famousEquationType="2DacousticTime", pointsInSpace=3, pointsInTime=3, supplementaryOrder=2,
    orderBspace=1, orderBtime=1)
    velocity = fill(Float64(velocity_value), shape)
    dt = cfl * dx / velocity_value
    delta = (dx, dx, dt)
    fieldItpl = (ptsSpace=1, ptsTime=1, offsetSpace=1, offsetTime=1, YorderBspace=0, YorderBtime=0)
    materItpl = (ptsSpace=1, ptsTime=1, offsetSpace=1, offsetTime=1, YorderBspace=0, YorderBtime=0)
    params = @strdict famousEquationType Δ=delta orderBtime orderBspace pointsInSpace pointsInTime supplementaryOrder fieldItpl materItpl
    optRec = Main.flexOPT.makeOPTsemiSymbolic(params)
    recette = optRec["recette"]
    modelPoints = Main.flexOPT.getModelPoints(velocity, pointsInTime, recette.numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
    modelFam = (models=[velocity], modelPoints=modelPoints, Δ=delta, modelName="toy_OPT_gaussian_drift")
    numParams = @strdict optRec=optRec modelFam=modelFam absorbingBoundaries=nothing maskedRegionInSpace=nothing representation="matrixfree"
    numOpt = Main.flexOPT.numericalOperatorConstruction(numParams)
    numOps = numOpt["numOperators"]
    prepared = Main.flexOPT.prepareLinearSystem(numOps)
    return (; velocity, delta, optRec, numOps, prepared)
end
