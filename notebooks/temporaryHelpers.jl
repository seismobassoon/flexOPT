
using LinearAlgebra
using SparseArrays
using Statistics
using DrWatson
using KernelAbstractions
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
        @error "A_unknown factorization failed" exception=(err, catch_backtrace()) matrix_report=implicit_matrix_report(A) diagonal_shift=diagonal_shift
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


function implicit_matrix_report(A::SparseMatrixCSC)
    A = sparse(A)
    nrow, ncol = size(A)
    rows, cols, vals = findnz(A)
    finite_mask = isfinite.(vals)
    finite_vals = vals[finite_mask]

    row_finite_nnz = zeros(Int, nrow)
    row_nonfinite_nnz = zeros(Int, nrow)
    row_abs_finite = zeros(Float64, nrow)
    for (r, v) in zip(rows, vals)
        if isfinite(v)
            row_finite_nnz[r] += 1
            row_abs_finite[r] += abs(Float64(v))
        else
            row_nonfinite_nnz[r] += 1
        end
    end

    diag_vals = diag(A)
    finite_diag_abs = abs.(Float64.(diag_vals[isfinite.(diag_vals)]))
    finite_abs = abs.(Float64.(finite_vals))
    return (
        size=size(A),
        nnz=nnz(A),
        stored_finite=length(finite_vals),
        stored_nonfinite=count(!isfinite, vals),
        stored_nan=count(isnan, vals),
        stored_inf=count(isinf, vals),
        zero_rows_finite=count(iszero, row_finite_nnz),
        rows_with_nonfinite=count(>(0), row_nonfinite_nnz),
        finite_abs_min=isempty(finite_abs) ? NaN : minimum(finite_abs),
        finite_abs_median=isempty(finite_abs) ? NaN : median(finite_abs),
        finite_abs_max=isempty(finite_abs) ? NaN : maximum(finite_abs),
        diag_finite_count=length(finite_diag_abs),
        diag_nonfinite_count=count(!isfinite, diag_vals),
        diag_abs_min=isempty(finite_diag_abs) ? NaN : minimum(finite_diag_abs),
        diag_abs_median=isempty(finite_diag_abs) ? NaN : median(finite_diag_abs),
        diag_abs_max=isempty(finite_diag_abs) ? NaN : maximum(finite_diag_abs),
        finite_row_abs_min=isempty(row_abs_finite) ? NaN : minimum(row_abs_finite),
        finite_row_abs_max=isempty(row_abs_finite) ? NaN : maximum(row_abs_finite),
    )
end

implicit_matrix_report(preparedLin) = implicit_matrix_report(sparse(preparedLin.A_template))

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

function build_toy_opt_prepared(; velocity_value=2600.0, shape=(201,201), dx=100.0, cfl=0.45, dt=nothing,
    famousEquationType="2DacousticTime", pointsInSpace=3, pointsInTime=3, supplementaryOrder=2,
    orderBspace=1, orderBtime=1, YorderBspace=-1, YorderBtime=-1, recipe_backend=CPU())
    velocity = fill(Float64(velocity_value), shape)
    dt_value = dt === nothing ? cfl * dx / velocity_value : Float64(dt)
    delta = (dx, dx, dt_value)
    fieldItpl = (ptsSpace=1, ptsTime=1, offsetSpace=1, offsetTime=1, YorderBspace=YorderBspace, YorderBtime=YorderBtime)
    materItpl = (ptsSpace=1, ptsTime=1, offsetSpace=1, offsetTime=1, YorderBspace=YorderBspace, YorderBtime=YorderBtime)
    params = @strdict famousEquationType Δ=delta orderBtime orderBspace pointsInSpace pointsInTime supplementaryOrder fieldItpl materItpl
    old_backend = Main.flexOPT.backend
    if recipe_backend !== nothing
        Main.flexOPT.backend = recipe_backend
    end
    optRec = try
        Main.flexOPT.makeOPTsemiSymbolic(params)
    finally
        Main.flexOPT.backend = old_backend
    end
    recette = optRec["recette"]
    modelPoints = Main.flexOPT.getModelPoints(velocity, pointsInTime, recette.numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
    material_model = occursin("Homo", famousEquationType) ? Float64(velocity_value) : velocity
    modelFam = (models=[material_model], modelPoints=modelPoints, Δ=delta, modelName="toy_OPT_gaussian_drift")
    numParams = @strdict optRec=optRec modelFam=modelFam absorbingBoundaries=nothing maskedRegionInSpace=nothing representation="matrixfree"
    numOpt = Main.flexOPT.numericalOperatorConstruction(numParams)
    numOps = numOpt["numOperators"]
    prepared = Main.flexOPT.prepareLinearSystem(numOps)
    return (; velocity, delta, optRec, numOps, prepared)
end


function time_role_from_slot(iT, activeTimePoints)
    iT == activeTimePoints && return :future
    iT == activeTimePoints - 1 && return :present
    return Symbol("past_", activeTimePoints - iT)
end

function recipe_local_index_audit(optRec; geometry=1)
    recette = optRec isa AbstractDict ? optRec["recette"] : optRec.recette
    nodes = recette.nodes[geometry]
    centre_linear = recette.centresIndices[geometry]
    centre_point = nodes[centre_linear]
    activeTimePoints = size(nodes, ndims(nodes))
    LI = LinearIndices(nodes)
    rows = NamedTuple[]
    for I in CartesianIndices(nodes)
        p = nodes[I]
        offset = Tuple(Int.(p .- centre_point))
        iT = Int(p[end])
        push!(rows, (
            linear=LI[I],
            cartesian_index=Tuple(I),
            local_point=Tuple(Int.(p)),
            offset=offset,
            spatial_offset=offset[1:end-1],
            time_slot=iT,
            time_role=time_role_from_slot(iT, activeTimePoints),
        ))
    end
    return rows
end

function operator_stencil_at_point(op, point::CartesianIndex; iExpr=1, iField=1, atol=0.0)
    geometry = op.geometry
    spaceShape = Tuple(geometry.wholeRegionPointsSpace)
    NpointsSpace = length(geometry.νWhole)
    activeTimePoints = geometry.activeTimePoints
    Nexpr = div(op.size[1], NpointsSpace)
    Nfield = div(op.size[2], NpointsSpace * activeTimePoints)
    pointLinear = LinearIndices(spaceShape)
    fieldTimeSpace = CartesianIndices((Nfield, activeTimePoints, spaceShape...))
    row = LinearIndices((Nexpr, NpointsSpace))[iExpr, pointLinear[point]]

    rows = NamedTuple[]
    for k in eachindex(op.table.vals)
        Int(op.table.rows[k]) == row || continue
        ci = fieldTimeSpace[Int(op.table.cols[k])]
        Int(ci[1]) == iField || continue
        coef = op.table.vals[k]
        abs(coef) > atol || continue
        jPoint = CartesianIndex(Tuple(ci)[3:end])
        offset = Tuple(jPoint - point)
        iT = Int(ci[2])
        push!(rows, (
            time_slot=iT,
            time_role=time_role_from_slot(iT, activeTimePoints),
            offset=offset,
            coef=coef,
            abscoef=abs(coef),
            col=Int(op.table.cols[k]),
        ))
    end
    sort!(rows, by = r -> (r.time_slot, r.offset...))
    return rows
end

function operator_stencil_at_point(numOps::NamedTuple, point::CartesianIndex; which=:left, iExpr=1, iField=1, atol=0.0)
    op = which == :left ? numOps.numericalOperators.left :
         which == :right ? numOps.numericalOperators.right :
         which == :residual ? numOps.numericalOperators.residual :
         error("which must be :left, :right, or :residual")
    return operator_stencil_at_point(op, point; iExpr=iExpr, iField=iField, atol=atol)
end

function stencil_time_summary(stencil_rows)
    out = NamedTuple[]
    for iT in sort(unique(r.time_slot for r in stencil_rows))
        rs = filter(r -> r.time_slot == iT, stencil_rows)
        push!(out, (
            time_slot=iT,
            time_role=rs[1].time_role,
            n=length(rs),
            sumcoef=sum(r.coef for r in rs),
            sumabs=sum(r.abscoef for r in rs),
            maxabs=maximum(r.abscoef for r in rs),
        ))
    end
    return out
end


function stencil_matrices_by_time(stencil_rows; radius=1)
    slots = sort(unique(r.time_slot for r in stencil_rows))
    mats = NamedTuple[]
    for iT in slots
        rs = filter(r -> r.time_slot == iT, stencil_rows)
        M = zeros(Float64, 2radius + 1, 2radius + 1)
        for r in rs
            ox, oy = r.offset
            if -radius <= ox <= radius && -radius <= oy <= radius
                M[ox + radius + 1, oy + radius + 1] = Float64(real(r.coef))
            end
        end
        push!(mats, (time_slot=iT, time_role=rs[1].time_role, matrix=M))
    end
    return mats
end

function fd2d_acoustic_reference_stencil(v, dx, dt)
    cx = (v / dx)^2
    ct = 1 / dt^2
    mats = [zeros(Float64, 3, 3) for _ in 1:3]
    # residual: u_tt - v^2(u_xx + u_yy)
    mats[1][2,2] = ct
    mats[2][2,2] = -2ct + 4cx
    mats[2][1,2] = -cx
    mats[2][3,2] = -cx
    mats[2][2,1] = -cx
    mats[2][2,3] = -cx
    mats[3][2,2] = ct
    return [
        (time_slot=1, time_role=:past_2, matrix=mats[1]),
        (time_slot=2, time_role=:present, matrix=mats[2]),
        (time_slot=3, time_role=:future, matrix=mats[3]),
    ]
end

function compare_stencil_scale_to_fd(stencil_rows, v, dx, dt; radius=1)
    opt = stencil_matrices_by_time(stencil_rows; radius=radius)
    fd = fd2d_acoustic_reference_stencil(v, dx, dt)
    rows = NamedTuple[]
    for (o, f) in zip(opt, fd)
        optmax = maximum(abs, o.matrix)
        fdmax = maximum(abs, f.matrix)
        push!(rows, (
            time_slot=o.time_slot,
            time_role=o.time_role,
            opt_maxabs=optmax,
            fd_maxabs=fdmax,
            opt_over_fd=fdmax == 0 ? NaN : optmax / fdmax,
            opt_sum=sum(o.matrix),
            fd_sum=sum(f.matrix),
        ))
    end
    return rows
end


function taylor_C_moment_report(optRec, delta; side=:lhs, geometry=1, mu_index=1)
    recette = optRec isa AbstractDict ? optRec["recette"] : optRec.recette
    C = side == :lhs ? recette.Cˡη : recette.CˡηForce
    nodes = recette.nodes[geometry]
    center = nodes[recette.centresIndices[geometry]]
    nEta, nL, nMu = size(C)
    nDim = length(center)
    Lside = round(Int, nL^(1 / nDim))
    prod(fill(Lside, nDim)) == nL || error("cannot infer tensor Taylor order shape from nL=$nL and nDim=$nDim")
    Linds = CartesianIndices(ntuple(_ -> Lside, nDim))

    A = zeros(Float64, nL, nEta)
    for (iEta, p) in enumerate(vec(nodes))
        dist = (Float64.(p .- center)) .* Float64.(delta)
        for J in Linds
            linJ = LinearIndices(Linds)[J]
            orders = Tuple(J) .- 1
            A[linJ, iEta] = prod(dist .^ orders) / prod(factorial.(orders))
        end
    end

    interesting = [ntuple(_ -> 1, nDim)]
    for d in 1:nDim
        idx = ntuple(i -> i == d ? 3 : 1, nDim)
        push!(interesting, idx)
    end

    rows = NamedTuple[]
    for idx in interesting
        lin = LinearIndices(Linds)[CartesianIndex(idx)]
        weights = C[:, lin, mu_index]
        moments = A * weights
        push!(rows, (
            order=Tuple(idx .- 1),
            linear_order=lin,
            weights_maxabs=maximum(abs, weights),
            target_moment=moments[lin],
            moments_maxabs=maximum(abs, moments),
            leakage_maxabs=maximum(abs, [moments[i] for i in eachindex(moments) if i != lin]),
        ))
    end
    return rows
end



function prepare_fd2d_acoustic_fd5_baseline(velocity, delta; source_scale=:dt2)
    spaceShape = size(velocity)
    length(spaceShape) == 2 || error("FD5 acoustic helper currently expects a 2D velocity model")
    dx, dz, dt = Float64.(delta)
    nPoints = prod(spaceShape)
    LI = LinearIndices(spaceShape)

    rowsA = Int[]; colsA = Int[]; valsA = Float64[]
    rowsL = Int[]; colsL = Int[]; valsL = Float64[]
    rowsR = Int[]; colsR = Int[]; valsR = Float64[]
    col_known(point, it) = point + (it - 1) * nPoints
    wx = Dict(-2 => -1 / (12dx^2), -1 => 4 / (3dx^2), 0 => -5 / (2dx^2), 1 => 4 / (3dx^2), 2 => -1 / (12dx^2))
    wz = Dict(-2 => -1 / (12dz^2), -1 => 4 / (3dz^2), 0 => -5 / (2dz^2), 1 => 4 / (3dz^2), 2 => -1 / (12dz^2))

    for I in CartesianIndices(spaceShape)
        p = LI[I]
        vdt2 = Float64(velocity[I])^2 * dt^2
        i, j = Tuple(I)

        push!(rowsA, p); push!(colsA, p); push!(valsA, 1.0)
        push!(rowsL, p); push!(colsL, col_known(p, 1)); push!(valsL, 1.0)
        push!(rowsL, p); push!(colsL, col_known(p, 2)); push!(valsL, -2.0 - vdt2 * (wx[0] + wz[0]))

        for ox in (-2, -1, 1, 2)
            ii = i + ox
            1 <= ii <= spaceShape[1] || continue
            q = LI[CartesianIndex(ii, j)]
            push!(rowsL, p); push!(colsL, col_known(q, 2)); push!(valsL, -vdt2 * wx[ox])
        end
        for oz in (-2, -1, 1, 2)
            jj = j + oz
            1 <= jj <= spaceShape[2] || continue
            q = LI[CartesianIndex(i, jj)]
            push!(rowsL, p); push!(colsL, col_known(q, 2)); push!(valsL, -vdt2 * wz[oz])
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
        is_fd2d_acoustic_fd5=true,
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

function prepare_fd2d_elastic_homogeneous_baseline(rho, lambda, mu, delta; source_scale=:dt2)
    spaceShape = size(rho)
    size(lambda) == spaceShape || error("lambda should have size $spaceShape")
    size(mu) == spaceShape || error("mu should have size $spaceShape")
    length(spaceShape) == 2 || error("elastic FD helper currently expects 2D arrays")
    maximum(abs, rho .- rho[1]) <= sqrt(eps(Float64)) * max(abs(rho[1]), 1.0) || error("this helper is for homogeneous rho")
    maximum(abs, lambda .- lambda[1]) <= sqrt(eps(Float64)) * max(abs(lambda[1]), 1.0) || error("this helper is for homogeneous lambda")
    maximum(abs, mu .- mu[1]) <= sqrt(eps(Float64)) * max(abs(mu[1]), 1.0) || error("this helper is for homogeneous mu")

    dx, dz, dt = Float64.(delta)
    ρ = Float64(rho[1]); λ = Float64(lambda[1]); μ = Float64(mu[1])
    ax = dt^2 / ρ
    nPoints = prod(spaceShape)
    NField = 2
    LI = LinearIndices(spaceShape)
    rowidx(field, point) = point + (field - 1) * nPoints
    col_known(point, field, it) = point + (field - 1) * nPoints + (it - 1) * nPoints * NField
    col_force(point, field, it) = point + (field - 1) * nPoints + (it - 1) * nPoints * NField

    rowsA = Int[]; colsA = Int[]; valsA = Float64[]
    rowsL = Int[]; colsL = Int[]; valsL = Float64[]
    rowsR = Int[]; colsR = Int[]; valsR = Float64[]
    addL(row, point, field, it, val) = (push!(rowsL, row); push!(colsL, col_known(point, field, it)); push!(valsL, val))
    addR(row, point, field, it, val) = (push!(rowsR, row); push!(colsR, col_force(point, field, it)); push!(valsR, val))

    for I in CartesianIndices(spaceShape)
        p = LI[I]
        i, j = Tuple(I)
        for fld in 1:NField
            r = rowidx(fld, p)
            push!(rowsA, r); push!(colsA, r); push!(valsA, 1.0)
            addL(r, p, fld, 1, 1.0)
            addL(r, p, fld, 2, -2.0)
            rscale = source_scale == :dt2 ? dt^2 / ρ : Float64(source_scale)
            addR(r, p, fld, 2, rscale)
        end

        ux_row = rowidx(1, p)
        uz_row = rowidx(2, p)
        for (di, coeff) in ((-1, 1/dx^2), (0, -2/dx^2), (1, 1/dx^2))
            ii = i + di
            1 <= ii <= spaceShape[1] || continue
            q = LI[CartesianIndex(ii, j)]
            addL(ux_row, q, 1, 2, -ax * (λ + 2μ) * coeff)
            addL(uz_row, q, 2, 2, -ax * μ * coeff)
        end
        for (dj, coeff) in ((-1, 1/dz^2), (0, -2/dz^2), (1, 1/dz^2))
            jj = j + dj
            1 <= jj <= spaceShape[2] || continue
            q = LI[CartesianIndex(i, jj)]
            addL(ux_row, q, 1, 2, -ax * μ * coeff)
            addL(uz_row, q, 2, 2, -ax * (λ + 2μ) * coeff)
        end
        for di in (-1, 1), dj in (-1, 1)
            ii = i + di; jj = j + dj
            1 <= ii <= spaceShape[1] && 1 <= jj <= spaceShape[2] || continue
            q = LI[CartesianIndex(ii, jj)]
            coeff = di * dj / (4dx * dz)
            addL(ux_row, q, 2, 2, -ax * (λ + μ) * coeff)
            addL(uz_row, q, 1, 2, -ax * (λ + μ) * coeff)
        end
    end

    nRows = nPoints * NField
    A_unknown = sparse(rowsA, colsA, valsA, nRows, nRows)
    L_known = sparse(rowsL, colsL, valsL, nRows, nRows * 2)
    R_force = sparse(rowsR, colsR, valsR, nRows, nRows * 3)
    b_template = zeros(Float64, nRows)
    known_lhs_template = zeros(Float64, nPoints, NField, 2)
    known_rhs_template = zeros(Float64, nPoints, NField, 3)

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
        is_fd2d_elastic_homogeneous=true,
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
        NField=NField,
        NForceField=NField,
        timePointsUsedForOneStep=3,
    )
end

function frame_difference_report(frames_a, frames_b)
    n = min(length(frames_a), length(frames_b))
    rows = NamedTuple[]
    for k in 1:n
        A = frames_a[k]
        B = frames_b[k]
        d = A .- B
        push!(rows, (
            frame=k,
            max_a=maximum(abs, A),
            max_b=maximum(abs, B),
            maxdiff=maximum(abs, d),
            reldiff=maximum(abs, d) / max(maximum(abs, A), maximum(abs, B), eps(Float64)),
        ))
    end
    return rows
end

function timed_propagation(preparedLin, Nt; kwargs...)
    elapsed = @elapsed frames = propagate_linear_frames_with_source(preparedLin, Nt; kwargs...)
    return (; elapsed, frames, report=wavefield_snapshot_report(component_frames(frames, 1)))
end

function source_time_samples(Nt, dt, timePointsUsedForOneStep; wavelet=nothing, t0=12dt, f0=0.04)
    ntime = Nt + timePointsUsedForOneStep - 1
    t = (0:ntime-1) .* Float64(dt)
    signal = wavelet === nothing ? Main.flexOPT.Ricker.(t, Float64(t0), Float64(f0)) : wavelet.(t)
    return Float64.(signal)
end

function make_source_full(preparedLin, weights::AbstractVector, timeSignal::AbstractVector; iForceField=1, amplitude=1.0)
    length(weights) == preparedLin.NforcePoints || error("weights length should be NforcePoints=$(preparedLin.NforcePoints)")
    1 <= iForceField <= preparedLin.NForceField || error("iForceField should be between 1 and $(preparedLin.NForceField)")
    sourceFull = zeros(Float64, preparedLin.NforcePoints, preparedLin.NForceField, length(timeSignal))
    for it in eachindex(timeSignal)
        sourceFull[:, iForceField, it] .= amplitude .* weights .* Float64(timeSignal[it])
    end
    return sourceFull
end

function point_source_full(preparedLin, point::CartesianIndex, timeSignal; iForceField=1, amplitude=1.0, normalise=:none)
    weights = point_source_weights(preparedLin, point; normalise=normalise)
    return make_source_full(preparedLin, weights, timeSignal; iForceField=iForceField, amplitude=amplitude)
end

function propagate_linear_frames_with_source(preparedLin, Nt; initialPast=nothing, initialPresent=nothing,
    sourceFull=nothing, store_every=1, blowup_limit=Inf, diagonal_shift=0.0)
    NField = preparedLin.NField
    NpointsSpace = preparedLin.NpointsSpace
    NForceField = preparedLin.NForceField
    timePointsUsedForOneStep = preparedLin.timePointsUsedForOneStep
    NknownTime = max(timePointsUsedForOneStep - 1, 0)
    NknownTime == 2 || error("source propagation helper expects exactly two known time levels; got $NknownTime")

    zeroFrame = zeros(Float64, preparedLin.spaceShape..., NField)
    initialPast = initialPast === nothing ? zeroFrame : Float64.(initialPast)
    initialPresent = initialPresent === nothing ? zeroFrame : Float64.(initialPresent)
    size(initialPast) == size(zeroFrame) || error("initialPast should have size $(size(zeroFrame))")
    size(initialPresent) == size(zeroFrame) || error("initialPresent should have size $(size(zeroFrame))")

    if sourceFull === nothing
        sourceFull = zeros(Float64, preparedLin.NforcePoints, NForceField, Nt + timePointsUsedForOneStep - 1)
    end
    size(sourceFull, 1) == preparedLin.NforcePoints || error("sourceFull first dimension should be NforcePoints=$(preparedLin.NforcePoints)")
    size(sourceFull, 2) == NForceField || error("sourceFull second dimension should be NForceField=$NForceField")
    size(sourceFull, 3) >= Nt + timePointsUsedForOneStep - 1 || error("sourceFull has too few time samples")

    knownField = zeros(Float64, size(preparedLin.known_lhs_template))
    knownField[:, :, 1] .= reshape(initialPast, NpointsSpace, NField)
    knownField[:, :, 2] .= reshape(initialPresent, NpointsSpace, NField)
    knownForce = zero(preparedLin.known_rhs_template)
    unknownField = zeros(Float64, NpointsSpace, NField)

    A = sparse(preparedLin.A_template)
    if diagonal_shift != 0.0
        A = A + diagonal_shift * I
    end
    factor = lu(A)
    b = copy(preparedLin.b_template)

    frames = Vector{Array{Float64}}()
    push!(frames, copy(initialPresent))
    for it in 1:Nt
        knownForce .= sourceFull[:, :, it:it+timePointsUsedForOneStep-1]
        knownInputs = vcat(vec(knownField), vec(knownForce))
        preparedLin.b_fun!(b, knownInputs)
        u = factor \ b
        unknownField .= reshape(real.(u), NpointsSpace, NField)

        umax = maximum(abs, unknownField)
        if !isfinite(umax) || umax > blowup_limit
            @warn "Stopping because wavefield blew up" it umax blowup_limit
            push!(frames, reshape(copy(unknownField), preparedLin.spaceShape..., NField))
            break
        end

        if it % store_every == 0 || it == Nt
            push!(frames, reshape(copy(unknownField), preparedLin.spaceShape..., NField))
        end

        knownField[:, :, 1] .= knownField[:, :, 2]
        knownField[:, :, 2] .= unknownField
    end
    return frames
end

function component_frames(frames, iField::Integer=1)
    return [Array(selectdim(frame, ndims(frame), iField)) for frame in frames]
end

function cfl_diagnostics(model_velocity, delta; cfl_safety=0.45, ppw=10, samples_per_period=20)
    v = Float64.(vec(model_velocity))
    v = v[isfinite.(v) .& (v .> 0)]
    dxmin = minimum(Float64.(delta[1:end-1]))
    vmax = maximum(v)
    vmed = median(v)
    suggested_dt_2D = cfl_safety * dxmin / (sqrt(2) * vmax)
    suggested_f0 = vmed / (ppw * dxmin)
    suggested_dt_from_f0 = 1 / (samples_per_period * suggested_f0)
    return (; vmax, vmed, dxmin, current_dt=Float64(delta[end]), suggested_dt_2D, suggested_f0, suggested_dt_from_f0)
end

function build_opt_prepared(famousEquationType, models, delta; pointsInSpace=3, pointsInTime=3,
    supplementaryOrder=2, orderBspace=1, orderBtime=1, YorderBspace=-1, YorderBtime=-1,
    modelName="OPT_model", recipe_backend=CPU())
    fieldItpl = (ptsSpace=1, ptsTime=1, offsetSpace=1, offsetTime=1, YorderBspace=YorderBspace, YorderBtime=YorderBtime)
    materItpl = (ptsSpace=1, ptsTime=1, offsetSpace=1, offsetTime=1, YorderBspace=YorderBspace, YorderBtime=YorderBtime)
    params = @strdict famousEquationType Δ=delta orderBtime orderBspace pointsInSpace pointsInTime supplementaryOrder fieldItpl materItpl
    old_backend = Main.flexOPT.backend
    if recipe_backend !== nothing
        Main.flexOPT.backend = recipe_backend
    end
    optRec = try
        Main.flexOPT.makeOPTsemiSymbolic(params)
    finally
        Main.flexOPT.backend = old_backend
    end
    recette = optRec["recette"]
    modelPoints = Main.flexOPT.getModelPoints(models[1], pointsInTime, recette.numbersOfTheSystem.numbersOfTheSystemL.timeMarching)
    modelFam = (models=models, modelPoints=modelPoints, Δ=delta, modelName=modelName)
    numParams = @strdict optRec=optRec modelFam=modelFam absorbingBoundaries=nothing maskedRegionInSpace=nothing representation="matrixfree"
    numOpt = Main.flexOPT.numericalOperatorConstruction(numParams)
    numOps = numOpt["numOperators"]
    prepared = Main.flexOPT.prepareLinearSystem(numOps)
    return (; optRec, numOps, prepared, modelFam)
end


function source_rhs_diagnostics(preparedLin, sourceFull; it=1, initialPast=nothing, initialPresent=nothing)
    NField = preparedLin.NField
    NpointsSpace = preparedLin.NpointsSpace
    NForceField = preparedLin.NForceField
    timePointsUsedForOneStep = preparedLin.timePointsUsedForOneStep
    NknownTime = max(timePointsUsedForOneStep - 1, 0)

    zeroFrame = zeros(Float64, preparedLin.spaceShape..., NField)
    initialPast = initialPast === nothing ? zeroFrame : Float64.(initialPast)
    initialPresent = initialPresent === nothing ? zeroFrame : Float64.(initialPresent)

    size(sourceFull, 1) == preparedLin.NforcePoints || error("sourceFull first dimension should be NforcePoints=$(preparedLin.NforcePoints)")
    size(sourceFull, 2) == NForceField || error("sourceFull second dimension should be NForceField=$NForceField")
    size(sourceFull, 3) >= it + timePointsUsedForOneStep - 1 || error("sourceFull has too few time samples for it=$it")

    knownField = zeros(Float64, size(preparedLin.known_lhs_template))
    if NknownTime >= 1
        knownField[:, :, 1] .= reshape(initialPast, NpointsSpace, NField)
    end
    if NknownTime >= 2
        knownField[:, :, 2] .= reshape(initialPresent, NpointsSpace, NField)
    end
    knownForce = zero(preparedLin.known_rhs_template)
    knownForce .= sourceFull[:, :, it:it+timePointsUsedForOneStep-1]

    b = copy(preparedLin.b_template)
    preparedLin.b_fun!(b, vcat(vec(knownField), vec(knownForce)))
    bmat = reshape(real.(b), preparedLin.spaceShape..., NField)
    return (
        it=it,
        force_max=maximum(abs, knownForce),
        force_sum=sum(knownForce),
        force_nonzero=count(!iszero, knownForce),
        b_norm=norm(b),
        b_max=maximum(abs, b),
        b_sum=sum(b),
        b_nonzero=count(!iszero, b),
        b_argmax=CartesianIndex(Tuple(argmax(abs.(bmat)))[1:ndims(bmat)-1]),
        b_frame=bmat,
    )
end

function source_rhs_scan(preparedLin, sourceFull; its=1:min(10, size(sourceFull, 3) - preparedLin.timePointsUsedForOneStep + 1))
    rows = NamedTuple[]
    for it in its
        d = source_rhs_diagnostics(preparedLin, sourceFull; it=it)
        push!(rows, (
            it=it,
            force_max=d.force_max,
            force_nonzero=d.force_nonzero,
            b_norm=d.b_norm,
            b_max=d.b_max,
            b_sum=d.b_sum,
            b_nonzero=d.b_nonzero,
            b_argmax=d.b_argmax,
        ))
    end
    return rows
end

function rhs_stencil_at_source(numOps, point::CartesianIndex; iExpr=1, iForceField=1, atol=0.0)
    return operator_stencil_at_point(numOps, point; which=:right, iExpr=iExpr, iField=iForceField, atol=atol)
end

function elastic_lame_from_rho_vp_vs(rho, vp, vs; rho_scale=1000.0, velocity_scale=1000.0)
    rho_mks = Float64.(rho) .* rho_scale
    vp_mks = Float64.(vp) .* velocity_scale
    vs_mks = Float64.(vs) .* velocity_scale
    mu = rho_mks .* vs_mks.^2
    lambda = rho_mks .* vp_mks.^2 .- 2 .* mu
    return rho_mks, lambda, mu, vp_mks, vs_mks
end

function downsample_center_crop(A, shape::Tuple{Int,Int}; step=1)
    B = A[1:step:end, 1:step:end]
    nx, nz = size(B)
    sx, sz = shape
    ix0 = max(1, cld(nx - sx, 2) + 1)
    iz0 = max(1, cld(nz - sz, 2) + 1)
    return B[ix0:min(ix0+sx-1, nx), iz0:min(iz0+sz-1, nz)]
end
