using Pkg
using Metal

cd(@__DIR__)
Pkg.activate("../")
ParamFile = "../config/testparam.csv"  # maybe GeoPoints and planet1D should be fusioned

# batchGPU should be at this level (I have not made it as a module yet, since the choice of Metal/CUDA should be done in a manual way)
include("../src/batchFiles/batchGPU.jl")


include("../src/commonBatchs.jl")
include("../src/planet1D.jl")
include("../src/GeoPoints.jl")

include("../src/flexOPT.jl")

using .commonBatchs, .planet1D, .GeoPoints, .flexOPT


nameConfigs = NamedTuple[]

push!(nameConfigs, (
    name = "convFD3",
    orderBtime = 1,
    orderBspace = -1,
    pointsInSpace = 3,
    pointsInTime = 1,
    supplementaryOrder = 0,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 1.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 1.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "convFD5",
    orderBtime = 1,
    orderBspace = -1,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 0,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "convFD3_sup2",
    orderBtime = 1,
    orderBspace = -1,
    pointsInSpace = 3,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 1.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 1.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "convFD5_sup2",
    orderBtime = 1,
    orderBspace = -1,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT3",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 3,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 1.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 1.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT3_lin",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 3,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 0.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 0.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT3_lin_stag",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 3,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 0.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 2, ptsTime = 1,
        offsetSpace = 0.5, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT4",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 4,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 1.5, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 1.5, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT4_lin",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 4,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 1.5, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 1.5, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT4_lin_stag",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 4,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 1.5, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 4, ptsTime = 1,
        offsetSpace = 0.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT4_2",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 4,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 1.5, offsetTime = 1,
        YorderBspace = 2, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 1.5, offsetTime = 1,
        YorderBspace = 2, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT4_2_stag",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 4,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 3, ptsTime = 1,
        offsetSpace = 1.5, offsetTime = 1,
        YorderBspace = 2, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 4, ptsTime = 1,
        offsetSpace = 0.0, offsetTime = 1,
        YorderBspace = 2, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5_sup4",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 4,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5_lin",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 2,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5_sup4_lin",
    orderBtime = 1,
    orderBspace = 1,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 4,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5_sup4_secon",
    orderBtime = 1,
    orderBspace = 2,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 4,
    fieldItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 1, ptsTime = 1,
        offsetSpace = 2.0, offsetTime = 1,
        YorderBspace = -1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5_sup4_secon_lin",
    orderBtime = 1,
    orderBspace = 2,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 4,
    fieldItpl = (
        ptsSpace = 4, ptsTime = 1,
        offsetSpace = 0.5, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 4, ptsTime = 1,
        offsetSpace = 0.5, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5_sup4_secon_2",
    orderBtime = 1,
    orderBspace = 2,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 4,
    fieldItpl = (
        ptsSpace = 4, ptsTime = 1,
        offsetSpace = 0.5, offsetTime = 1,
        YorderBspace = 2, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 4, ptsTime = 1,
        offsetSpace = 0.5, offsetTime = 1,
        YorderBspace = 2, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5_sup4_secon_lin_stag",
    orderBtime = 1,
    orderBspace = 2,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 4,
    fieldItpl = (
        ptsSpace = 5, ptsTime = 1,
        offsetSpace = 0.0, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 4, ptsTime = 1,
        offsetSpace = 0.5, offsetTime = 1,
        YorderBspace = 1, YorderBtime = 1,
    ),
))

push!(nameConfigs, (
    name = "OPT5_sup4_secon_2_stag",
    orderBtime = 1,
    orderBspace = 2,
    pointsInSpace = 5,
    pointsInTime = 1,
    supplementaryOrder = 4,
    fieldItpl = (
        ptsSpace = 5, ptsTime = 1,
        offsetSpace = 0.0, offsetTime = 1,
        YorderBspace = 2, YorderBtime = 1,
    ),
    materItpl = (
        ptsSpace = 4, ptsTime = 1,
        offsetSpace = 0.5, offsetTime = 1,
        YorderBspace = 2, YorderBtime = 1,
    ),
))


nConfigurations=length(nameConfigs)

using Symbolics,CairoMakie,LinearAlgebra
CairoMakie.activate!()

logsOfHinverse = [0.5*i for i in 0:6] #[1.0*i for i in 0:4]

cases=[]
#prefix="B"*string(tmpOrderBspace)*"_"*"w"*string(tmpWorderBspace)*"_"*string(tmpSupplementaryOrder)*"_"
prefix=""
L = 10.0*π
@variables x
∂ = Differential(x)

cases = push!(cases,(name=prefix*"homogeneous",u=cos(x),β=1.0))
cases = push!(cases,(name=prefix*"sameλ",u=cos(x),β=sin(x)+2))
cases = push!(cases,(name=prefix*"twiceλ",u=cos(x),β=sin(x/2) + 2))
cases = push!(cases,(name=prefix*"sameλ_shifted_π_3",u=cos(x),β=sin(x+π/3) + 2))
cases = push!(cases,(name=prefix*"λ_2",u=cos(x),β=cos(x).^2 + 1))
cases = push!(cases,(name=prefix*"quadratic",u=cos(x),β=x^2+ 1))


misfit = Array{Float64,3}(undef,length(logsOfHinverse),length(cases),nConfigurations)


fig = Figure()
ax = Axis(fig[1, 1], xlabel="x", ylabel="model")
iH=5
for iCase ∈ eachindex(cases)
    name,_,β = cases[iCase]
    ΔxTry = exp(-logsOfHinverse[iH])
    Nx = Int(L÷ΔxTry) +1
    Δx = L/(Nx-1)
    X = [Δx * (i-1) for i ∈ range(1,Nx)]
    model=[Symbolics.value(substitute(β,Dict(x=>X[i]))) for i ∈ range(1,Nx)]
    lines!(ax, X, model)
end
ylims!(ax, 0, 5)
display(fig)
Nx=nothing
Δx=nothing
iH=nothing
modelFamily=Array{Any,3}(undef,length(logsOfHinverse),length(cases),nConfigurations)
forceFamily=Array{Any,3}(undef,length(logsOfHinverse),length(cases),nConfigurations)
for iConfigUsed ∈ eachindex(nameConfigs), iCase ∈ eachindex(cases), iH ∈ eachindex(logsOfHinverse)
    name,T,β = cases[iCase]
    iExperiment = (iH=iH,iCase=iCase,iPointsUsed=iConfigUsed)

    ΔxTry = exp(-logsOfHinverse[iH])
    Nx = Int(L÷ΔxTry) +1
    Δx = L/(Nx-1)
    X = [Δx * (i-1) for i ∈ range(1,Nx)]
    modelName = name*string(Nx)
    models=[]
    model=[Symbolics.value(substitute(β,Dict(x=>X[i]))) for i ∈ range(1,Nx)]
    models=push!(models, model)
    modelPoints = (Nx)
    
    
    q = mySimplify(β*∂(T))
    qₓ = mySimplify(∂(q))
    symbols=(Nx=Nx,Δx=Δx,X=X,q=q,qₓ=qₓ,T=T,β=β)
    
    tmpModel = (models=models, modelName=modelName, modelPoints=modelPoints,Δ=(Δx),symbols=symbols)
    modelFamily[iH,iCase,iConfigUsed]=tmpModel

    force = [Symbolics.value(substitute(qₓ,Dict(x=>X[i]))) for i ∈ range(1,Nx)]
    sourceFull=reshape(force,Nx,1,1)
    forceFamily[iH,iCase,iConfigUsed]=sourceFull
   
end


famousEquationType="1DpoissonHetero" 


using JLD2, LinearAlgebra

checkpoint_file = "misfit_checkpoint_smaller.jld2"
figdir = joinpath(pwd(), "tmp", "figures")
mkpath(figdir)


if isfile(checkpoint_file)
    data = load(checkpoint_file)
    misfit = data["misfit"]
    done = data["done"]
else
    misfit = fill(NaN, length(logsOfHinverse), length(cases), nConfigurations)
    done = falses(length(logsOfHinverse), length(cases), nConfigurations)
end

for iConfigUsed in eachindex(nameConfigs), iCase in eachindex(cases), iH in eachindex(logsOfHinverse)
    if done[iH, iCase, iConfigUsed]
        continue
    end
    #if iH*iCase >1
    #    continue
    #end

    


    try
        name, T, β = cases[iCase]
        iExperiment = (iH=iH, iCase=iCase, iPointsUsed=iConfigUsed)

        @unpack Nx, Δx, X, q, qₓ, T, β = modelFamily[iH, iCase, iConfigUsed].symbols
        cfg = nameConfigs[iConfigUsed]
        configNameTmp=cfg.name
        @unpack orderBtime, orderBspace, pointsInSpace, pointsInTime, supplementaryOrder, fieldItpl, materItpl = cfg
        #@unpack orderBtime, orderBspace, pointsInSpace, pointsInTime, supplementaryOrder, fieldItpl, materItpl = all_Configcases[iConfigUsed]

        concreteParametersForOPTConstruction = @strdict famousEquationType Δ=modelFamily[iH,iCase,iConfigUsed].Δ orderBtime orderBspace pointsInSpace pointsInTime supplementaryOrder fieldItpl materItpl
        optRec = myProduceOrLoad(makeOPTsemiSymbolic, concreteParametersForOPTConstruction, "semiSymbolic")

        params = @strdict optRec=optRec modelFam=modelFamily[iExperiment.iH,iExperiment.iCase,iExperiment.iPointsUsed] absorbingBoundaries=nothing maskedRegionInSpace=nothing
        numOpt = numericalOperatorConstruction(params)

        numOps = numOpt["numOperators"]
        #prepared = prepareNumericalOperators(numOps)
        preparedLin = prepareLinearSystem(numOps)


        sourceFull = forceFamily[iExperiment.iH, iExperiment.iCase, iExperiment.iPointsUsed]        

        syntheticData = timeMarchingSchemeLinear(
            preparedLin,
            1,
            modelFamily[iExperiment.iH,iExperiment.iCase,iExperiment.iPointsUsed].Δ,
            modelFamily[iExperiment.iH,iExperiment.iCase,iExperiment.iPointsUsed].modelName;
            videoMode=false,
            sourceType="Explicit",
            sourceFull=sourceFull,
            iExperiment=iExperiment,
            boundaryConditionForced=true,
        )



        #syntheticData = timeMarchingSchemePrepared(
        #    prepared,
        #    1,
        #    modelFamily[iExperiment.iH,iExperiment.iCase,iExperiment.iPointsUsed].Δ,
        #    modelFamily[iExperiment.iH,iExperiment.iCase,iExperiment.iPointsUsed].modelName;
        #    videoMode=false,
        #    sourceType="Explicit",
        #    sourceFull=sourceFull,
        #    iExperiment=iExperiment,
        #    boundaryConditionForced=true,
        #)

        syntheticData = reduce(vcat, syntheticData)
        analyticalData = [Symbolics.value(substitute(T, Dict(x => X[i]))) for i in 1:Nx]

        misfit[iH, iCase, iConfigUsed] = norm(syntheticData - analyticalData) / Nx
        done[iH, iCase, iConfigUsed] = true

        fig = Figure()
        ax = Axis(
            fig[1, 1],
            title = "model=$(cases[iCase].name), $(configNameTmp)",
            xlabel = "x",
            ylabel = "solution",
        )
        lines!(ax, X, analyticalData, color=:blue, label="analytical")
        scatter!(ax, X, syntheticData, color=:red, marker=:circle, label="synthetic")
        axislegend(ax)

        figfile = joinpath(
            figdir,
            "cmp_iH$(iH)_iCase$(iCase)_iConfig$(iConfigUsed)_obs$(orderBspace)_obt$(orderBtime)_pts$(pointsInSpace)_supp$(supplementaryOrder).png",
        )
        save(figfile, fig)

        jldsave(checkpoint_file; misfit=misfit, done=done)

    catch err
        @warn "Failed at (iH=$iH, iCase=$iCase, iConfigUsed=$iConfigUsed)" exception=(err, catch_backtrace())
        jldsave(checkpoint_file; misfit=misfit, done=done)
    end
end