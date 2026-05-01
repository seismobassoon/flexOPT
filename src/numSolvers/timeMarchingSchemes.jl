
using DrWatson,JLD2
using GLMakie: Figure, Axis, heatmap!, Colorbar, record

function timeMarchingScheme(opt, Nt, Δnum,modelName;videoMode=true,sourceType="Ricker",t₀=50,f₀=0.04,initialCondition=0.0,sourceFull=nothing, iExperiment=nothing)

    #region reading data and source time functions
    if !isdir(datadir("fieldResults"))
        mkdir(datadir("fieldResults"))
    end
    @unpack iH, iCase, iPointsUsed = iExperiment
    experiment_name = "$(iH)_$(iCase)_$(iPointsUsed)"
    @show sequentialFileName=datadir("fieldResults", savename((Nt,Δnum...,modelName,sourceType,experiment_name),"jld2"))
    # filename that can be given by DrWatson
    @show compactFileName=datadir("fieldResults", savename("compact",(Nt,Δnum...,modelName,sourceType,experiment_name),"jld2"))
 

    operators = opt["numOperators"]
    costfunctions,fieldLHS,fieldRHS,champsLimité = operators
    costfunctions=reduce(vcat,costfunctions) # we know that costfunctions is a 2D array

    #@show size(costfunctions),size(fieldLHS[1,1]),size(fieldRHS)

    timePointsUsedForOneStep = size(fieldLHS)[2]
    NField = size(fieldLHS)[1]

    pointsFieldSpace = size(fieldLHS[1,1]) # the list of points in space which is useful for reshaping
    NpointsSpace = length(fieldLHS[1,1]) # the full number of points

    itVec=collect(1:1:Nt)
    t=(itVec.-1).*Δnum[end] # time vector # if it's not time marching t will give you just 0.0 regardless of Δnum[end]
    sourceTime = nothing
    #@show costfunctions[45] # check if we can see the source terms


    if sourceType === "Ricker"
        # making a small source

        # if you want to plot:
        # myRicker(x)=Ricker(x,20,0.03)
        # scatter(0:1:100,myRicker) or lines 

       
        myRicker(x)=Ricker(x,t₀,f₀)
        
        sourceTime = myRicker.(t)
    else sourceType == "Explicit"
        #@show size(sourceFull), (pointsFieldSpace...,NField,Nt)
        if size(sourceFull) !== (pointsFieldSpace...,NField,Nt)
            @error "Oh, this options really demands you a full description of the force for space and time for the points ! Cheers!"
        else
            sourceFull=reshape(sourceFull,NpointsSpace,NField,Nt)
        end
    end

    #endregion

    #region 1Dvectorisation of knowns and unknowns
    
    symbKnownField = Array{Num,3}(undef,NpointsSpace,NField,timePointsUsedForOneStep-1)
    knownField = Array{Float64,3}(undef,NpointsSpace,NField,timePointsUsedForOneStep-1)
    Rcoord = CartesianIndices(pointsFieldSpace)
    for iT in 1:timePointsUsedForOneStep-1
        for iField in 1:NField
            for j in Rcoord
                linearJ = LinearIndices(Rcoord)[j]
                symbKnownField[linearJ,iField,iT]=fieldLHS[iField,iT][j]
                knownField[linearJ,iField,iT] = initialCondition
            end
        end
    end

    symbUnknownField = Array{Num,2}(undef,NpointsSpace,NField)
    unknownField = Array{Float64,2}(undef,NpointsSpace,NField)
    for iField in 1:NField
        for j in Rcoord
            linearJ = LinearIndices(Rcoord)[j]
            symbUnknownField[linearJ,iField] = fieldLHS[iField,end][j]
            unknownField[linearJ,iField] = initialCondition
        end
    end

    symbKnownForce = nothing
    knownForce =nothing
    if champsLimité === nothing
        symbKnownForce=Array{Num,3}(undef,NpointsSpace,NField,timePointsUsedForOneStep)
        knownForce=Array{Float64,3}(undef,NpointsSpace,NField,timePointsUsedForOneStep)
        for iT in 1:timePointsUsedForOneStep
            for iField in 1:NField
                for j in Rcoord
                    linearJ = LinearIndices(Rcoord)[j]
                    symbKnownForce[linearJ,iField,iT]=fieldRHS[iField,iT][j]
                    knownForce[linearJ,iField,iT] = initialCondition
                end
            end
        end
    else
        NpointsLimités = length(champsLimité[1,1])
        symbKnownForce=Array{Num,3}(undef,NpointsLimités,NField,timePointsUsedForOneStep)
        knownForce=Array{Float64,3}(undef,NpointsLimités,NField,timePointsUsedForOneStep)
        for iT in 1:timePointsUsedForOneStep
            for iField in 1:NField
                for j in eachindex(NpointsLimités)
                    symbKnownForce[j,iField,iT]=champsLimité[iField,iT][j]
                    knownForce[j,iField,iT] = initialCondition
                end
            end
        end
    end
   
    #endregion

    #region numerical operators

    f= buildNumericalFunctions(costfunctions,symbUnknownField,symbKnownField,symbKnownForce)

    #end


    #region sparse matrix colouring 

    

    J,cache=sparseColouring(f,unknownField,knownField,knownForce)
 
    #endregion

    #region preparation for time marching

    if !isfile(sequentialFileName)


        if sourceType === "Ricker" && timePointsUsedForOneStep !== 1
        # source time function will be shifted with timePointsUsedForOneStep - 1

            prepend!(sourceTime,zeros(timePointsUsedForOneStep))
        else sourceType === "Explicit" && timePointsUsedForOneStep !==1
            # prepend zeros for kindness but maybe it's aweful to code it
        end


        unknownField .= initialCondition

        
        fieldFile = jldopen(sequentialFileName, "w") # file open

        #endregion

        #region make figure
        # reshape
        newField = reshape(unknownField,NField,pointsFieldSpace...)
        # Construct the full slice dynamically
        ndim = ndims(newField)
        slice = (1, ntuple(i -> Colon(), ndim - 1)...)  # (1, :, :, ..., :)
        # Apply the slice
        slice_for_field = newField[1,:,:]

        if videoMode
            fig = Figure()
            ax = Axis(fig[1, 1])
            hm = heatmap!(ax, Float32.(slice_for_field), colormap = :deep, colorrange = (-1e-5, 1e-5))
            display(fig)
        end
        #endregion


        #region time marching scheme
   
        for it in itVec
            if sourceType === "Ricker"
                knownForce[1:end,1:end,1:timePointsUsedForOneStep] .= sourceTime[it:it+timePointsUsedForOneStep-1]
            else sourceType === "Explicit"
                knownForce[1:end,1:end,1:timePointsUsedForOneStep] .= sourceFull[1:end,1:end,it:it+timePointsUsedForOneStep-1]
            end
            #field to be shifted from the past
            
            # this is not true!!! Just for debugging!
            #knownForce[1:timePointsUsedForOneStep] .= 1.0

            timeStepOptimisation!(f,unknownField,knownField,knownForce,J,cache,NpointsSpace,NField)
            @show maximum(unknownField)
            # update
            if length(itVec) > 1
                knownField[:,:,1:end-1] .= knownField[:,:,2:end]
                knownField[:,:,end] .= unknownField[:,:]
            end
   

            # reshape
            newField = reshape(unknownField,NField,pointsFieldSpace...)

            fieldFile["timestep_$it"] = newField

            if videoMode

                slice_for_field = Float32.(newField[1, :, :])

                # properly update the heatmap plot
                hm[1] = slice_for_field
            
                @show it

            end
        end
        close(fieldFile)

    end

    #endregion


    #region making the file compact

    # Reconstruct array from separate datasets

    file = jldopen(sequentialFileName, "r")
    a = [file["timestep_$it"] for it in itVec]
    close(file)

    

    if videoMode
    fig = Figure()
    ax = Axis(fig[1, 1])
        
    hm=heatmap!(ax,Float32.(a[1])[1,:,:],colormap = :plasma, colorrange = (-1e-6, 1e-6))
    Colorbar(fig[1, 2], hm)
    display(fig)
    end
    #endregion

    #region compact file
    # Convert to single 3D array: a[N, L, K]
    #a = permutedims(reshape(reduce(hcat, a), NField,pointsField...,Nt), (3, 1, 2))
    
    a=reduce(hcat,a)

    # Now save as compact single dataset
    @save compactFileName a compress=true

    #endregion

    # this is only for benchmark 1d
    return a
    # only only for benchmark 1d

end