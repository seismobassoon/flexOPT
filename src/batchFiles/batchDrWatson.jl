using DrWatson, JLD2
using SHA
using FilePathsBase: mkpath

function hash_parameters(params)
    # Serialize as string and hash
    filtered = Dict(k => v for (k,v) in params)
    str = repr(filtered)
    return bytes2hex(sha1(str))[1:8]  # Short hash
end

function myProduceOrLoad(functionName,paramDict,directoryName::String)
    output=myProduceOrLoad(functionName,paramDict,directoryName,directoryName)
    return output
end

function myProduceOrLoad(functionName,paramDict,directoryName::String,prefixName::String)

    # this is myFunction strategy for DrWatson

    hash_id = hash_parameters(paramDict)
    

    newDict = Dict{String,Any}(paramDict)
    newDict["hash_id"] = hash_id
    
    output, _ = produce_or_load(functionName,newDict,datadir(directoryName);filename = config -> "$(prefixName)_$(newDict["hash_id"])")
    
    return output

end


# lazy version of myProduceOrLoad

function lazyProduceOrLoad(x, f, args...; folder="./tmp", kwargs...)
    # Ensure the folder exists
    mkpath(folder)
    
    # Build filename
    xString = string(x)
    filename = joinpath(folder, xString * ".jld2")
    
    if isfile(filename)
        # Load saved result
        println("Loading from ", filename)
        @load filename output
        return output
    else
        # Compute and save
        println("Computing ", xString)
        output = f(args...;kwargs)  # call the function with parameters
        @save filename output
        return output
    end
end
