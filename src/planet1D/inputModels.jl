using CSV, DataFrames,BenchmarkTools
@doc raw"""

    test1DModelType: test the model type and return the model type

    readConvertedModel: cll DSM1DPSVmodel constructor to convert the model to DSM1D format
    Note that the max radius of the solid model should be treated as the planetary radius for DSM/MINEOS while the max depth of the TauP model should be treated as the planetary radius

    read1DModel(modelFile::String):This function will read the input model from a file in DSM native format or MINEOS format 
    TauP format can be read and will interpreted with some a priori information

    readTopographyModel(modelFile::String): This function will read the topography model in XXX format
    """

function readConvertedModel(rawArray::Vector{Any})
    
    # first we need to know the model type
    modelType=test1DModelType(rawArray)

    if DSM1D.input.averagedPlanetRadius === nothing
        DSM1D.input.averagedPlanetRadius=0.e0
    end
    if modelType=="DSM"|| modelType == "DSM-SH" || modelType == "DSM-Q-Complete"
        array=convert(Array{Float64,1}, rawArray[2:end])
        tmpDSM1Dmodel= DSM1DPSVmodel(Int64(rawArray[1]),modelType,DSM1D.input.averagedPlanetRadius,array)
        #println(tmpDSM1Dmodel)
    elseif modelType=="MINEOS" || modelType=="TauP"
        array=convert(Array{Float64,1}, rawArray[1:end])
        tmpDSM1Dmodel=DSM1DPSVmodel(modelType,array,DSM1D.dsm1Dconfig.傾き許容度,DSM1D.dsm1Dconfig.eps)
    else
        @error("Model type $modelType not recognized")
        # construct first DSM1D PSV modelthen we find the way to convert it to DSM1D format 
    end
    return tmpDSM1Dmodel
end 

function test1DModelType(array::Vector{Any})

    ModelType="unknown"
    nzoneMaybe=Int64(array[1])

    if rem(length(array),4) == 0
        ModelType = "TauP"
    else
        if length(array)==nzoneMaybe*28+1
            ModelType = "DSM"
        elseif length(array)==nzoneMaybe*15+1
            ModelType = "DSM-SH"
        elseif length(array)==nzoneMaybe*42+1
            ModelType = "DSM-Q-Complete"
        else
            nzoneMaybe=Int64(array[4])
            if length(array)==nzoneMaybe*9+7
                ModelType = "MINEOS"
            end
        end
    end
    return ModelType
end

function read1DModel(modelFile::String)

    # this function will read the model and grid it
    try        
        # honestly the dinosaur way to read a file but it works
        a=readlines(modelFile)
        b=[]
        for i in eachindex(a)
            
            m = match(r"^\s*[+-]?\d+", a[i]) # match the first digit

            if m === nothing
                # if the line starts from non numeric, then it is a comment
            else
                pat=r"[+-]?\d*\.{0,1}[Ee]?\d+"
                numbers = map(eachmatch(pat, a[i])) do mm
                    parse(Float64, mm.match)
                end
                append!(b,numbers)
            end 
        end
        myDSMmodel=readConvertedModel(b)
        return myDSMmodel
    catch
         @error("Model file $modelFile not found or not in the right format")
    end
end

function readTopographyModel(modelFile::String,type="free surface")
    # type can be "free surface" or "solid on liquid"
    try        
        a=readdlm(IOBuffer(replace(read(modelFile), UInt8('\t') => UInt8(' '))), ' '; header=false, comments=true, comment_char='#')
        #data=CSV.File(modelFile; header=false, delim = ' ', comments=true, comment_char='#')
        #println(data)
        a=vec(a)
        a=filter(!=(""),a)
        # filter out no digit characters
        #a=filter(x->all(isnumber.(x)),a)
        # filter all the floats and integers
        #a=filter(x->typeof(x)!=String,a)
      
    catch
        @error("Model file $modelFile not found")
    end
    
    return
end 



# Below are the execution commands





# @btime readTopographyModel(TopographyFile1,"solid on liquid")
# 3D model (perturbation file)