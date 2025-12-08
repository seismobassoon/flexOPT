using YAML,StructTypes,CSV,DataFrames,Unitful

@doc raw"""
    planet1Dconfig

    This structure contains the parameters used in the planet1D code. 
    The parameters are read from a file called planet1Dconfig.txt. 
    If the file is not found, the default values are used. 
    The file must be in the same directory as the planet1D executable.

    planet1Dconfig.re: relative error (see GT95 eq. 6.2: you can control the quality of synthetics)
    planet1Dconfig.ratc: ampratio using in grid cut-off (1.d-10 is recommended)
    planet1Dconfig.ratl: ampratio using in l-cutoff
    planet1Dconfig.omegaiTlen: wrap-around attenuation for omegaiTlen (usually 1.d-2 is used)
    planet1Dconfig.maxlmax: maximum lmax

    Example of file:
    re = 1.e-2
    ratc = 1.e-10
    ratl = 1.e-5
    omegaiTlen = 1.e-2
    maxlmax = 80000

    initial values in case the variables are not found
"""
input=Input()
metaInfo=MetaInfo()
dsm1Dconfig=planet1Dconfig()

try 
    data = YAML.load_file(dirname(@__FILE__) *"/../../dataInput/config.yaml"; dicttype=Dict{Symbol, Any})
    test_dict = data[:classicDSM][1]
    StructTypes.StructType(::Type{MetaInfo}) = StructTypes.Mutable()
    global metaInfo = StructTypes.constructfrom(MetaInfo, test_dict)
catch
    @error("planet1D configuration file not found or not readable", dirname(@__FILE__) *"/../../dataInput/config.yaml")
end

try 
    data = YAML.load_file(dirname(@__FILE__) *"/../"* metaInfo.planet1DconfigFile; dicttype=Dict{Symbol, Any})
    test_dict = data[:global][1]
    StructTypes.StructType(::Type{planet1Dconfig}) = StructTypes.Mutable()
    global dsm1Dconfig = StructTypes.constructfrom(planet1Dconfig, test_dict)
catch
    @warn("planet1D configuration file not found")
    # see mainStructures.jl for the default values
    dsm1Dconfig.re = 1.e-2
    dsm1Dconfig.ratc = 1.e-10
    dsm1Dconfig.ratl = 1.e-5
    dsm1Dconfig.omegaiTlen = 1.e-2
    dsm1Dconfig.maxlmax = 80000
    dsm1Dconfig.傾き許容度 = 2.0
    dsm1Dconfig.eps = 1.5e-3
    global dsm1Dconfig.modelFolder = "../../dataInput/models"
end



try 
    paramFile= dirname(@__FILE__) *"/../"*Main.ParamFile
    data = CSV.read(paramFile,DataFrame; header=false, comment="#")
    dic=Dict(Pair.(data.Column1, data.Column2))
    input.modelFile=strip(dic["modelFile"])
    if haskey(dic,"averagedPlanetRadius")
        input.averagedPlanetRadius=uparse(strip(dic["averagedPlanetRadius"]))
    else
        input.averagedPlanetRadius=0.e0
    end

   


    if haskey(dic, "timeWindowMinimum")
        timeWindowMinimum=uparse(strip(dic["timeWindowMinimum"]))
    else
        timeWindowMinimum=0.e0
    end
    
    
    if haskey(dic, "timeWindow")
        input.timeWindow=uparse(strip(dic["timeWindow"]))
    else
        input.timeWindow=0.e0
    end


    if timeWindowMinimum == 0.e0 && input.timeWindow == 0.e0
        input.timeWindow=3276.8e0
    else
        if input.timeWindow==0.e0
            tw=1.e-1
            while tw < timeWindowMinimum
                tw*=2.e0
            end
            input.timeWindow=tw
        end
    end
   
    input.ωᵢ=-log(dsm1Dconfig.omegaiTlen)/input.timeWindow

    if haskey(dic, "maxFrequencyMin")
        maxFrequencyMin=uparse(strip(dic["maxFrequencyMin"]))
    else
        @error("maxFrequencyMin not found in ParamFile $paramFile")
    end

    if haskey(dic, "minFrequencyMax")
        minFrequencyMax=uparse(strip(dic["minFrequencyMax"]))
    else
        minFrequencyMax=0.e0
    end

    # here we compute the angular frequency array

    dω=2.e0*π/input.timeWindow
    iEnd=2^trunc(Int64,log(maxFrequencyMin*input.timeWindow)/log(2.e0))
    if minFrequencyMax>0.e0
        2^trunc(Int64,log(minFrequencyMax*input.timeWindow)/log(2.e0))
    else
        iStart=1
    end
    input.ωᵣ=range(dω*iStart,stop=dω*iEnd,length=iEnd-iStart+1)
    

    if haskey(dic, "maxFrequencyMin")
        maxFrequencyMin=uparse(strip(dic["maxFrequencyMin"]))
    else
        @error("maxFrequencyMin not found in ParamFile $paramFile")
    end

    if haskey(dic, "minFrequencyMax")
        minFrequencyMax=uparse(strip(dic["minFrequencyMax"]))
    else
        minFrequencyMax=0.e0
    end

    # here we compute the angular frequency array



    if haskey(dic,"GUIoption")
        
        if strip(dic["GUIoption"])=="on"
            input.GUIoption=true
        else
            input.GUIoption=false
        end
    else
        input.GUIoption=false
    end
catch
    @error("ParamFile file not found or not readable", dirname(@__FILE__) *"/../"*Main.ParamFile)
end



