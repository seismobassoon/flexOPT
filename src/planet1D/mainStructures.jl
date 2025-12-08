using Printf


# VERY SERIOUS PROBLEM!! :

# I don't remember why but all the polynomial coefficients are stored in the weird order [izone, iOrder]
# I am reluctant to modify this for the moment but it should be done soon (10/07/2025)
# NF

mutable struct Input
    modelFile::String
    averagedPlanetRadius::Float64
    timeWindow::Float64
    ωᵢ::Float64
    ωᵣ::Array{Float64,1}
    GUIoption::Bool
    Input() = new()
end

mutable struct MetaInfo 
    name::String
    info::String
    planet1DconfigFile::String
    MetaInfo() = new()
end

mutable struct planet1Dconfig 
    re::Float64 # (re) relative error (see GT95 eq. 6.2: you can control the quality of synthetics)
    ratc::Float64 # ampratio using in grid cut-off (1.d-10 is recommended)
    ratl::Float64 # ampratio using in l-cutoff
    omegaiTlen::Float64 # wrap-around attenuation for omegaiTlen (usually 1.d-2 is used)
    # input.omegai=-log(dsm1Dconfig.omegaiTlen)/input.timeWindow
    maxlmax::Int64 # maximum lmax
    modelFolder::String
    傾き許容度::Float64 # tolerance for the slope change to be considered for a fusion
    eps::Float64 # tolerance for the error in interpolation
    planet1Dconfig() = new()
end


mutable struct planet1DPSVmodel 
    nzone::Int64
    averagedPlanetRadiusInKilometer::Float64
    averagedPlanetCMBInKilometer::Float64
    solid_or_fluid::Array{String,1} # U: unknown, S: solid, F: fluid
    bottomRadius::Array{Float64,1}
    topRadius::Array{Float64,1}
    normalisedBottomRadius::Array{Float64,1} # normalised to averagedPlanetRadiusInKilometer
    normalisedTopRadius::Array{Float64,1}
    normalisedCMB::Float64
    C_ρ     ::Array{Float64,2}
    C_Vpv   ::Array{Float64,2}
    C_Vph   ::Array{Float64,2}
    C_Vsv   ::Array{Float64,2}
    C_Vsh   ::Array{Float64,2}
    C_Qκ    ::Array{Float64,2}
    C_Qμ    ::Array{Float64,2}
    C_QκPower    ::Array{Float64,2}
    C_QμPower    ::Array{Float64,2}
    C_η     ::Array{Float64,2}
    
    function planet1DPSVmodel()
        # This constructs an empty struct with the right dimensions
        return new(1,0.e0,["U"],
        [0.e0],[0.e0], [0.e0],[0.e0],
        zeros(1,4),zeros(1,4),zeros(1,4),
        zeros(1,4),zeros(1,4),zeros(1,4),
        zeros(1,4),zeros(1,4),zeros(1,4),
        zeros(1,4))
    end
    function planet1DPSVmodel(nzone::Int64,averagedPlanetRadius::Float64)
        # This constructs an empty struct with the right dimensions
        averagedPlanetCMBInKilometer::Float64 = 0.0
        normalisedCMB::Float64 = 0.0
        return new(nzone,averagedPlanetRadius,averagedPlanetCMBInKilometer,Array{String,1}(undef,nzone),
        Array{Float64,1}(undef,nzone),Array{Float64,1}(undef,nzone),
        Array{Float64,1}(undef,nzone),Array{Float64,1}(undef,nzone),
        normalisedCMB,
        Array{Float64,2}(undef,nzone,4),Array{Float64,2}(undef,nzone,4),Array{Float64,2}(undef,nzone,4),
        Array{Float64,2}(undef,nzone,4),Array{Float64,2}(undef,nzone,4),Array{Float64,2}(undef,nzone,4),
        Array{Float64,2}(undef,nzone,4),Array{Float64,2}(undef,nzone,4),Array{Float64,2}(undef,nzone,4),
        Array{Float64,2}(undef,nzone,4))
    end


    

    function planet1DPSVmodel(modelType::String,rawArray::Array{Float64,1},傾き許容度::Float64,eps::Float64)
        # This subroutine constructs a planet1D model from a MINEOS model

        if modelType == "MINEOS" || modelType== "TauP" # for MINEOS
            if modelType == "MINEOS"
                tmpRadius = rawArray[8:9:end]*1.e-3
                tmpDensity = rawArray[9:9:end]*1.e-3
                tmpVpv = rawArray[10:9:end]*1.e-3
                tmpVsv = rawArray[11:9:end]*1.e-3
                tmpQκ = rawArray[12:9:end]
                tmpQμ = rawArray[13:9:end]
                tmpVph = rawArray[14:9:end]*1.e-3
                tmpVsh = rawArray[15:9:end]*1.e-3
                tmpη = rawArray[16:9:end]
                tmpaveragedPlanetRadius=0.e0
            elseif modelType == "TauP"
                nLayers=convert(Int64,length(rawArray)/4)
                tmpaveragedPlanetRadius=rawArray[end-3]
                tmpRadius=tmpaveragedPlanetRadius.-rawArray[end-3:-4:1]
                tmpVp=rawArray[end-2:-4:2]
                tmpVs=rawArray[end-1:-4:3]
                tmpDensity=rawArray[end:-4:4]
                tmpVpv=tmpVp
                tmpVsv=tmpVs
                tmpQκ=zeros(nLayers)
                tmpQμ=zeros(nLayers)
                tmpVph=tmpVp
                tmpVsh=tmpVs
                tmpη=ones(nLayers)
                tmpaveragedPlanetRadius=0.e0
            end
            
            for i in eachindex(tmpRadius)    
                if tmpVsv[i] != 0.0
                    tmpaveragedPlanetRadius=tmpRadius[i]
                end
            end
            normalisedRadius = tmpRadius./tmpaveragedPlanetRadius
            # we don't change the number of layers for each parameter, we will fusion at the end
            nLayers=0
            normalisedtmpBottomRadius=Float64[]
            normalisedtmpTopRadius=Float64[]
            tmpBottomRadius=Float64[]
            tmpTopRadius=Float64[]
            tmpCoefs=(Density=[], Vpv=[], Vsv=[], Qκ=[], Qμ=[], Vph=[], Vsh=[], η=[])

            for i in 1:length(normalisedRadius)-1
                if tmpRadius[i] != tmpRadius[i+1]
                    nLayers+=1
                    # Density
                    a=(tmpDensity[i+1]-tmpDensity[i])/(normalisedRadius[i+1]-normalisedRadius[i])
                    b=tmpDensity[i]-a*normalisedRadius[i]
                    push!(tmpCoefs.Density,[b,a])
                    push!(normalisedtmpBottomRadius,normalisedRadius[i])
                    push!(normalisedtmpTopRadius,normalisedRadius[i+1])
                    push!(tmpBottomRadius,tmpRadius[i])
                    push!(tmpTopRadius,tmpRadius[i+1])
                    # Vpv
                    a=(tmpVpv[i+1]-tmpVpv[i])/(normalisedRadius[i+1]-normalisedRadius[i])
                    b=tmpVpv[i]-a*normalisedRadius[i]
                    push!(tmpCoefs.Vpv,[b,a])
                    # Vsv
                    a=(tmpVsv[i+1]-tmpVsv[i])/(normalisedRadius[i+1]-normalisedRadius[i])
                    b=tmpVsv[i]-a*normalisedRadius[i]
                    push!(tmpCoefs.Vsv,[b,a])
                    # Qκ
                    a=(tmpQκ[i+1]-tmpQκ[i])/(normalisedRadius[i+1]-normalisedRadius[i])
                    b=tmpQκ[i]-a*normalisedRadius[i]
                    push!(tmpCoefs.Qκ,[b,a])
                    # Qμ
                    a=(tmpQμ[i+1]-tmpQμ[i])/(normalisedRadius[i+1]-normalisedRadius[i])
                    b=tmpQμ[i]-a*normalisedRadius[i]
                    push!(tmpCoefs.Qμ,[b,a])
                    # Vph
                    a=(tmpVph[i+1]-tmpVph[i])/(normalisedRadius[i+1]-normalisedRadius[i])
                    b=tmpVph[i]-a*normalisedRadius[i]
                    push!(tmpCoefs.Vph,[b,a])
                    # Vsh
                    a=(tmpVsh[i+1]-tmpVsh[i])/(normalisedRadius[i+1]-normalisedRadius[i])
                    b=tmpVsh[i]-a*normalisedRadius[i]
                    push!(tmpCoefs.Vsh,[b,a])
                    # η
                    a=(tmpη[i+1]-tmpη[i])/(normalisedRadius[i+1]-normalisedRadius[i])
                    b=tmpη[i]-a*normalisedRadius[i]
                    push!(tmpCoefs.η,[b,a])
                end
            end
 
            Coefs=planet1DPSVmodel(nLayers,tmpaveragedPlanetRadius)

            Coefs.C_ρ =modifyLinearCoefs(nLayers, tmpCoefs.Density, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
            Coefs.C_Vpv=modifyLinearCoefs(nLayers, tmpCoefs.Vpv, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
            Coefs.C_Vsv=modifyLinearCoefs(nLayers, tmpCoefs.Vsv, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
            Coefs.C_Qκ=modifyLinearCoefs(nLayers, tmpCoefs.Qκ, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
            Coefs.C_Qμ=modifyLinearCoefs(nLayers, tmpCoefs.Qμ, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
            Coefs.C_Vph=modifyLinearCoefs(nLayers, tmpCoefs.Vph, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
            Coefs.C_Vsh=modifyLinearCoefs(nLayers, tmpCoefs.Vsh, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
            Coefs.C_η=modifyLinearCoefs(nLayers, tmpCoefs.η, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
            Coefs.C_QκPower=zeros(nLayers,4)
            Coefs.C_QμPower=zeros(nLayers,4)
            Coefs.topRadius=tmpTopRadius
            Coefs.bottomRadius=tmpBottomRadius
            Coefs.normalisedTopRadius=normalisedtmpTopRadius
            Coefs.normalisedBottomRadius=normalisedtmpBottomRadius

            # We then find out the doublons

            Coefs=reduceLayersDSM(Coefs)
            #arrayRadius, arrayParams = compute1DseismicParamtersFromPolynomialCoefficients(Coefs, 10)
            for i in 1:Coefs.nzone
                if Coefs.C_Vsv[i,:] == zeros(4)
                    Coefs.solid_or_fluid[i] = "F"
                else
                    Coefs.solid_or_fluid[i] = "S"
                end
            end
        else 
            @error "modelType $modelType is not proper to use this MINEOS-DSM constructor"
        end

        return Coefs
    end

    function planet1DPSVmodel(nzone::Int64,modelType::String,averagedPlanetRadius::Float64,array::Array{Float64,1})
        # This is for DSM classic / DSM-SH / DSM-Q-Complete 
        
        Csolid_or_fluid=Array{String,1}(undef,nzone)
        if modelType == "DSM" # DSM classic
            CbottomRadius = array[1:28:end]
            CtopRadius = array[2:28:end]
            CC_ρ=reduce(hcat,[array[3:28:end],array[4:28:end],array[5:28:end],array[6:28:end]])
            CC_Vpv=reduce(hcat,[array[7:28:end],array[8:28:end],array[9:28:end],array[10:28:end]])
            CC_Vph=reduce(hcat,[array[11:28:end],array[12:28:end],array[13:28:end],array[14:28:end]])
            CC_Vsv=reduce(hcat,[array[15:28:end],array[16:28:end],array[17:28:end],array[18:28:end]])
            CC_Vsh=reduce(hcat,[array[19:28:end],array[20:28:end],array[21:28:end],array[22:28:end]])
            CC_η=reduce(hcat,[array[23:28:end],array[24:28:end],array[25:28:end],array[26:28:end]])
            CC_Qμ=reduce(hcat,[array[27:28:end],zeros(nzone),zeros(nzone),zeros(nzone)])
            CC_Qκ=reduce(hcat,[array[28:28:end],zeros(nzone),zeros(nzone),zeros(nzone)])
            CC_QμPower=zeros(nzone,4)
            CC_QκPower=zeros(nzone,4)
        elseif modelType == "DSM-SH" # DSM-SH
            CbottomRadius = array[1:15:end]
            CtopRadius = array[2:15:end]
            CC_ρ=reduce(hcat,[array[3:15:end],array[4:15:end],array[5:15:end],array[6:15:end]])
            CC_Vpv=zeros(nzone,4)
            CC_Vph=zeros(nzone,4)
            CC_Vsv=reduce(hcat,[array[7:15:end],array[8:15:end],array[9:15:end],array[10:15:end]])
            CC_Vsh=reduce(hcat,[array[11:15:end],array[12:15:end],array[13:15:end],array[14:15:end]])
            CC_η=reduce(hcat,[ones(nzone),zeros(nzone),zeros(nzone),zeros(nzone)])
            CC_Qμ=reduce(hcat,[array[15:15:end],zeros(nzone),zeros(nzone),zeros(nzone)])
            CC_Qκ=zeros(nzone,4)    
            CC_QμPower=zeros(nzone,4)
            CC_QκPower=zeros(nzone,4)
        elseif modelType == "DSM-Q-Complete" # DSM-Q-Complete
            CbottomRadius = array[1:42:end]
            CtopRadius = array[2:42:end]
            CC_ρ=reduce(hcat,[array[3:42:end],array[4:42:end],array[5:42:end],array[6:42:end]])
            CC_Vpv=reduce(hcat,[array[7:42:end],array[8:42:end],array[9:42:end],array[10:42:end]])
            CC_Vph=reduce(hcat,[array[11:42:end],array[12:42:end],array[13:42:end],array[14:42:end]])
            CC_Vsv=reduce(hcat,[array[15:42:end],array[16:42:end],array[17:42:end],array[18:42:end]])
            CC_Vsh=reduce(hcat,[array[19:42:end],array[20:42:end],array[21:42:end],array[22:42:end]])
            CC_η=reduce(hcat,[array[23:42:end],array[24:42:end],array[25:42:end],array[26:42:end]])
            CC_Qμ=reduce(hcat,[array[27:42:end],array[28:42:end],array[29:42:end],array[30:42:end]])
            CC_Qκ=reduce(hcat,[array[31:42:end],array[32:42:end],array[33:42:end],array[34:42:end]])
            CC_QμPower=reduce(hcat,[array[35:42:end],array[36:42:end],array[37:42:end],array[38:42:end]])
            CC_QκPower=reduce(hcat,[array[39:42:end],array[40:42:end],array[41:42:end],array[42:42:end]])
        else
            @error "modelType $modelType is not proper to use this constructor"
        end

        tmpaveragedPlanetRadius=0.e0
        tmpaveragedPlanetCMB=0.e0
        for i in 1:nzone
            if CC_Vsv[i,:] == zeros(4)
                Csolid_or_fluid[i] = "F"
                if i < nzone && CC_Vsv[i+1,:] !== zeros(4)
                    tmpaveragedPlanetCMB=CtopRadius[i]
                end
            else
                Csolid_or_fluid[i] = "S"
                tmpaveragedPlanetRadius=CtopRadius[i]
            end
        end
     
        if averagedPlanetRadius == 0.e0
            averagedPlanetRadius=tmpaveragedPlanetRadius
        end


        CnormalisedBottomRadius = CbottomRadius/averagedPlanetRadius
        CnormalisedTopRadius = CtopRadius/averagedPlanetRadius
        CnormalisedCMB = tmpaveragedPlanetCMB/averagedPlanetRadius
     

        return new(nzone,averagedPlanetRadius,tmpaveragedPlanetCMB,Csolid_or_fluid,
        CbottomRadius,CtopRadius,CnormalisedBottomRadius,CnormalisedTopRadius,CnormalisedCMB,
        CC_ρ,CC_Vpv,CC_Vph,CC_Vsv,CC_Vsh,CC_Qκ,CC_Qμ,CC_QκPower,CC_QμPower,CC_η)
    end
end


mutable struct VerticalGridStructure
    # This should contain the vertical grids for the DSM computation (SH and P-SV)
    # This structure can differ for different localωMax and locallCritical
    # The choice of locallCritical (max for big l and min for small l) is described in Kawai et al. 2006
    # However, neither Kawai version (EPS, UTokyo) nor Takeuchi version (ERI, UTokyo) uses this criterion.
    # In fact, Kawai's calgrid can potentially handle this but they keep locallCritical=0
    # On the top of this, there is a rapid decay for big l: thus this vertical griding should be revisited
    numberOfGrids::Int64
    numberOfSolidGrids::Int64
    RadiusInKm::Array{Float64,1}
    RadiusSolidInKm::Array{Float64,1} 
    normalisedRadius::Array{Float64,1}
    normalisedSolidRadius::Array{Float64,1}
    solid_or_fluid::Array{String,1}


    function VerticalGridStructure(Coefs::planet1DPSVmodel,relativeError::Float64,localωMax::Float64,locallCritical::Int64=0)
        # We scan the model by using the equation 6.2 of Geller & Takeuchi 1995
        tmpRadius=Float64[]
        tmpNormalisedRadius=Float64[]
        tmpSolidRadius=Float64[]
        tmpSolidNormalisedRadius=Float64[]
        tmpSolid_or_fluid=String[]
        numberOfGrids=0


        #@show Coefs.normalisedTopRadius
        for i in 1:Coefs.nzone

            #@show Coefs.normalisedBottomRadius[i], Coefs.normalisedTopRadius[i],Coefs.C_Vpv[i,:]
            if Coefs.solid_or_fluid[i] == "S"
                vmin=min(getParameterForOnePoint(Coefs.C_Vsv[i,:],Coefs.normalisedBottomRadius[i]),getParameterForOnePoint(Coefs.C_Vsv[i,:],Coefs.normalisedTopRadius[i])) 
            else
                vmin=min(getParameterForOnePoint(Coefs.C_Vpv[i,:],Coefs.normalisedBottomRadius[i]),getParameterForOnePoint(Coefs.C_Vpv[i,:],Coefs.normalisedTopRadius[i]))
            end


            rmax=Coefs.topRadius[i]
            kₓ=(locallCritical+5.e-1)/rmax
            Δz = sqrt(abs(relativeError/3.3e0*4.e0*π^2/(localωMax^2/vmin^2-kₓ^2)))
            zWidth=Coefs.topRadius[i]-Coefs.bottomRadius[i]
 
            localNumberOfGrids=max(trunc(Int64,zWidth/Δz)+1,5)
            numberOfGrids+=localNumberOfGrids
            Δz= (Coefs.topRadius[i]-Coefs.bottomRadius[i])/convert(Float64,localNumberOfGrids-1)
            normalisedΔz=Δz/Coefs.averagedPlanetRadiusInKilometer
            for j in 1:localNumberOfGrids
                push!(tmpRadius,Coefs.bottomRadius[i]+Δz*convert(Float64,j-1))
                push!(tmpNormalisedRadius,Coefs.normalisedBottomRadius[i]+normalisedΔz*convert(Float64,j-1))
                push!(tmpSolid_or_fluid,Coefs.solid_or_fluid[i])
            end
        end
   
        # We then find only the last continuous solid layers
    
        a=findall(x->x=="S",tmpSolid_or_fluid)   
        n=length(a)
        iStartSolid=n
        for i in n:-1:2
            if(a[i]-a[i-1]>1)
                iStartSolid=i
                break
            else
                iStartSolid=i-1
            end
        end
        numberOfSolidGrids=n-iStartSolid+1
        for i in a[iStartSolid:n]
            push!(tmpSolidRadius,tmpRadius[i])
            push!(tmpSolidNormalisedRadius,tmpNormalisedRadius[i])
        end

        # it's not yet done! Think about how to do this for SH waves
        return new(numberOfGrids,numberOfSolidGrids,tmpRadius,tmpSolidRadius,tmpNormalisedRadius,tmpSolidNormalisedRadius,tmpSolid_or_fluid)
    end

end


function modifyLinearCoefs(nLayers, tmpCoefs, normalisedtmpBottomRadius, normalisedtmpTopRadius, 傾き許容度,eps)
    # This is a function that takes the coefficients of the linear interpolation and tries to reduce the number of layers by fusing them up to 3rd order polynomials
    # 傾き許容度 is the tolerance for the slope change to be considered for a fusion
    # eps is the tolerance for the error in interpolation
    coefs=zeros(Float64, nLayers, 4)
    i=1 
    while i <= nLayers
    
        # we don't do anything if there is a jump
        if i+1 <= nLayers && isapprox((tmpCoefs[i+1][1]-tmpCoefs[i][1]) + (tmpCoefs[i+1][2]-tmpCoefs[i][2]) * normalisedtmpTopRadius[i] , 0.0 ; atol=eps)
           
        # we change only if the slope changes a little bit
            if i+1 <= nLayers && tmpCoefs[i][2] != 0.e0  && abs(tmpCoefs[i+1][2]-tmpCoefs[i][2])/abs(tmpCoefs[i][2]) < 傾き許容度 
                # normally this will make a second-order polynomial

                # but first we test if we can regroup one more layers above
                if i+2 <= nLayers && isapprox((tmpCoefs[i+2][1]-tmpCoefs[i+1][1]) + (tmpCoefs[i+2][2]-tmpCoefs[i+1][2]) * normalisedtmpTopRadius[i] , 0.0 ; atol=1.e-4) &&tmpCoefs[i+1][2] != 0.e0 && abs(tmpCoefs[i+2][2]-tmpCoefs[i+1][2])/abs(tmpCoefs[i+1][2]) < 傾き許容度 
                    # create an optimised Vandermonde matrix and get the coefs for a third-order polynomial
                    x₀=normalisedtmpBottomRadius[i]
                    x₁=normalisedtmpTopRadius[i]
                    x₂=normalisedtmpTopRadius[i+1]
                    x₃=normalisedtmpTopRadius[i+2]

                    A=[x₁-x₀  x₁^2-x₀^2  x₁^3-x₀^3;
                    x₂-x₁  x₂^2-x₁^2  x₂^3-x₁^3;
                    x₃-x₂  x₃^2-x₂^2  x₃^3-x₂^3]
                    #c=[c₀; c₁; c₂]
        
                    y₀=tmpCoefs[i][1]+tmpCoefs[i][2]*x₀ 
                    y₁=tmpCoefs[i][1]+tmpCoefs[i][2]*x₁
                    y₂=tmpCoefs[i+1][1]+tmpCoefs[i+1][2]*x₂
                    y₃=tmpCoefs[i+2][1]+tmpCoefs[i+2][2]*x₃

                    y=[tmpCoefs[i][2]*(x₁-x₀); tmpCoefs[i+1][2]*(x₂-x₁); tmpCoefs[i+2][2]*(x₃-x₂)]
                
                    invA=inv(A)

                    xc=invA*y
                    c₁=xc[1]
                    c₂=xc[2]
                    c₃=xc[3]
                    c₀ = y₀-c₁*x₀-c₂*x₀^2-c₃*x₀^3

                    coefs[i,:]=[c₀,c₁,c₂,c₃]
                    coefs[i+1,:]=[c₀,c₁,c₂,c₃]
                    coefs[i+2,:]=[c₀,c₁,c₂,c₃]


                    i+=1 # we skip the next layer because it is already included in the polynomial (twice since i+=1 reads at the outer loop as well)

                else
                    x₀=normalisedtmpBottomRadius[i]
                    x₁=normalisedtmpTopRadius[i]
                    x₂=normalisedtmpTopRadius[i+1]


                    A=[x₁-x₀  x₁^2-x₀^2 ;
                    x₂-x₁  x₂^2-x₁^2 ]
                    #c=[c₀; c₁; c₂]
        
                    y₀=tmpCoefs[i][1]+tmpCoefs[i][2]*x₀ 
                    y₁=tmpCoefs[i][1]+tmpCoefs[i][2]*x₁
                    y₂=tmpCoefs[i+1][1]+tmpCoefs[i+1][2]*x₂


                    y=[tmpCoefs[i][2]*(x₁-x₀); tmpCoefs[i+1][2]*(x₂-x₁)]
                
                    invA=inv(A)

                    xc=invA*y
                    c₁=xc[1]
                    c₂=xc[2]

                    c₀ = y₀-c₁*x₀-c₂*x₀^2

                    coefs[i,:]=[c₀,c₁,c₂,0.e0]
                    coefs[i+1,:]=[c₀,c₁,c₂,0.e0]

                end
                i+=1 # we skip the next layer because it is already included in the polynomial
            else
                coefs[i,:]=[tmpCoefs[i][1],tmpCoefs[i][2],0.0,0.0]
            end
        else
            coefs[i,:]=[tmpCoefs[i][1],tmpCoefs[i][2],0.0,0.0]  
        end

        i+=1 # normal incrementation
    end
    
    # The coefficients above should be 'exact' for each point 
   
    # We then scan from the bottom to the top again to see if we can further reduce the independent coefficients

    i=1
    while i <= nLayers

        if i+1 <= nLayers && coefs[i,:] != coefs[i+1,:]
            # we test the lower coefficients to see if the next layer can be fused
            x₂ = normalisedtmpTopRadius[i+1]
            realValue= coefs[i+1,1]+coefs[i+1,2]*x₂+coefs[i+1,3]*x₂^2+coefs[i+1,4]*x₂^3
            extrapolatedValueFromLowerLayer= coefs[i,1]+coefs[i,2]*x₂+coefs[i,3]*x₂^2+coefs[i,4]*x₂^3
            if abs(realValue-extrapolatedValueFromLowerLayer)/abs(realValue) < eps
                # we can fuse the two layers
                coefs[i+1,:]=coefs[i,:]
            end
        end
        i+=1
    end 
    return coefs
end 

function reduceLayersDSM(tmpCoefs::planet1DPSVmodel)
    nLayers=tmpCoefs.nzone
    j_start=[]
    j_end=[]
    i=1
    while i <= nLayers
        push!(j_start,i)
        while i+1 <= nLayers && tmpCoefs.C_ρ[i,:] == tmpCoefs.C_ρ[i+1,:] && tmpCoefs.C_Vpv[i,:] == tmpCoefs.C_Vpv[i+1,:] && tmpCoefs.C_Vph[i,:] == tmpCoefs.C_Vph[i+1,:] && tmpCoefs.C_Vsv[i,:] == tmpCoefs.C_Vsv[i+1,:] && tmpCoefs.C_Vsh[i,:] == tmpCoefs.C_Vsh[i+1,:] && tmpCoefs.C_Qκ[i,:] == tmpCoefs.C_Qκ[i+1,:] && tmpCoefs.C_Qμ[i,:] == tmpCoefs.C_Qμ[i+1,:] && tmpCoefs.C_QκPower[i,:] == tmpCoefs.C_QκPower[i+1,:] && tmpCoefs.C_QμPower[i,:] == tmpCoefs.C_QμPower[i+1,:] && tmpCoefs.C_η[i,:] == tmpCoefs.C_η[i+1,:]
            i+=1
        end
        push!(j_end,i)
        i+=1
    end
    nLayers=length(j_start)
    Coefs=planet1DPSVmodel(nLayers,tmpCoefs.averagedPlanetRadiusInKilometer)
    for i in 1:nLayers
        Coefs.topRadius[i]=tmpCoefs.topRadius[j_end[i]]
        Coefs.bottomRadius[i]=tmpCoefs.bottomRadius[j_start[i]]
        Coefs.normalisedTopRadius[i]=tmpCoefs.normalisedTopRadius[j_end[i]]
        Coefs.normalisedBottomRadius[i]=tmpCoefs.normalisedBottomRadius[j_start[i]]
        Coefs.C_ρ[i,:]=tmpCoefs.C_ρ[j_start[i],:]
        Coefs.C_Vpv[i,:]=tmpCoefs.C_Vpv[j_start[i],:]
        Coefs.C_Vph[i,:]=tmpCoefs.C_Vph[j_start[i],:]
        Coefs.C_Vsv[i,:]=tmpCoefs.C_Vsv[j_start[i],:]
        Coefs.C_Vsh[i,:]=tmpCoefs.C_Vsh[j_start[i],:]
        Coefs.C_Qκ[i,:]=tmpCoefs.C_Qκ[j_start[i],:]
        Coefs.C_Qμ[i,:]=tmpCoefs.C_Qμ[j_start[i],:]
        Coefs.C_QκPower[i,:]=tmpCoefs.C_QκPower[j_start[i],:]
        Coefs.C_QμPower[i,:]=tmpCoefs.C_QμPower[j_start[i],:]
        Coefs.C_η[i,:]=tmpCoefs.C_η[j_start[i],:]
    end
    return Coefs
end

function getParameterForOnePoint(oneCoefs::Array{Float64,1}, oneNormalisedRadius::Float64)
    # this function will compute the seismic parameters from the polynomial coefficients and the normalised radii
    return oneCoefs[1]+oneCoefs[2]*oneNormalisedRadius+oneCoefs[3]*oneNormalisedRadius^2+oneCoefs[4]*oneNormalisedRadius^3
end


function getParameterDerivativeForOnePoint(oneCoefs::Array{Float64,1}, oneNormalisedRadius::Float64)
    # this function will compute the seismic parameters from the polynomial coefficients and the normalised radii
    return oneCoefs[2]+2.e0*oneCoefs[3]*oneNormalisedRadius+oneCoefs[4]*3.e0*oneNormalisedRadius^2
end

function compute1DseismicParametersFromPolynomialCoefficients(nLayers, coefs, normalisedtmpBottomRadius, normalisedtmpTopRadius, normalisedRadius::Float64)
    # this function will compute the seismic parameters from the polynomial coefficients and the normalised radii
   
    for i in 1:nLayers
        if normalisedRadius >= normalisedtmpBottomRadius[i] && normalisedRadius <= normalisedtmpTopRadius[i]
            return getParameterForOnePoint(coefs[i,:],normalisedRadius) 
        end
    end
    @error("Radius is not in the model")
end



function compute1DseismicParamtersFromPolynomialCoefficients(Coefs::planet1DPSVmodel, numPointsPerLayer::Int64)
    # this function will compute the seismic parameters from the polynomial coefficients and the normalised radii
    nLayers=Coefs.nzone
    normalisedtmpBottomRadius=Coefs.bottomRadius./Coefs.averagedPlanetRadiusInKilometer
    normalisedtmpTopRadius=Coefs.topRadius./Coefs.averagedPlanetRadiusInKilometer

    if numPointsPerLayer < 3
        numPointsPerLayer= 3
    end
    radius=Float64[]
    parameter=(ρ=Float64[],Vpv=Float64[],Vph=Float64[],Vsv=Float64[],Vsh=Float64[],Qμ=Float64[],Qκ=Float64[],QμPower=Float64[],QκPower=Float64[],η=Float64[])
    
   
    for i in 1:nLayers
        radiusInterval = (normalisedtmpTopRadius[i]-normalisedtmpBottomRadius[i])/convert(Float64,numPointsPerLayer-1)
        for j in 1:numPointsPerLayer
            x = normalisedtmpBottomRadius[i]+radiusInterval*convert(Float64,j-1)
            push!(radius,x)
            push!(parameter.ρ,Coefs.C_ρ[i,1]+Coefs.C_ρ[i,2]*x+Coefs.C_ρ[i,3]*x^2+Coefs.C_ρ[i,4]*x^3)
            push!(parameter.Vpv,Coefs.C_Vpv[i,1]+Coefs.C_Vpv[i,2]*x+Coefs.C_Vpv[i,3]*x^2+Coefs.C_Vpv[i,4]*x^3)
            push!(parameter.Vph,Coefs.C_Vph[i,1]+Coefs.C_Vph[i,2]*x+Coefs.C_Vph[i,3]*x^2+Coefs.C_Vph[i,4]*x^3)
            push!(parameter.Vsv,Coefs.C_Vsv[i,1]+Coefs.C_Vsv[i,2]*x+Coefs.C_Vsv[i,3]*x^2+Coefs.C_Vsv[i,4]*x^3)
            push!(parameter.Vsh,Coefs.C_Vsh[i,1]+Coefs.C_Vsh[i,2]*x+Coefs.C_Vsh[i,3]*x^2+Coefs.C_Vsh[i,4]*x^3)
            push!(parameter.Qμ,Coefs.C_Qμ[i,1]+Coefs.C_Qμ[i,2]*x+Coefs.C_Qμ[i,3]*x^2+Coefs.C_Qμ[i,4]*x^3)
            push!(parameter.Qκ,Coefs.C_Qκ[i,1]+Coefs.C_Qκ[i,2]*x+Coefs.C_Qκ[i,3]*x^2+Coefs.C_Qκ[i,4]*x^3)
            push!(parameter.QμPower,Coefs.C_QμPower[i,1]+Coefs.C_QμPower[i,2]*x+Coefs.C_QμPower[i,3]*x^2+Coefs.C_QμPower[i,4]*x^3)
            push!(parameter.QκPower,Coefs.C_QκPower[i,1]+Coefs.C_QκPower[i,2]*x+Coefs.C_QκPower[i,3]*x^2+Coefs.C_QκPower[i,4]*x^3)
            push!(parameter.η,Coefs.C_η[i,1]+Coefs.C_η[i,2]*x+Coefs.C_η[i,3]*x^2+Coefs.C_η[i,4]*x^3)
        end
    end

    return radius, parameter
end

function compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(Coefs::planet1DPSVmodel, tmpRadius::Float64)
    #this is the scalar version of the lengthy function below

    nLayers=Coefs.nzone
    x = tmpRadius/Coefs.averagedPlanetRadiusInKilometer
    zoneIndex = 0

    for i in 1:nLayers
        if Coefs.bottomRadius[i] < tmpRadius <= Coefs.topRadius[i]
            zoneIndex = i
        end
    end


    if tmpRadius === 0.0 #center of the Earth
        zoneIndex = 1
    end




    i = zoneIndex
    ρ =0.0
    if zoneIndex === 0 
        ρ = 0.0
    else
        i = zoneIndex
        ρ = Coefs.C_ρ[i,1]+Coefs.C_ρ[i,2]*x+Coefs.C_ρ[i,3]*x^2+Coefs.C_ρ[i,4]*x^3
    end 
    return ρ
end




function compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(Coefs::planet1DPSVmodel, tmpRadiiInKilometer::Array{Float64,1}, belowOrAbove::String)


    # this function gives back an array of parameters for an array of radii
    # in order to accelerate, we sort the radii array from the centre to the surface

    # belowOrAbove : which to prefer when the radius coincides with the discontinuity 
    #                (no choice for the outmost boundary: below)
    # if we repeat the same radius for the discontinuity, we take below and above values

    nLayers=Coefs.nzone
    normalisedtmpBottomRadius=Coefs.bottomRadius./Coefs.averagedPlanetRadiusInKilometer
    normalisedtmpTopRadius=Coefs.topRadius./Coefs.averagedPlanetRadiusInKilometer

    #radius=Float64[]

    # here we sort the radii

    tmpRadiiInKilometer = sort(tmpRadiiInKilometer)
    numberRadii = length(tmpRadiiInKilometer)

    tmpNormalisedRadii = tmpRadiiInKilometer./Coefs.averagedPlanetRadiusInKilometer


    #parameter=(ρ=Float64[],Vpv=Float64[],Vph=Float64[],Vsv=Float64[],Vsh=Float64[],Qμ=Float64[],Qκ=Float64[],QμPower=Float64[],QκPower=Float64[],η=Float64[])


    parameter=(ρ=zeros(Float64,numberRadii),Vpv=zeros(Float64,numberRadii),Vph=zeros(Float64,numberRadii),Vsv=zeros(Float64,numberRadii),Vsh=zeros(Float64,numberRadii),Qμ=zeros(Float64,numberRadii),Qκ=zeros(Float64,numberRadii),QμPower=zeros(Float64,numberRadii),QκPower=zeros(Float64,numberRadii),η=zeros(Float64,numberRadii))

    N深さインデックス=size(tmpNormalisedRadii)[1]
 
    深さインデックス=1
    i = 1

    twoPointsAtTheSameRadius = false
    twoPointsAtTheSameRadiusOnceDone  = false

    while i <= nLayers && 深さインデックス <=  N深さインデックス
        distanceFromBottom = tmpNormalisedRadii[深さインデックス]-normalisedtmpBottomRadius[i]
        distanceFromTop = tmpNormalisedRadii[深さインデックス]-normalisedtmpTopRadius[i]
        negativeWhenItsIn=distanceFromBottom*distanceFromTop

        if 深さインデックス < N深さインデックス
            if tmpNormalisedRadii[深さインデックス] == tmpNormalisedRadii[深さインデックス+1]
                twoPointsAtTheSameRadius = true
                #@show tmpRadiiInKilometer[深さインデックス]
            else
                twoPointsAtTheSameRadius = false
            end
        end


        if negativeWhenItsIn < 0 || (belowOrAbove=="below" && distanceFromTop==0.) || (belowOrAbove=="above" && distanceFromBottom == 0.) ||
            (twoPointsAtTheSameRadius && distanceFromTop == 0.) || (twoPointsAtTheSameRadiusOnceDone && distanceFromBottom == 0.) || 
            (i == 1 && distanceFromBottom == 0.) || (i == nLayers && distanceFromTop == 0.)
            
            twoPointsAtTheSameRadiusOnceDone  = false

            x = tmpNormalisedRadii[深さインデックス]

            #push!(parameter.ρ,Coefs.C_ρ[i,1]+Coefs.C_ρ[i,2]*x+Coefs.C_ρ[i,3]*x^2+Coefs.C_ρ[i,4]*x^3)
            #push!(parameter.Vpv,Coefs.C_Vpv[i,1]+Coefs.C_Vpv[i,2]*x+Coefs.C_Vpv[i,3]*x^2+Coefs.C_Vpv[i,4]*x^3)
            #push!(parameter.Vph,Coefs.C_Vph[i,1]+Coefs.C_Vph[i,2]*x+Coefs.C_Vph[i,3]*x^2+Coefs.C_Vph[i,4]*x^3)
            #push!(parameter.Vsv,Coefs.C_Vsv[i,1]+Coefs.C_Vsv[i,2]*x+Coefs.C_Vsv[i,3]*x^2+Coefs.C_Vsv[i,4]*x^3)
            #push!(parameter.Vsh,Coefs.C_Vsh[i,1]+Coefs.C_Vsh[i,2]*x+Coefs.C_Vsh[i,3]*x^2+Coefs.C_Vsh[i,4]*x^3)
            #push!(parameter.Qμ,Coefs.C_Qμ[i,1]+Coefs.C_Qμ[i,2]*x+Coefs.C_Qμ[i,3]*x^2+Coefs.C_Qμ[i,4]*x^3)
            #push!(parameter.Qκ,Coefs.C_Qκ[i,1]+Coefs.C_Qκ[i,2]*x+Coefs.C_Qκ[i,3]*x^2+Coefs.C_Qκ[i,4]*x^3)
            #push!(parameter.QμPower,Coefs.C_QμPower[i,1]+Coefs.C_QμPower[i,2]*x+Coefs.C_QμPower[i,3]*x^2+Coefs.C_QμPower[i,4]*x^3)
            #push!(parameter.QκPower,Coefs.C_QκPower[i,1]+Coefs.C_QκPower[i,2]*x+Coefs.C_QκPower[i,3]*x^2+Coefs.C_QκPower[i,4]*x^3)
            #push!(parameter.η,Coefs.C_η[i,1]+Coefs.C_η[i,2]*x+Coefs.C_η[i,3]*x^2+Coefs.C_η[i,4]*x^3)

            parameter.ρ[深さインデックス]=Coefs.C_ρ[i,1]+Coefs.C_ρ[i,2]*x+Coefs.C_ρ[i,3]*x^2+Coefs.C_ρ[i,4]*x^3
            parameter.Vpv[深さインデックス]=Coefs.C_Vpv[i,1]+Coefs.C_Vpv[i,2]*x+Coefs.C_Vpv[i,3]*x^2+Coefs.C_Vpv[i,4]*x^3
            parameter.Vph[深さインデックス]=Coefs.C_Vph[i,1]+Coefs.C_Vph[i,2]*x+Coefs.C_Vph[i,3]*x^2+Coefs.C_Vph[i,4]*x^3
            parameter.Vsv[深さインデックス]=Coefs.C_Vsv[i,1]+Coefs.C_Vsv[i,2]*x+Coefs.C_Vsv[i,3]*x^2+Coefs.C_Vsv[i,4]*x^3
            parameter.Vsh[深さインデックス]=Coefs.C_Vsh[i,1]+Coefs.C_Vsh[i,2]*x+Coefs.C_Vsh[i,3]*x^2+Coefs.C_Vsh[i,4]*x^3
            parameter.Qμ[深さインデックス]=Coefs.C_Qμ[i,1]+Coefs.C_Qμ[i,2]*x+Coefs.C_Qμ[i,3]*x^2+Coefs.C_Qμ[i,4]*x^3
            parameter.Qκ[深さインデックス]=Coefs.C_Qκ[i,1]+Coefs.C_Qκ[i,2]*x+Coefs.C_Qκ[i,3]*x^2+Coefs.C_Qκ[i,4]*x^3
            parameter.QμPower[深さインデックス]=Coefs.C_QμPower[i,1]+Coefs.C_QμPower[i,2]*x+Coefs.C_QμPower[i,3]*x^2+Coefs.C_QμPower[i,4]*x^3
            parameter.QκPower[深さインデックス]=Coefs.C_QκPower[i,1]+Coefs.C_QκPower[i,2]*x+Coefs.C_QκPower[i,3]*x^2+Coefs.C_QκPower[i,4]*x^3
            parameter.η[深さインデックス]=Coefs.C_η[i,1]+Coefs.C_η[i,2]*x+Coefs.C_η[i,3]*x^2+Coefs.C_η[i,4]*x^3



            if 深さインデックス < N深さインデックス
                深さインデックス += 1
                distanceFromBottom = tmpNormalisedRadii[深さインデックス]-normalisedtmpBottomRadius[i]
                distanceFromTop = tmpNormalisedRadii[深さインデックス]-normalisedtmpTopRadius[i]
                negativeWhenItsIn=distanceFromBottom*distanceFromTop
            else
                negativeWhenItsIn = 100000.
                i += nLayers
            end

            if 深さインデックス < N深さインデックス
                if tmpNormalisedRadii[深さインデックス] == tmpNormalisedRadii[深さインデックス+1]
                    twoPointsAtTheSameRadius = true
                else
                    twoPointsAtTheSameRadius = false
                end
            end

            if distanceFromTop == 0.
                twoPointsAtTheSameRadiusOnceDone = true
                negativeWhenItsIn = 100000.
            end

        else
            
            i += 1
            
        end
    end

    return tmpRadiiInKilometer, parameter
end



function writeClassicplanet1DPSVmodel(Coefs::planet1DPSVmodel,filename::String)
    # this needs 'using Printf'
    fmtDepth(x) = @sprintf("%.2f",x)
    fmtCoefs(x) = @sprintf("%.5f",x)
    fmtInt(x) = @sprintf("%d",x)
    io = open(filename,"w")
    write(io,fmtInt(Coefs.nzone)*'\n')
    
    for i in 1:Coefs.nzone
        write(io,fmtDepth(Coefs.bottomRadius[i])*' '*fmtDepth(Coefs.topRadius[i])*'\n')
        write(io,fmtCoefs(Coefs.C_ρ[i,1])*' '*fmtCoefs(Coefs.C_ρ[i,2])*' '*fmtCoefs(Coefs.C_ρ[i,3])*' '*fmtCoefs(Coefs.C_ρ[i,4])*'\n')
        write(io,fmtCoefs(Coefs.C_Vpv[i,1])*' '*fmtCoefs(Coefs.C_Vpv[i,2])*' '*fmtCoefs(Coefs.C_Vpv[i,3])*' '*fmtCoefs(Coefs.C_Vpv[i,4])*'\n')
        write(io,fmtCoefs(Coefs.C_Vph[i,1])*' '*fmtCoefs(Coefs.C_Vph[i,2])*' '*fmtCoefs(Coefs.C_Vph[i,3])*' '*fmtCoefs(Coefs.C_Vph[i,4])*'\n')
        write(io,fmtCoefs(Coefs.C_Vsv[i,1])*' '*fmtCoefs(Coefs.C_Vsv[i,2])*' '*fmtCoefs(Coefs.C_Vsv[i,3])*' '*fmtCoefs(Coefs.C_Vsv[i,4])*'\n')
        write(io,fmtCoefs(Coefs.C_Vsh[i,1])*' '*fmtCoefs(Coefs.C_Vsh[i,2])*' '*fmtCoefs(Coefs.C_Vsh[i,3])*' '*fmtCoefs(Coefs.C_Vsh[i,4])*'\n')
        write(io,fmtCoefs(Coefs.C_η[i,1])*' '*fmtCoefs(Coefs.C_η[i,2])*' '*fmtCoefs(Coefs.C_η[i,3])*' '*fmtCoefs(Coefs.C_η[i,4])*'\n')
        write(io,fmtCoefs(Coefs.C_Qμ[i,1])*' 'fmtCoefs(Coefs.C_Qκ[i,1])*'\n')
    end
    write(io,"end"*'\n')
    close(io)
end


#export writeClassicplanet1DPSVmodel