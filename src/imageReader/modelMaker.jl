
function defineModel(file::String;modelDefinitionMethod="2DimageFile",colormap="jet")

    model = nothing

    if modelDefinitionMethod !== nothing
            
        #region Model input - option i) Model domain definition

        if modelDefinitionMethod === "ToyModel"
            DomainWindow=(DomainWindowT=1.0,DomainWindowX=1.0,DomainWindowY=1.0,DomainWindowZ=1.0)
            ModelSizeTXYZ=(ModelSizeT=101,ModelSizeX=101,ModelSizeY=101,ModelSizeZ=0)
        end

        #endregion

        #region Model input - option ii) Read a file (2D or 3D) and define Δs

        if modelDefinitionMethod === "2DimageFile"
            imagefile = file
            #colormap = "jet" #colormap can be RGB vector or predefined colormap

            #model=read2DimageModel(imagefile,colormap;Nwidth=10,Nheight=10,showRecoveredImage=false)
            model=read2DimageModel(imagefile,colormap;showRecoveredImage=false)
        end
        #endregion

        #region Model input - option iii) Read a file (1D spherical planet models)

        if modelDefinitionMethod ==="1DsphericalPlanet"
            # use some programmes that are developed during Xmas 2023
            # inputModels.jl
        end

    #endregion
    end
    return model
end

