function myChoiceColormap(fieldname)
    if fieldname ==="temperature"
        colormap=cgrad(:seismic)

    elseif fieldname === "rho" || fieldname === "composition"
        colormap=cgrad(:viridis)

    elseif fieldname === "wtr"
        colormap=cgrad(:blues)

    else 
        colormap=cgrad(:viridis)

    end
end 

function myPlot2DConvectionModel(iTime, fieldname, filename)
#only if the field in DIVandrun is the same as in readStagYYFiles
    file = filename[iTime]
    field, Xnode, Ynode, rcmb = readStagYYFiles(file)
    extendToCoreWithρ!(field, Xnode, Ynode, rcmb, dR, iCheckCoreModel=false)
    #quarterDiskExtrapolationRawGrid!(field, Xnode, Ynode)
    fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
    
    diam = maxX - minX
    x = range(0, diam, length=size(fi)[1])
    y = range(0, diam, length=size(fi)[2])

    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1)
    colormap = myChoiceColormap(fieldname)
    hm=heatmap!(ax,x, y, fi, colormap=colormap)
    Colorbar(fig[:,2], hm, label="Density (kg/m3)")

    return fig, ax, fi
end


# below is for making a video
function myAnimation(step, iTime, fieldname, filename)
    fi_list=[]
    for file in filename[1:step:iTime]
        field, Xnode, Ynode= readStagYYFiles(file)
        #extendToCoreWithρ!(field, Xnode, Ynode, rcmb, dR)
        #quarterDiskExtrapolationRawGrid!(field, Xnode, Ynode)
        fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
        push!(fi_list, fi)
    end

    fig = Figure()
    ax = Axis(fig[1,1],aspect = 1)
    data = Observable(fi_list[1])
    colormap = myChoiceColormap(fieldname)
    hm=heatmap!(ax, data, colormap=colormap)#, colorrange=()) if needed
    Colorbar(fig[:, 2], hm, label="Temperature (K)")
    record(fig, "animation2D.mp4", 1:length(fi_list); framerate = 2) do i
        data[] = fi_list[i]
    end
end