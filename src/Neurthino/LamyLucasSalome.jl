# these are some developments that NF has been doing to deal with 2D/3D heterogeneous Earth models 
# with the aid of/inpired by Estelle Salomé, Mahé Lucas, Baptiste Lamy 
# during the summer 2025 and the winter 2026





function ZOverAwithWaterWithThreeLayeredZA1Dmodel(wtrfieldCartesien,xCartesien,yCartesien)#avec le tableau de Lucas (p.4 Unveiling the Outer core with not)
    ZoverAwtr=similar(wtrfieldCartesien)
    ZoverAwtr.= 0.0
    ZoverA=similar(wtrfieldCartesien)
    ZoverA.= 0.0
    for i in eachindex(xCartesien), j in eachindex(yCartesien)
        x=xCartesien[i]
        y=yCartesien[j]
        r=sqrt(x*x+y*y)
        if 0.0 ≤ r < 1221000.5
            ZoverA[i,j]=0.466
        elseif  1221000.5 <= r < 3480000.0
            ZoverA[i,j]=0.466
        elseif 3480000.0 <= r < 6400000.0
            ZoverA[i,j]=0.496
        end
        #if 0.0 ≤ r < 6368000.0 #a rajouter si on veut que l'espace ai un Z/A de 0
            ZoverAwtr[i,j] = wtrfieldCartesien[i,j]*(5.0/9.0)+ (1.0-wtrfieldCartesien[i,j])*ZoverA[i,j]
        #end
     end

    
    return ZoverAwtr,ZoverA    
end




function lineDensityElectron2D(n_pts, effectiveρModel, positionDetector, NeutrinoSource, colorname, ax1, dR)
    #draw a line between positionDetector and NeutrinoSource (coordinates) and give the density/distance profile
    #dependencies : Makie

    #fig, ax, fi = myPREMPlot2DConvectionModel(iTime, "rho", rhoFiles)
    fi = effectiveρModel

    x_phys = range(positionDetector[1], NeutrinoSource[1], length=n_pts)
    y_phys = range(positionDetector[2], NeutrinoSource[2], length=n_pts)  
    
    #lines!(ax, x_phys,y_phys, color=colorname)  # (x,y)_phys in m
    #display(fig)

    x_grid = x_phys ./dR
    y_grid = y_phys ./dR
    itp = interpolate(fi, BSpline(Linear()), OnGrid())

    densGrids = Float64[]
    for i in eachindex(x_grid)
        x = x_grid[i]
        y = y_grid[i]
        push!(densGrids, itp(x,y)*1e-3) #g/cm3
    end

    dens=Float64[]
    for i in eachindex(densGrids)[1:end-1]
        push!(dens, 0.5*(densGrids[i]+densGrids[i+1]))
    end

    segmentLengthInKm = sqrt((x_phys[2]-x_phys[1])^2 + (y_phys[2]-y_phys[1])^2) * 1.e-3
    sections = segmentLengthInKm .* ones(Float64,n_pts-1) 
    dist = segmentLengthInKm*collect(0:1:n_pts-1) #km

    lines!(ax1, dist, densGrids, color=colorname)
    return dens, sections
end


function creationPaths(n_vectors, pos,zposition,effectiveρModel;center = [6.5e6, 6.5e6])

    densities_list, sections_list = vectorsFromDetector(n_vectors, pos,zposition,effectiveρModel; center=center) 
    paths = Vector{Path}(undef, n_vectors)  

    for i in eachindex(paths)
        paths[i]= Path(densities_list[i],sections_list[i])
    end

    return paths
end

function creationPaths(n_vectors, pos,zposition,effectiveρModel;center = [6.5e6, 6.5e6])

    densities_list, sections_list = vectorsFromDetector(n_vectors, pos,zposition,effectiveρModel; center=center) 
    paths = Vector{Path}(undef, n_vectors)  

    for i in eachindex(paths)
        paths[i]= Path(densities_list[i],sections_list[i])
    end

    return paths
end



function vectorsFromDetector(n_vectors, pos, zposition,effectiveρModel ;center = [6.5e6, 6.5e6])
    #draw n_vectors (diff θ) from a detector (placed by interaction) and return density profiles for each vector through the Earth
    #dependencies : GLMakie, Makie, Colors


    
    #@show pos
    x, y = pos[1], pos[2]
    new_x, new_y,zposition = correctedPosition(x,y, zposition) 
    XY = sourcePosition((center[1], center[2]), (new_x, new_y), n_vectors, zposition)

    segments_pts = []
    for source in XY
        push!(segments_pts, (new_x, new_y))
        push!(segments_pts, (source[1], source[2]))
    end

    
    CairoMakie.activate!()
    fig1 = Figure()
    ax1 = Axis(fig1[1,1], xlabel="Path (km)", ylabel="Density (g/cm3)")


    densities_list = []
    sections_list = []
    for i in eachindex(XY)
        colorname = rand(collect(keys(Colors.color_names)))
        detector = new_x, new_y
        source = XY[i][1], XY[i][2]
        dens, section = lineDensityElectron2D(n_pts, effectiveρModel, detector,source, colorname, ax1, dR)

        push!(densities_list, dens)
        push!(sections_list, section)

    end

    display(fig1)
    return densities_list, sections_list
end
