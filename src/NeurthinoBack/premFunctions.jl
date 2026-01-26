# to have a PREM version to compare with the model
#PREM version of everything

function extendWithρ!(ρfield, Xnode, Ynode, dR;rcmb=6.371e6, dθ=2*π/360.0, iCheckCoreModel=true)
    # local function here: this requires planet1D.jl, testparam.csv

    premCMB = planet1D.my1DDSMmodel.averagedPlanetCMBInKilometer * 1.e3

    arrayRadius = collect(0:dR:rcmb)
    if arrayRadius[end] != rcmb
        arrayRadius = push!(arrayRadius,rcmb)
    end


    _, arrayParams  = planet1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(planet1D.my1DDSMmodel, arrayRadius.*1e-3, "below")
    #planet1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(planet1D.my1DDSMmodel, arrayRadius.*1.e-3, "above")
    tmpDensity=arrayParams.ρ


    tmpXnode = [(tmpRadius*cos(tmpθ)) for tmpRadius in arrayRadius for tmpθ in collect(0:dθ:2π)]
    tmpYnode = [(tmpRadius*sin(tmpθ)) for tmpRadius in arrayRadius for tmpθ in collect(0:dθ:2π)]
    tmpValue = [(1.e3*tmpDensity[iRadius]) for iRadius in eachindex(arrayRadius) for tmpθ in collect(0:dθ:2π)]
    #@show minimum(tmpXnode),maximum(tmpXnode)

    if iCheckCoreModel
        f=Figure()
    
        lines(f[1,1],arrayRadius, arrayParams.ρ,color=:red)

        display(f)
    end


    Xnode=append!(Xnode,tmpXnode)
    Ynode=append!(Ynode,tmpYnode)
    ρfield=append!(ρfield,tmpValue)

    return
    
end




function myPREMPlot2DConvectionModel(iTime, fieldname, filename)
#only if the field in DIVandrun is the same as in readStagYYFiles
    file = filename[iTime]
    field, Xnode, Ynode, _ = readStagYYFiles(file)
    extendWithρ!(field, Xnode, Ynode, dR, iCheckCoreModel=false)
    #quarterDiskExtrapolationRawGrid!(field, Xnode, Ynode)
    fi,_ = DIVAndrun(mask,(pm,pn),(xi,yi),(Xnode,Ynode),field,correlationLength,epsilon2);
    
    diam = maxX - minX
    x = range(0, diam, length=size(fi)[1])
    y = range(0, diam, length=size(fi)[2])

    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1)
    colormap = myChoiceColormap(fieldname)
    hm=heatmap!(ax,x, y, fi, colormap=colormap)
    Colorbar(fig[:,2], hm, label="wtr (wt.%)")

    return fig, ax, fi
end




function PREMlineDensityElectron2D(n_pts, iTime, positionDetector, NeutrinoSource, colorname, ax1, dR)
    #draw a line between positionDetector and NeutrinoSource (coordinates) and give the density/distance profile
    #dependencies : Makie

    fig, ax, fi = myPREMPlot2DConvectionModel(iTime, "rho", rhoFiles)
    x_phys = range(positionDetector[1], NeutrinoSource[1], length=n_pts)
    y_phys = range(positionDetector[2], NeutrinoSource[2], length=n_pts)  
    
    lines!(ax, x_phys,y_phys, color=colorname)  # (x,y)_phys in m
    display(fig)

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

function PREMcorrectedPosition(x,y, zposition; center = [6.5e6, 6.5e6], earth_radius = 6.371e6)
    #to place the detector precisely on the surface, then the detector is buried for zposition

    dx = x - center[1]
    dy = y - center[2]

    real_pos = earth_radius - zposition #m
    dist_radiale = sqrt(dx^2 + dy^2)
    new_x = center[1] + real_pos*dx/dist_radiale
    new_y = center[2] + real_pos*dy/dist_radiale
    return new_x, new_y, zposition
end


function solveQuadraticEquation(a,b,c)
    Δ = b^2 - 4*a*c

    if Δ>0
        x1 = ((-b - sqrt(Δ))/(2*a))
        x2 = ((-b + sqrt(Δ))/(2*a))
    else
        x1 = -b/(2 *a)
        x2 = x1
    end
    return x1, x2

end

function posOrNeg(cos_θ, sign = :positive)
    #to choose if we want a positive or a negative angle
    if sign == :positive
        θ = acos.(cos_θ)
    else
        θ = .- acos.(cos_θ)
    end
    return θ
end


function PREMsourcePosition(center, positionDetector, n_vectors, zposition; earthRadius = 6.371e6)
    #to get the position of the different sources

    (xc, yc) = center[1], center[2] #m 
    (xd, yd) = positionDetector[1], positionDetector[2] #m
    XY = []

    cos_θ = range(-1, 0, length = n_vectors)
    θ = posOrNeg(cos_θ, :positive)
    cos_epi = cos.(2 .*θ .- π)
    sin_epi = sin.(2 .*θ .- π)
    rotation = cos_epi .+ im .* sin_epi


    for i in eachindex(cos_θ)
        equ = ((xd - xc) + (yd-yc)*im) * rotation[i]
        X = real(equ)+xc
        Y = imag(equ)+yc
        # this is the position with zposition below

        newX = nothing
        newY = nothing
        if cos_θ[i] !== 0.0
            if X-xd != 0.0
                slope = (Y-yd)/(X-xd)

                a = 1+slope^2
                b = -2*xc + 2*Y*slope-2*slope^2*X-2*slope*yc
                c = xc^2 + Y^2 - 2*Y*slope*X + slope^2*X^2 -2*Y*yc + 2*slope*X*yc +yc^2 - earthRadius^2
                sol1,sol2 = solveQuadraticEquation(a,b,c)

                if (sol1-xd)*(X-xd)>0.0
                    newX = sol1
                    newY = Y + slope* (newX - X)
                else
                    newX = sol2
                    newY = Y + slope* (newX - X)
                end

            else
                slope = (X-xd)/(Y-yd)
                
                a = 1+slope^2
                b = -2*yc + 2*X*slope-2*slope^2*Y-2*slope*xc
                c = yc^2 + X^2 - 2*X*slope*Y + slope^2*Y^2 -2*X*xc + 2*slope*Y*xc +xc^2 - earthRadius^2
                sol1,sol2 = solveQuadraticEquation(a,b,c)
                
                if (sol1-yd)*(Y-yd)>0.0
                    newY = sol1
                    newX = Y + slope* (newY - Y)
                else
                    newY = sol2
                    newX = Y + slope* (newY - Y)
                end

            end

        else
            segmentfromDtoS = sqrt(earthRadius^2-(earthRadius-zposition)^2)
            if θ[i] > 0
                newX = xd - (yd-yc)/(earthRadius-zposition)*segmentfromDtoS
                newY = yd + (xd-xc)/(earthRadius-zposition)*segmentfromDtoS
            else 
                newX = xd + (yd-yc)/(earthRadius-zposition)*segmentfromDtoS
                newY = yd - (xd-xc)/(earthRadius-zposition)*segmentfromDtoS
            end
        end
        push!(XY, (newX,newY))

    end
    return XY

end

function PREMvectorsFromDetector(n_vectors, zposition ;center = [6.5e6, 6.5e6])
    #draw n_vectors (diff θ) from a detector (placed by interaction) and return density profiles for each vector through the Earth
    #dependencies : GLMakie, Makie, Colors


    pos = Float32[1.8340824f6, 7.7562185f6] #you have to put manually the position of the detector to compare
    @show pos
    x, y = pos[1], pos[2]
    new_x, new_y,zposition = PREMcorrectedPosition(x,y, zposition) 
    XY = PREMsourcePosition((center[1], center[2]), (new_x, new_y), n_vectors, zposition)

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
        dens, section = PREMlineDensityElectron2D(n_pts, iTime, detector,source, colorname, ax1, dR)

        push!(densities_list, dens)
        push!(sections_list, section)

    end

    display(fig1)
    return densities_list, sections_list
end

function roundExt(x,step)
    #to round the values of H and U
    if x isa Complex
        xreal = round(real(x)/step)*step
        ximag = round(imag(x)/step)*step
        return complex(xreal, ximag)
    else
        return round(x/step)*step
    end
end


function PREMcreationPaths(n_vectors, zposition)

    densities_list, sections_list = PREMvectorsFromDetector(n_vectors, zposition) 
    paths = Vector{Path}(undef, n_vectors)  

    for i in eachindex(paths)
        paths[i]= Path(densities_list[i],sections_list[i])
    end

    return paths
end


function linkWithNeurthinoPREM()
    osc = OscillationParameters(3)
    setθ!(osc, 1=>2, 0.59)
    setθ!(osc, 1=>3, 0.15)
    setθ!(osc, 2=>3, 0.84)
    setδ!(osc, 1=>3, 3.86)
    setΔm²!(osc, 2=>3, -2.523e-3)
    setΔm²!(osc, 1=>2, -7.39e-5)
    U = PMNSMatrix(osc)
    H = Hamiltonian(osc)
    Uround = roundExt.(U, 0.01)
    Hround = roundExt.(H, 0.00001)
    cos_θ = range(-1, 0, length = n_vectors)

    paths = PREMcreationPaths(n_vectors, zposition)
    energies = 10 .^ range(0, stop=2, length=n_vectors)
    probs2 = Pνν(Uround, Hround, energies, paths)[:,:,1,2]
    matprobs2=parent(probs2)

    fig = Figure()
    ax = Axis(fig[1,1], aspect = 1, xscale=log10, xlabel="Energy (GeV)", ylabel="cos(θ)")
    hm=heatmap!(ax, energies, cos_θ, matprobs2, colormap=cgrad(:inferno))
    Colorbar(fig[:,2], hm, label="Probability")
    display(fig)

    return probs2
end


