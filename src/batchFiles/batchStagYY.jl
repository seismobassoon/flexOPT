


fieldcenStagYY=("t    ","vp   ","eta  ","c    ","f    ","fc   ","df   ","vm   ","ed   ","str  ", "rho  ","k    ","vd   ","ah   ","H    ","nrh  ","ph   ","ly   ","bs   ","hz   ", "TTG  ","bascr","TTGcr","cc   ","air  ","prm  ","hpe  ","imp  ","wtr  ","crb  ","ntg  ","mtl  ","FeO  ","age  ","disl ","gbs  ","Peier","plast","gs   ","defmd","Cp   ","expan","Pdyn ","tpc  ","trho ","Wsol ","Bvel ","Pvel ","Svel ","phs11","phs21","phs31","phs41","phs51","phs12","phs22","phs32","phs42","phs52")
fieldnameStagYY=("Temperature         ","Pressure            ","Viscosity           ","Composition         ","Melt_fraction       ","Melt_composition    ","dmelt/dt            ","Melt_velocity       ","Strain_rate         ","Stress              ","Density             ","Thermal_conductivity","Viscous_dissipation ","Adiabatic_heating   ","Internal_heating    ","Density-normalised  ","Phase               ","Lyapunov_exponent   ","Basalt              ","Harzburgite         ","TTG                 ","Basaltic_crust      ","TTG_crust           ","Continental_crust   ","Air                 ","Primordial          ","HPE                 ","Impactor            ","Water               ","Carbon              ", "Nitrogen            ","Metal               ","FeO                 ","Age                 ","Disl_creep_fraction ","GBS_creep_fraction  ","Peierls_creep_fract.","Plasticity_fraction ","Grain_size          ","Deformation_mode    ","Heat_capacity       ","Thermal_expansivity ","Dynamic_Pressure    ","Tracers_per_cell    ","Tracer-based density","Water_solubility    ","Bulk_sound_velocity ","P-wave_velocity     ","S-wave_velocity     ","Phase_Ol_1          ","Phase_Ol_2          ","Phase_Ol_3          ","Phase_Ol_4          ","Phase_Ol_5          ","Phase_PxGt_1        ", "Phase_PxGt_2        ","Phase_PxGt_3        ","Phase_PxGt_4        ","Phase_PxGt_5        ")


fileExt2FieldNameStagYY = Dict(fieldcenStagYY .=> fieldnameStagYY)
fieldName2FileExtStagYY = Dict(fieldnameStagYY .=> fieldcenStagYY)


function binRead(io,dataTypeTmp::DataType)

    value=read(io,dataTypeTmp)
        
    return value
end

function binRead(io,dataTypeTmp::DataType,numberOfElements)

    values=[]
    for i in 1:numberOfElements
        value=binRead(io,dataTypeTmp)
        values=push!(values,value)
    end
    return values
end


function myListDir(dir;pattern=nothing)
    list=readdir(dir)
    newlist=[]
    if pattern != nothing
        for element in list
            if match(pattern,element) != nothing
                newlist = push!(newlist,dir*"/"*element)
            end
        end
    else 
        newlist=list
    end
    
    return newlist
end

function read_magic(filename;closeFile = false)

    inputILEN = 4
    byte_reverse_in = false

    # Open file in binary mode
    io = open(filename, "r")

    # Read the first 4 bytes as a 32-bit integer
    magic = binRead(io, Int32)

    if !valid_magic(magic) && magic > 8000  # Try 64-bit
        magic -= 8000
        if valid_magic(magic)
            binRead(io, Int32)  # Read an extra 4 bytes (dummy)
            inputILEN = 8
        end
    end

    

    if !valid_magic(magic)  # Try 32-bit with byte reversal
        byte_reverse_in = true
        close(io)
        io = open(filename, "r")
        magic = bswap(binRead(io, Int32))  # Swap bytes for endianness
    end

    if !valid_magic(magic) && magic > 8000  # Try 64-bit with byte reversal
        magic -= 8000
        byte_reverse_in = true
        if valid_magic(magic)
            binRead(io, Int32)   # Read an extra 4 bytes (dummy)
            inputILEN = 8
        end
    end

    if closeFile
        close(io)
    end
    
    return magic, inputILEN, byte_reverse_in, io
end








# Dummy function for validation (replace with real logic)
function valid_magic(m)
    magicLogic = false
    if (m>=1 && m<=12) || (m>=301 && m<=312) || (m>=401 && m<=412) || (m>=101 && m<=112) 
        magicLogic=true
    end

    return magicLogic   # Replace with real valid values
end

function evaluate_nval_from_magicNumber(magic)
    nval=1
    if magic > 400
        nval = 4     #  v+p (new files)
    elseif magic > 300
        nval = 3     # v only (old version)
    end
    magic = mod(magic,100)
    return magic, nval
end

function read_header(io,inputILEN,magic)
    # assume integer and real have same precision in input file
    if inputILEN == 4
        intDataType=Int32
        floatDataType=Float32
    elseif inputILEN == 8
        intDataType=Int64
        floatDataType=Float64
    end
    inputILEN
    

    nx,ny,nz=binRead(io,intDataType,3)
    
    nb=1
    if magic>=7
        nb= binRead(io,intDataType)
    end


    asp_ratio=binRead(io,floatDataType,2) # asp_ratio[1] = asp_ratio[2]/ny 

    nnx,nny,nnz=binRead(io,intDataType,3) # partition numbers

    nnb=1
    if magic>=8
        nnb = binRead(io,intDataType)
    end
    
    

    if magic >=2
        tmpzc=binRead(io,floatDataType,2*(nz+1)-1)
        tmpzc=push!(tmpzc,0.0) # 0.0 is a dummy # who did this??

        zc = reshape(tmpzc,(2,nz+1)) # 2 columns because of the staggered grids, we use the first for plotting (150 km of atmosphere)
    else # normally this should not be the case
        dz = 1.0 / nz
        zc=zeros(floatDataType,2,nz+1)
        for iz in 1:nz
            zc[1,iz] = iz*dz
            zc[2,iz] = dz*(iz+0.5)
        end 
    end
    
    if magic >=7 # core radius for cylindrical or spherical geometry
        rcmb= binRead(io,floatDataType)

        if rcmb==-1.
            boolSpherical= false  #cartesian
        else
            boolSpherical= true  # spherical
        end 
    end 


    if boolSpherical
        dx = asp_ratio[1]/nx
        dy = asp_ratio[2]/ny
    else
        dx = asp_ratio[1]/nx * (zc[1,nz+1]-zc[1,1])
        dy = asp_ratio[2]/ny * (zc[1,nz+1]-zc[1,1])
    end 


    if magic >= 3
        iStep = binRead(io,intDataType)
        time =  binRead(io,floatDataType)
    else
        istep = 0
        time = 0.0
    end

    erupta = zeros(floatDataType,2)

    if magic >= 5 
        erupta[1] = binRead(io,floatDataType)
    end
    if magic >= 12
        erupta[2] = binRead(io,floatDataType)
        intruda = binRead(io,floatDataType,2)
        TTGmass = binRead(io,floatDataType,3)
    end
    if magic >= 6
        botT_val = binRead(io,floatDataType)
    end
    T_core=botT_val
    if magic >= 10
        T_core = binRead(io,floatDataType)
    end
    if magic >= 11
        water_budget = binRead(io,floatDataType)
    end
    if magic >= 4
        binRead(io,floatDataType,nx+ny+nz)
    end

    #usefulHeaders = intDataType,floatDataType,dx,dy,nx,ny,nz,nnx,nny,nnz,nb,rcmb,iStep,time,zc,boolSpherical

    return intDataType,floatDataType,dx,dy,nx,ny,nz,nb,nnx,nny,nnz,nnb,rcmb,iStep,time,zc,boolSpherical
end



function readField(io,floatDataType,nval,nx,ny,nz,nb,nnx,nny,nnz,nnb; xyghost = 0)
    
    
    # xyghost should be controled for vector field

    nvali=nval
    nxpn = nx ÷ nnx
    nypn = ny ÷ nny
    nzpn = nz ÷ nnz
    nbpn = nb ÷ nnb
    nppn = nvali * (nxpn+xyghost) * (nypn+xyghost) * nzpn * nbpn
    
    #wk4=zeros(nvali,nxpn+xyghost,nypn+xyghost,nzpn,nbpn) 
    # original: allocate(wk4(nvali,0:nxpn+xyghost-1,0:nypn+xyghost-1,0:nzpn-1,nbpn))
    
    A=zeros(floatDataType,nvali,nx+xyghost+2,ny+xyghost+2,nz+2,nb)
    # original: allocate(a(nval,-1:nx+xyghost,-1:ny+xyghost,-1:nz,nb)) ! nval first, unlike in StagYY

    
    for ibc in 0:nnb-1             
        ib0i = ibc*nbpn
        for izc in 0:nnz-1
            iz0i = izc*nzpn
            for iyc in 0:nny-1
                iy0i = iyc*nypn
                for ixc in 0:nnx-1
                    ix0i = ixc*nxpn
                    tmp1Darray=binRead(io,floatDataType,nppn)
                    wk4 = reshape(tmp1Darray,nvali,nxpn+xyghost,nypn+xyghost,nzpn,nbpn)
                    xygi = xyghost
                    for ib in 1:nbpn #! loop over input points
                        ibl = ib + ib0i
                        for iz in 0:nzpn-1 
                            izl = iz + iz0i
                            for iy in 0:nypn+xygi-1
                                iyl = iy + iy0i
                                for ix in 0:nxpn+xygi-1
                                    ixl = ix + ix0i
                                    A[1:nvali,ixl+2,iyl+2,izl+2,ibl]=wk4[1:nvali,ix+1,iy+1,iz+1,ib]
                                end
                            end
                        end
                    end
                end
                
            end
        end
    end
    #original: (ix=0:nx-1,iy=0:ny-1,iz=0:nz-1,ib=1:1)
    #original: fields(ix,iy,iz,ib,ifi) = a(1,ix,iy,iz,ib)
    newField = A[1:nval,2:nx+1,2:ny+1,2:nz+1,1:nb]
                         

    return newField
end




function quarterDiskExtrapolationRawGrid!(fi, Xnode, Ynode)
    fiO=copy(fi)
    XnodeO=copy(Xnode)
    YnodeO=copy(Ynode)
    append!(fi,fiO)
    append!(Xnode,-YnodeO)
    append!(Ynode,XnodeO)

    append!(fi,fiO)
    append!(Xnode,-XnodeO)
    append!(Ynode,-YnodeO)
    
    append!(fi,fiO)
    append!(Xnode,YnodeO)
    append!(Ynode,-XnodeO)
    return 
end


function quarterDiskExtrapolation(fi,nX,nY)
#it depends but  the discontinuities are sometimes visible
    halfnX = (nX-1) ÷ 2
    halfnY = (nY-1) ÷ 2

    fi[1:halfnX+1,1:halfnY+1]=transpose(fi[nX:-1:halfnX+1,1:halfnY+1])
    fi[halfnX+1:nX,halfnY+1:nY]= transpose(fi[halfnX+1:nX,halfnY+1:-1:1])
    fi[1:halfnX+1,halfnY+1:nY]=fi[nX:-1:halfnX+1,halfnY+1:-1:1]

    return fi
end




function readStagYYFiles(file)
    magic, inputILEN, byte_reverse_in, io=read_magic(file)
    magic,nval = evaluate_nval_from_magicNumber(magic)
    intDataType,floatDataType,dx,dy,nx,ny,nz,nb,nnx,nny,nnz,nnb,rcmb,iStep,time,zc,boolSpherical=read_header(io,inputILEN,magic)
   
    nodes=nothing
    boolFlat = true

    Xnode=[]
    Ynode=[]    
    Znode=[]


   
    if nx == 1
       boolFlat = true
    end
 


    if boolSpherical
        # these are the grid points including the boundaries
        r = rcmb .+ zc[1,1:nz+1]
        
        theta = (collect(1:1:nx+1) .- (nx+1)) .* dx .+ 0.5*π
        phi = (collect(1:1:ny+1) .- (ny+1)) .* dy
        

        # these are the midpoints 
        r_mid = rcmb .+ zc[2,1:nz]
        theta_mid = (collect(1:1:nx) .- (nx+1)) .* dx .+ 0.5*π .+ 0.5*dx
        phi_mid = phi = (collect(1:1:ny) .- (ny+1)) .* dy .+ 0.5*dy
        
        # interpolation should be done with the midpoints + end points

        r_new = prepend!(r_mid,r[1])
        r_new = push!(r_mid,r[end])
        minR = r[1]
        maxR = r[end]
        rC=r_new

        theta_new = prepend!(theta_mid,theta[1])
        theta_new = push!(theta_new,theta[end])
        minθ = theta[1]
        maxθ = theta[end]
        phi_new = prepend!(phi_mid,phi[1])
        phi_new = push!(phi_new,phi[end])
        minϕ = phi[1]
        maxϕ = phi[end]
        if boolFlat
            minθ = 0.
            maxθ = 0.
            theta_new = (π/2)  # if flat, we are looking at the 
            nodes=(theta_new,phi_new,r_new)
        else   
            nodes=(theta_new,phi_new,r_new)
        end

    
    else
        # do not know if this works (certainly not! try to mimic the code above for the spherical case)
        x=collect(1:1:nx+1) .* dx
        y=collect(1:1:ny+1) .* dy
        z=zc
        nodes=(x,y,z)
    end

    

    theta_new_piCoef = theta_new ./ π
    phi_new_piCoef = phi_new ./ π

    if boolSpherical
        for rTemp in r_new
            for phiTemp in phi_new_piCoef
                cosPhiPi = cospi(phiTemp)
                sinPhiPi = sinpi(phiTemp)
               
                for thetaTemp in theta_new_piCoef
                    cosThetaPi = cospi(thetaTemp)
                    sinThetaPi = sinpi(thetaTemp) 
                    Xnode = push!(Xnode,rTemp*sinThetaPi*cosPhiPi)
                    Ynode = push!(Ynode,rTemp*sinThetaPi*sinPhiPi)
                    Znode = push!(Znode,rTemp*cosThetaPi)

                end
            end
        end
    else
        @error "cartesian format: not yet developed"
    end

    


   #coord[1,ix,iy,iz] = r[iz]*sin(theta)*cos(phi)
   #coord[2,ix,iy,iz] = r[iz]*sin(theta)*sin(phi)
   #coord[3,ix,iy,iz] = r[iz]*cos(theta)
   
   
   

   #field = binRead(io,floatDataType,(nx)*(ny)*(nz))
   #field = binRead(io,Float64,nx*ny*nz)
   #rawField =field
   #newField = nothing

   rawField = readField(io,floatDataType,nval,nx,ny,nz,nb,nnx,nny,nnz,nnb)

   
    if boolFlat
        #field=reshape(field, (ny,nz)) 

        field = zeros(floatDataType,ny,nz)
        field[1:ny,1:nz]=rawField[1,1,1:ny,1:nz,1]


        newField = zeros(floatDataType,ny+2,nz+2)

        # the interior

        newField[2:ny+1,2:nz+1] = field[1:end,1:end]

        # 4 'surfaces'

        newField[1,2:nz+1] = field[1,1:nz]
        newField[end,2:nz+1] = field[end,1:nz]
        newField[2:ny+1,1] = field[1:ny,1]
        newField[2:ny+1,end] = field[1:ny,end]

        # 4 'endpoints'

        newField[1,1,1]=field[1,1,1]
        newField[1,end,1]=field[1,end,1]
        newField[1,1,end]=field[1,1,end]
        newField[1,end,end]=field[1,end,end]


    else
        #field=reshape(field, (nx,ny,nz)) 
        field = zeros(floatDataType,nx,ny,nz)
        field[1:ny,1:nz]=rawField[1,1:nx,1:ny,1:nz,1]

        newField = zeros(floatDataType,nx+2,ny+2,nz+2)

        # the interior

        newField[2:nx+1,2:ny+1,2:nz+1] = field[1:nx,1:ny,1:nz]

        # 6 'surfaces'
        
        newField[1,2:ny+1,2:nz+1] = field[1,1:ny,1:nz]
        newField[end,2:ny+1,2:nz+1] = field[end,1:ny,1:nz]
        newField[2:nx+1,1,2:nz+1] = field[1:nx,1,1:nz]
        newField[2:nx+1,end,2:nz+1] = field[1:nx,end,1:nz]
        newField[2:nx+1,2:ny+1,1] = field[1:nx,1:ny,1]
        newField[2:nx+1,2:ny+1,end] = field[1:nx,1:ny,end]

        # 8 'endpoints'
        newField[1,1,1]=field[1,1,1]
        newField[end,1,1]=field[end,1,1]
        newField[1,end,1]=field[1,end,1]
        newField[end,end,1]=field[end,end,1]
        newField[1,1,end]=field[1,1,end]
        newField[end,1,end]=field[end,1,end]
        newField[1,end,end]=field[1,end,end]
        newField[end,end,end]=field[end,end,end]

    end

    if boolFlat
        newField=reshape(newField,(ny+2)*(nz+2))
        return newField, Xnode, Ynode, rcmb
    else
        newField=reshape(newField,(nx+2)*(ny+2)*(nz+2))
        return newField, Xnode, Ynode, Znode, rcmb
    end
   #fieldInterpolated=interpolate(nodes,newField,Gridded(Linear()))

end




function readStagYYFilesAverage(file)
    magic, inputILEN, byte_reverse_in, io=read_magic(file)
    magic,nval = evaluate_nval_from_magicNumber(magic)
    intDataType,floatDataType,dx,dy,nx,ny,nz,nb,nnx,nny,nnz,nnb,rcmb,iStep,time,zc,boolSpherical=read_header(io,inputILEN,magic)
   
    nodes=nothing
    boolFlat = true

    Xnode=[]
    Ynode=[]    
    Znode=[]


   
    if nx == 1
       boolFlat = true
    end
 


    if boolSpherical
        # these are the grid points including the boundaries
        r = rcmb .+ zc[1,1:nz+1]
        
        theta = (collect(1:1:nx+1) .- (nx+1)) .* dx .+ 0.5*π
        phi = (collect(1:1:ny+1) .- (ny+1)) .* dy
        

        # these are the midpoints 
        r_mid = rcmb .+ zc[2,1:nz]
        theta_mid = (collect(1:1:nx) .- (nx+1)) .* dx .+ 0.5*π .+ 0.5*dx
        phi_mid = phi = (collect(1:1:ny) .- (ny+1)) .* dy .+ 0.5*dy
        
        # interpolation should be done with the midpoints + end points

        r_new = prepend!(r_mid,r[1])
        r_new = push!(r_mid,r[end])
        minR = r[1]
        maxR = r[end]
        rC=r_new

        theta_new = prepend!(theta_mid,theta[1])
        theta_new = push!(theta_new,theta[end])
        minθ = theta[1]
        maxθ = theta[end]
        phi_new = prepend!(phi_mid,phi[1])
        phi_new = push!(phi_new,phi[end])
        minϕ = phi[1]
        maxϕ = phi[end]
        if boolFlat
            minθ = 0.
            maxθ = 0.
            theta_new = (π/2)  # if flat, we are looking at the 
            nodes=(theta_new,phi_new,r_new)
        else   
            nodes=(theta_new,phi_new,r_new)
        end

    
    else
        # do not know if this works (certainly not! try to mimic the code above for the spherical case)
        x=collect(1:1:nx+1) .* dx
        y=collect(1:1:ny+1) .* dy
        z=zc
        nodes=(x,y,z)
    end

    

    theta_new_piCoef = theta_new ./ π
    phi_new_piCoef = phi_new ./ π

    if boolSpherical
        for rTemp in r_new
            for phiTemp in phi_new_piCoef
                cosPhiPi = cospi(phiTemp)
                sinPhiPi = sinpi(phiTemp)
               
                for thetaTemp in theta_new_piCoef
                    cosThetaPi = cospi(thetaTemp)
                    sinThetaPi = sinpi(thetaTemp) 
                    Xnode = push!(Xnode,rTemp*sinThetaPi*cosPhiPi)
                    Ynode = push!(Ynode,rTemp*sinThetaPi*sinPhiPi)
                    Znode = push!(Znode,rTemp*cosThetaPi)

                end
            end
        end
    else
        @error "cartesian format: not yet developed"
    end

    


   #coord[1,ix,iy,iz] = r[iz]*sin(theta)*cos(phi)
   #coord[2,ix,iy,iz] = r[iz]*sin(theta)*sin(phi)
   #coord[3,ix,iy,iz] = r[iz]*cos(theta)
   
   
   

   #field = binRead(io,floatDataType,(nx)*(ny)*(nz))
   #field = binRead(io,Float64,nx*ny*nz)
   #rawField =field
   #newField = nothing

   rawField = readField(io,floatDataType,nval,nx,ny,nz,nb,nnx,nny,nnz,nnb)

   
    if boolFlat
        #field=reshape(field, (ny,nz)) 

        field = zeros(floatDataType,ny,nz)
        field[1:ny,1:nz]=rawField[1,1,1:ny,1:nz,1]


        newField = zeros(floatDataType,ny+2,nz+2)

        # the interior

        newField[2:ny+1,2:nz+1] = field[1:end,1:end]

        # 4 'surfaces'

        newField[1,2:nz+1] = field[1,1:nz]
        newField[end,2:nz+1] = field[end,1:nz]
        newField[2:ny+1,1] = field[1:ny,1]
        newField[2:ny+1,end] = field[1:ny,end]

        # 4 'endpoints'

        newField[1,1,1]=field[1,1,1]
        newField[1,end,1]=field[1,end,1]
        newField[1,1,end]=field[1,1,end]
        newField[1,end,end]=field[1,end,end]


    else
        #field=reshape(field, (nx,ny,nz)) 
        field = zeros(floatDataType,nx,ny,nz)
        field[1:ny,1:nz]=rawField[1,1:nx,1:ny,1:nz,1]

        newField = zeros(floatDataType,nx+2,ny+2,nz+2)

        # the interior

        newField[2:nx+1,2:ny+1,2:nz+1] = field[1:nx,1:ny,1:nz]

        # 6 'surfaces'
        
        newField[1,2:ny+1,2:nz+1] = field[1,1:ny,1:nz]
        newField[end,2:ny+1,2:nz+1] = field[end,1:ny,1:nz]
        newField[2:nx+1,1,2:nz+1] = field[1:nx,1,1:nz]
        newField[2:nx+1,end,2:nz+1] = field[1:nx,end,1:nz]
        newField[2:nx+1,2:ny+1,1] = field[1:nx,1:ny,1]
        newField[2:nx+1,2:ny+1,end] = field[1:nx,1:ny,end]

        # 8 'endpoints'
        newField[1,1,1]=field[1,1,1]
        newField[end,1,1]=field[end,1,1]
        newField[1,end,1]=field[1,end,1]
        newField[end,end,1]=field[end,end,1]
        newField[1,1,end]=field[1,1,end]
        newField[end,1,end]=field[end,1,end]
        newField[1,end,end]=field[1,end,end]
        newField[end,end,end]=field[end,end,end]

    end

    if boolFlat
        avNewField = similar(newField)
        diffNewField = similar(newField)
        for iz in 1:nz+2
            average = sum(newField[:,iz])
            average /= Float64(ny+2)
            avNewField[:,iz] .= average
        end
        #@show size(avNewField)
        diffNewField = newField .- avNewField
        newField=reshape(newField,(ny+2)*(nz+2))
        avNewField=reshape(avNewField,(ny+2)*(nz+2))
        diffNewField=reshape(diffNewField,(ny+2)*(nz+2))
        return newField, avNewField, diffNewField, Xnode, Ynode, rcmb
    else
        newField=reshape(newField,(nx+2)*(ny+2)*(nz+2))
        return newField, Xnode, Ynode, Znode, rcmb
    end
   #fieldInterpolated=interpolate(nodes,newField,Gridded(Linear()))

end


function extendToCoreWithρ!(ρfield, Xnode, Ynode, rcmb, dR; dθ=2*π/360.0, iCheckCoreModel=true)
    # local function here: this requires planet1D.jl, testparam.csv
    #
    # This function will add the ρ field computed only for the core 

    # Xnode and Ynodes are the one that we give but 

    # the 1D core model will be given by specifying the file to use in testparam.csv 


    premCMB = planet1D.my1DDSMmodel.averagedPlanetCMBInKilometer * 1.e3

    arrayRadius = collect(0:dR:rcmb)
    if arrayRadius[end] != rcmb
        arrayRadius = push!(arrayRadius,rcmb)
    end
    arrayRadiusFaked = min.(arrayRadius,premCMB) .* 1.e-3



    _, arrayParams  = planet1D.compute1DseismicParamtersFromPolynomialCoefficientsWithGivenRadiiArray(planet1D.my1DDSMmodel, arrayRadiusFaked, "below")
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
