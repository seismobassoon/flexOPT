DEFAULT_PLANET = Ref(:Earth) #set_default_planet! can change this


# Planetary parameters (semi-major axis `a` in meters, flattening `f`)

function planet_ellipsoid(name::Symbol)
    if name == :Earth # WGS84
        a= 6378137.0
        f =1.0/298.257223563
    elseif name == :Mars
        a = 3396200.0
        f = 1 / 169.8
    elseif name == :Moon
        a = 1737400.0
        f = 0.0
    elseif name == :Venus
        a = 6051800.0
        f = 0.0
    elseif name == :Mercury
        a = 2439700.0
        f = 0.0
    else
        error("Unknown planet: $name")
    end
    e2 = 2f - f^2
    b = a * (1 - f)
    return Ellipsoid(a, b, f, e2, name)
end


const DEFAULT_ELLIPSOID = Ref(planet_ellipsoid(DEFAULT_PLANET[]))

# The planet can be changed by this function

function set_default_planet!(name::Symbol)
    DEFAULT_PLANET[]=name
    DEFAULT_ELLIPSOID[] = planet_ellipsoid(name)
end


struct GeoPoint
    lat::Float64 # in degree
    lon::Float64 # in degree
    alt::Float64 # in metre
    ecef::SVector{3,Float64}
    radius::Float64 # in metre
    #effectiveRadius::Float64 # in metre (useful to get values from 1D averaged model)
end

struct localCoord2D
    iXZ::SVector{2,Integer}
    xz::SVector{2,Float64}
    horizontalVector::SVector{2,Float64}
    normalVector::SVector{2,Float64}
end

struct localCoord3D
    iXYZ::SVector{3,Integer}
    xyz::SVector{3,Float64}
    horizontalVector1::SVector{3,Float64}
    horizontalVector2::SVector{3,Float64}
    normalVector::SVector{3,Float64}
end

# 2D <-> 3D conversion

function localCoord2D(p::localCoord3D)
    iXZ=SVector{2,Integer}(p.iXYZ[1],p.iXYZ[3])
    xz=SVector{2,Float64}(p.xyz[1],p.xyz[3])
    horizontalVector=SVector{2,Float64}(p.horizontalVector1[1],p.horizontalVector1[3])
    normalVector=SVector{2,Float64}(p.normalVector[1],p.normalVector[3])
    localCoord2D(iXZ,xz,horizontalVector,normalVector)
end

function localCoord3D(p::localCoord2D)
    iXYZ=SVector{3,Integer}(p.iXYZ[1],1,p.iXYZ[2])
    xyz=SVector{3,Float64}(p.xyz[1],0.0,p.xyz[2])
    horizontalVector1=SVector{3,Float64}(p.horizontalVector1[1],0.0,p.horizontalVector1[2])
    horizontalVector2=SVector{3,Float64}(0.0,1.0,0.0)
    normalVector=SVector{3,Float64}(p.normalVector[1],0.0,p.normalVector[2])
    localCoord3D(iXYZ,xyz,horizontalVector1,horizontalVector2,normalVector)
end

# for simple Cartesian coordinates

function localCoord2D(ix::Integer,iz::Integer,Δx::Float64,Δz::Float64; startX::Float64=-Δx,startZ::Float64=-Δz)
    # this gives simple Cartesian coordinates
    iXZ = SVector{2,Integer}(ix,iz)
    xz = iXZ .* SVector{2,Float64}(Δx,Δz) + SVector{2,Float64}(startX,startZ)
    horizontalVector=SVector{2,Float64}(1.0,0.0)
    normalVector=SVector{2,Float64}(0.0,1.0)
    localCoord2D(iXZ,xz,horizontalVector,normalVector)
end


function localCoord3D(ix::Integer,iy::Integer,iz::Integer,Δx::Float64,Δy::Float64,Δz::Float64; startX::Float64=-Δx,startY::Float64=-Δy,startZ::Float64=-Δz)
    # this gives simple Cartesian coordinates
    iXYZ = SVector{3,Integer}(ix,iy,iz)
    xyz = iXYZ .* SVector{3,Float64}(Δx,Δy,Δz) + SVector{3,Float64}(startX,startY,startZ)
    horizontalVector1=SVector{3,Float64}(1.0,0.0,0.0)
    horizontalVector2=SVector{3,Float64}(0.0,1.0,0.0)
    normalVector=SVector{3,Float64}(0.0,0.0,1.0)
    localCoord3D(iXYZ,xyz,horizontalVector1,horizontalVector2,normalVector)
end


# local Cartesian using normals that are defined by the point of interest

function localCoord3D(ix::Integer,iy::Integer,iz::Integer,Δx::Float64,Δy::Float64,Δz::Float64,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}; startX::Float64=-Δx,startY::Float64=-Δy,startZ::Float64=-Δz,pCentre::SVector{3,Float64}=SVector(0.0,0.0,0.0))
    iXYZ = SVector{3,Integer}(ix,iy,iz)
    xyz = iXYZ .* SVector{3,Float64}(Δx,Δy,Δz) + SVector{3,Float64}(startX,startY,startZ)
    
    # getting local vectors
    localVertical = p_local_to_ECEF(xyz,pOrigin,R)-pCentre
    if localVertical === 0.0
        localVertical = SVector(0.0,1.0,0.0)
    end
    normalGlobal=normalize(localVertical)
    normalLocal = R'*normalGlobal
    horizonY_tentative = SVector(0.0,1.0,0.0)
    horizontalVector1=normalize(cross(horizonY_tentative,normalLocal))
    horizontalVector2=normalize(cross(normalLocal,horizontalVector1))
    normalVector=normalLocal
    localCoord3D(iXYZ,xyz,horizontalVector1,horizontalVector2,normalVector)
end



function localCoord2D(ix::Integer,iz::Integer,Δx::Float64,Δz::Float64,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}; startX::Float64=-Δx,startZ::Float64=-Δz,pCentre::SVector{3,Float64}=SVector(0.0,0.0,0.0))
    
    iy = 0
    Δy = 0.0
    localCoord2D(localCoord3D(ix,iy,iz,Δx,Δy,Δz,pOrigin,R;startX=startX,startZ=startZ,pCentre=pCentre))

end





function GeoPoint(lat::Float64, lon::Float64; alt=0.0, ell=DEFAULT_ELLIPSOID[])
    lla = LLA(lat,lon, alt) # be careful LLA uses degrees by default!!
    ecef_coords = ECEF(lla,ell)
    radius = norm([ecef_coords.x,ecef_coords.y,ecef_coords.z])
    GeoPoint(lat, lon, alt, SVector(ecef_coords.x, ecef_coords.y, ecef_coords.z),radius)
end

function GeoPoint(x::Real,y::Real,z::Real; ell=DEFAULT_ELLIPSOID[])
    ecef = SVector{3,Float64}(x,y,z)
    GeoPoint(ecef;ell=ell)
end

function GeoPoint(ecef::SVector{3,Float64}; ell=DEFAULT_ELLIPSOID[])
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function +(a::GeoPoint,b::GeoPoint; ell=DEFAULT_ELLIPSOID[])
    ecef=a.ecef + b.ecef
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function -(a::GeoPoint,b::GeoPoint; ell=DEFAULT_ELLIPSOID[])
    ecef=a.ecef - b.ecef
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end

function /(a::GeoPoint,c::Real; ell=DEFAULT_ELLIPSOID[])
    ecef=a.ecef / c
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end


function *(a::GeoPoint,c::Real; ell=DEFAULT_ELLIPSOID[])
    ecef=a.ecef * c
    lla = LLA(ECEF(ecef...),ell)
    radius = norm(ecef)
    GeoPoint(lla.lat,lla.lon,lla.alt,ecef,radius)
end


function effectiveRadius(a::GeoPoint,r0::Float64; ell=DEFAULT_ELLIPSOID[])
    radiusPlanetHere = GeoPoint(a.lat,a.lon).radius 
    ratio = r0/radiusPlanetHere
    return a.radius*ratio
end

function makeLocalCoordinateUnitVectors(p1::GeoPoint,p2::GeoPoint;p0::SVector{3,Float64}, p2_1::SVector{3,Float64})

    makeLocalCoordinateUnitVectors(p1.ecef,p2.ecef;p0=p0,p2_1=p2_1)
end


function makeLocalCoordinateUnitVectors(p1::GeoPoint, p2::GeoPoint; p0::GeoPoint = (p1 + p2) / 2.0,p2_1::GeoPoint = p2 - p1)
    # Prepare ECEF equivalents first
    ecef_p1 = p1.ecef
    ecef_p2 = p2.ecef
    ecef_p0 = p0.ecef
    ecef_p2_1 = p2_1.ecef

    makeLocalCoordinateUnitVectors(ecef_p1, ecef_p2; p0 = ecef_p0, p2_1 = ecef_p2_1)
end






function makeLocalCoordinateUnitVectors(p1::SVector{3,Float64},p2::SVector{3,Float64}; p0::SVector{3,Float64}=(p1+p2)/2.0, p2_1::SVector{3,Float64}=p2-p1)

    # this function will define local (x,y,z) coordinates centred at pOrigin
    # 
    # z axis should be parallel to p0 
    # y axis is the normal to the plane that is defined by p0 and p2_1
    # x axis is the vector on the plane that is p0 and p2_1 which can be similar to p2_1 direction
    
    if norm(p2_1) === 0.0
        p2_1 = SVector(1.0,0.0,0.0)
    end

    if norm(p0) === 0.0
        p0 = SVector(0.0,0.0,1.0)
    end

    x_axis_tentative = normalize(p2_1)
    z_axis = normalize(p0)
    # y-axis: complete right-handed system
    y_axis = normalize(cross(z_axis, x_axis_tentative))
    x_axis = cross(y_axis, z_axis)  # now perfectly orthogonal

    #Rotation matrix
    R = SMatrix{3,3,Float64}(
        x_axis[1], x_axis[2], x_axis[3],
        y_axis[1], y_axis[2], y_axis[3],
        z_axis[1], z_axis[2], z_axis[3]
    )
    # x_axis = R[:,1] ... etc
    return R
end

p_ECEF_to_local(p_3D::SVector{3,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = R' * (p_3D - pOrigin)

p_local_to_ECEF(x_2D,z_2D,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(x_2D,0.e0,z_2D)

p_local_to_ECEF(vec2D::SVector{2,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(vec2D[1],0.e0,vec2D[2])

p_local_to_ECEF(x_3D,y_3D,z_3D,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(x_3D,y_3D,z_3D)

p_local_to_ECEF(vec3D::SVector{3,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(vec3D[1],vec3D[2],vec3D[3])




function constructLocalBox(p1::GeoPoint,p2::GeoPoint,Δx::Float64,Δz::Float64,altMin::Float64,altMax::Float64;leftLimit::Float64 = 0.0, rightLimit::Float64=(p2-p1).radius)

    # 2D
    Δy = 1.0 # in metre as a dummy
    奥行きMin=0.0 # y axis range
    奥行きMax=1.0
    allGridsInGeoPoints, allGridsInCartesian, effectiveRadii=constructLocalBox(p1,p2,Δx,Δy,Δz,奥行きMin,奥行きMax,altMin,altMax;leftLimit=leftLimit,rightLimit=rightLimit)
    return allGridsInGeoPoints[:,1,:], localCoord2D.(allGridsInCartesian[:,1,:]), effectiveRadii[:,1,:]
end



function constructLocalBox(p1::GeoPoint,p2::GeoPoint,Δx::Float64,Δy::Float64,Δz::Float64,奥行きMin::Float64,奥行きMax::Float64,altMin::Float64,altMax::Float64;leftLimit::Float64=0.0,rightLimit::Float64=(p2-p1).radius,centreOption="middle")

    # 3D

    R=makeLocalCoordinateUnitVectors(p1,p2) 
    pOrigin = nothing
    xOrigin = 0.0
    if centreOption === "p1"
        pOrigin = p1.ecef # p1 centred coordinates
    elseif centreOption === "middle"
        pMiddle = (p1+p2)/2.0
        pOrigin = GeoPoint(pMiddle.lat,pMiddle.lon).ecef # the middle point centred coordinates at the surface (altitude=0.0)
        xOrigin = -(p2-p1).radius/2.0
    end

    pCentre = SVector(0.0,0.0,0.0) # This is the default centre of the planet to measure the local vertical vectors

    Nx = Int64((rightLimit-leftLimit) ÷ Δx+1) 
    Ny = Int64((奥行きMax-奥行きMin) ÷ Δy +1)
    Nz = Int64((altMax-altMin) ÷ Δz + 1 ) 

    allGridsInGeoPoints=Array{GeoPoint,3}(undef,Nx,Ny,Nz)

    allGridsInCartesian=Array{localCoord3D,3}(undef,Nx,Ny,Nz)

    effectiveRadii=Array{Float64,3}(undef,Nx,Ny,Nz)


    for iXYZ in CartesianIndices(allGridsInGeoPoints)
        
        ix, iy, iz = Tuple(iXYZ)
        
        tmpLocalPoint=localCoord3D(ix,iy,iz,Δx,Δy,Δz,pOrigin,R;startX=xOrigin+leftLimit-Δx,startY=奥行きMin-Δy,startZ=altMin-Δz,pCentre=pCentre)
        tmpGeoPoint=GeoPoint(p_local_to_ECEF(tmpLocalPoint.xyz,pOrigin,R))
        
        allGridsInGeoPoints[iXYZ]=tmpGeoPoint
        allGridsInCartesian[iXYZ] = tmpLocalPoint

        effectiveRadii[iXYZ]=effectiveRadius(tmpGeoPoint,planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer*1.e3 )
    end

    return allGridsInGeoPoints, allGridsInCartesian, effectiveRadii
    
end


function contructLocalBox(Nx::Int,Ny::Int,Nz::Int,xStart::Float64,yStart::Float64,zStart::Float64,Δx::Float64,Δy::Float64,Δz::Float64;pOrigin::GeoPoint=GeoPoint(0.0,0.0),pushingSurfacetoSphere=false)
    # here it will make all necessary arrays based on a local 
    
end
