
DEFAULT_PLANET = Ref(:Earth) #set_default_planet! can change this


# Planetary parameters (semi-major axis `a` in meters, flattening `f`)

function planet_ellipsoid(name::Symbol)
    if name == :Earth # WGS84
        a= 6378137.0
        f =1.0/298.257223563
    elseif name == :SphericalEarth
        a = 6371000.0
        f = 0.0
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


function GeoPoint(lat::Float64,lon::Float64; alt=0.0, ell=DEFAULT_ELLIPSOID[])
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

function *(c::Real,a::GeoPoint; ell=DEFAULT_ELLIPSOID[])
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

function p_ECEF_to_local2D(p_3D::SVector{3,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) 
    coord3D=R' * (p_3D - pOrigin)
    coord2D = SVector(coord3D[1],coord3D[3])
    return coord2D
end

p_local_to_ECEF(x_2D,z_2D,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(x_2D,0.e0,z_2D)

p_local_to_ECEF(vec2D::SVector{2,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(vec2D[1],0.e0,vec2D[2])

p_local_to_ECEF(x_3D,y_3D,z_3D,pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(x_3D,y_3D,z_3D)

p_local_to_ECEF(vec3D::SVector{3,Float64},pOrigin::SVector{3,Float64},R::SMatrix{3,3,Float64}) = pOrigin+R*SVector(vec3D[1],vec3D[2],vec3D[3])


function constructLocalBox(arrayModel::AbstractArray{<:Any,2},altMin::Float64,altMax::Float64,leftLimit::Float64, rightLimit::Float64; p1::GeoPoint=GeoPoint(0.0,-1.0),p2::GeoPoint=GeoPoint(0.0,1,0),centreOption="middle")
    # this function builds very differently the local box with respect to the following constructLocalBox functions
    # most of the time the user will not care the p1 and p2 when talking about the very local Cartesian coordinates so the treatment is not propre for the moment
    # 2D
    
    奥行きMin=0.0 # y axis range
    奥行きMax=1.0
    arrayModel3D = zeros(typeof(arrayModel[1]),size(arrayModel)[1],1,size(arrayModel)[2])
    output=constructLocalBox(arrayModel3D,altMin,altMax,leftLimit, rightLimit,奥行きMin,奥行きMax; p1,p2,centreOption=centreOption)
    allGridsInGeoPoints = output.allGridsInGeoPoints[:,1,:]
    allGridsInCartesian = localCoord2D.(output.allGridsInCartesian[:,1,:])
    effectiveRadii = output.effectiveRadii[:,1,:]
    return (allGridsInGeoPoints=allGridsInGeoPoints, allGridsInCartesian=allGridsInCartesian, effectiveRadii=effectiveRadii, 
            Nx=output.Nx,Ny=output.Ny,Nz=output.Nz,Δx=output.Δx,Δy=output.Δy,Δz=output.Δz,pOriginECEF=output.pOriginECEF,rotationMatrix=output.rotationMatrix)
end

function constructLocalBox(arrayModel::AbstractArray{<:Any,3},altMin::Float64,altMax::Float64,leftLimit::Float64, rightLimit::Float64,奥行きMin::Float64,奥行きMax::Float64; p1::GeoPoint=GeoPoint(0.0,-1.0),p2::GeoPoint=GeoPoint(0.0,1.0),centreOption="nothing")
    # this function builds very differently the local box with respect to the following constructLocalBox functions
    # most of the time the user will not care the p1 and p2 when talking about the very local Cartesian coordinates so the treatment is not propre for the moment
    # 3D

    Nx = size(arrayModel)[1]
    Ny = size(arrayModel)[2]
    Nz = size(arrayModel)[3]

    Δx, Δy,Δz = 1.0,1.0,1.0
    if Nx != 1
        Δx = (rightLimit-leftLimit)/(Nx-1.0)
    end
    if Ny != 1
        Δy = (奥行きMax-奥行きMin)/(Ny-1.0)
    end
    if Nz != 1
        Δz = (altMax-altMin)/(Nz-1.0)
    end
    
    output=constructLocalBox(p1,p2,Δx,Δy,Δz,奥行きMin,奥行きMax,altMin,altMax;leftLimit=leftLimit,rightLimit=rightLimit,centreOption=centreOption)
    return output

end




function constructLocalBox(p0::GeoPoint,Δx::Float64,Δz::Float64,横行きMin::Float64,横行きMax::Float64,altMin::Float64,altMax::Float64;axis_angle_deg=90.0,centreOption="p0")
    # axis_angle_deg: 0=north, 90: east

    # 2D
    Δy = 1.0 # in metre as a dummy
    奥行きMin=0.0 # y axis range
    奥行きMax=1.0
    output =constructLocalBox(p0,Δx,Δy,Δz,横行きMin,横行きMax,奥行きMin,奥行きMax,altMin,altMax;axis_angle_deg=axis_angle_deg,centreOption=centreOption)
    allGridsInGeoPoints = output.allGridsInGeoPoints[:,1,:]
    allGridsInCartesian = localCoord2D.(output.allGridsInCartesian[:,1,:])
    effectiveRadii = output.effectiveRadii[:,1,:]
    return (allGridsInGeoPoints=allGridsInGeoPoints, allGridsInCartesian=allGridsInCartesian, effectiveRadii=effectiveRadii, 
            Nx=output.Nx,Ny=output.Ny,Nz=output.Nz,Δx=output.Δx,Δy=output.Δy,Δz=output.Δz,pOriginECEF=output.pOriginECEF,rotationMatrix=output.rotationMatrix)
end

function constructLocalBox(p0::GeoPoint,Δx::Float64,Δy::Float64,Δz::Float64,横行きMin::Float64,横行きMax::Float64,奥行きMin::Float64,奥行きMax::Float64,altMin::Float64,altMax::Float64;axis_angle_deg=90.0,centreOption="p0")

    # axis_angle_deg: 0=north, 90: east
    
    
    # 3D

    # define p1 and p2 from 横行き making just +/- on the east-west direction

    eps = 0.01
    u_east  = sind(axis_angle_deg)*eps
    u_north = cosd(axis_angle_deg)*eps

    axisVector = GeoPoint(p0.lat+u_north,p0.lon+u_east)-p0
    axisVector = axisVector/axisVector.radius

    p1ᶜ = p0 - axisVector
    p2ᶜ = p0 + axisVector

    axisVector = p2ᶜ - p1ᶜ
    axisVector = axisVector/axisVector.radius
    
    p1 = p0 + 横行きMin*axisVector
    p1 = GeoPoint(p1.lat,p1.lon)

    p2 = p0 + 横行きMax*axisVector
    p2 = GeoPoint(p2.lat,p2.lon)

    return constructLocalBox(p1,p2,Δx,Δy,Δz,奥行きMin,奥行きMax,altMin,altMax;leftLimit=横行きMin,rightLimit=横行きMax,centreOption=centreOption,p0=p0)
    
end



function constructLocalBox(p1::GeoPoint,p2::GeoPoint,Δx::Float64,Δz::Float64,altMin::Float64,altMax::Float64;leftLimit::Float64 = 0.0, rightLimit::Float64=(p2-p1).radius,centreOption="middle")

    # 2D
    Δy = 1.0 # in metre as a dummy
    奥行きMin=0.0 # y axis range
    奥行きMax=1.0
    output =constructLocalBox(p1,p2,Δx,Δy,Δz,奥行きMin,奥行きMax,altMin,altMax;leftLimit=leftLimit,rightLimit=rightLimit,centreOption=centreOption)
    allGridsInGeoPoints = output.allGridsInGeoPoints[:,1,:]
    allGridsInCartesian = localCoord2D.(output.allGridsInCartesian[:,1,:])
    effectiveRadii = output.effectiveRadii[:,1,:]
    return (allGridsInGeoPoints=allGridsInGeoPoints, allGridsInCartesian=allGridsInCartesian, effectiveRadii=effectiveRadii, 
            Nx=output.Nx,Ny=output.Ny,Nz=output.Nz,Δx=output.Δx,Δy=output.Δy,Δz=output.Δz,pOriginECEF=output.pOriginECEF,rotationMatrix=output.rotationMatrix)
end



function constructLocalBox(p1::GeoPoint,p2::GeoPoint,Δx::Float64,Δy::Float64,Δz::Float64,奥行きMin::Float64,奥行きMax::Float64,altMin::Float64,altMax::Float64;leftLimit::Float64=0.0,rightLimit::Float64=(p2-p1).radius,centreOption="middle",p0=nothing)
    @show centreOption

    # 3D

    R=makeLocalCoordinateUnitVectors(p1,p2) 
    pOrigin = nothing
    xOrigin = 0.0
    if centreOption == "p0"
        pOrigin = p0.ecef
    elseif centreOption == "p1"
        pOrigin = p1.ecef # p1 centred coordinates
    elseif centreOption == "middle"
        pMiddle = (p1+p2)/2.0
        pOrigin = GeoPoint(pMiddle.lat,pMiddle.lon).ecef # the middle point centred coordinates at the surface (altitude=0.0)
        xOrigin = -(p2-p1).radius/2.0
    elseif centreOption == "nothing"
        pOrigin = p1.ecef
        xOrigin = 0.0
    elseif centreOption == "centreOfPlanet"
        pOrigin = SVector(0.0,0.0,0.0)
        xOrigin = 0.0
        #xOrigin = -(p2-p1).radius/2.0
    end

    pCentre = SVector(0.0,0.0,0.0) # This is the default centre of the planet to measure the local vertical vectors

    Nx = Int64((rightLimit-leftLimit) ÷ Δx+1) 
    Ny = Int64((奥行きMax-奥行きMin) ÷ Δy +1)
    Nz = Int64((altMax-altMin) ÷ Δz + 1 ) 


    allGridsInGeoPoints = Array{GeoPoint,3}(undef, Nx, Ny, Nz)
    allGridsInCartesian = Array{localCoord3D,3}(undef, Nx, Ny, Nz)
    effectiveRadii = Array{Float64,3}(undef, Nx, Ny, Nz)

    averagePlanetRadius =
        planet1D.my1DDSMmodel.averagedPlanetRadiusInKilometer * 1.e3

    Threads.@threads for iz in 1:Nz
        for iy in 1:Ny
            for ix in 1:Nx
                iXYZ = CartesianIndex(ix, iy, iz)

                tmpLocalPoint = localCoord3D(
                    ix, iy, iz,
                    Δx, Δy, Δz,
                    pOrigin, R;
                    startX = xOrigin + leftLimit - Δx,
                    startY = 奥行きMin - Δy,
                    startZ = altMin - Δz,
                    pCentre = pCentre,
                )

                tmpGeoPoint = GeoPoint(
                    p_local_to_ECEF(tmpLocalPoint.xyz, pOrigin, R)
                )

                allGridsInGeoPoints[iXYZ] = tmpGeoPoint
                allGridsInCartesian[iXYZ] = tmpLocalPoint
                effectiveRadii[iXYZ] = effectiveRadius(tmpGeoPoint, averagePlanetRadius)
            end
        end
    end

    return (allGridsInGeoPoints=allGridsInGeoPoints, allGridsInCartesian=allGridsInCartesian, effectiveRadii=effectiveRadii, 
            Nx=Nx,Ny=Ny,Nz=Nz,Δx=Δx,Δy=Δy,Δz=Δz,pOriginECEF=pOrigin,rotationMatrix=R)
    
end

struct GeoPointSet
    name::String
    points::Vector{GeoPoint}
    color
    marker
end

function GeoPointSet(name::String, points::Vector{GeoPoint}; color=:black, marker=:circle)
    return GeoPointSet(name, points, color, marker)
end

function enuBasis(p0::GeoPoint)
    lat = deg2rad(p0.lat)
    lon = deg2rad(p0.lon)

    east = SVector(-sin(lon), cos(lon), 0.0)

    north = SVector(
        -sin(lat) * cos(lon),
        -sin(lat) * sin(lon),
         cos(lat),
    )

    up = SVector(
        cos(lat) * cos(lon),
        cos(lat) * sin(lon),
        sin(lat),
    )

    return east, north, up
end

"""
    rotation_between_geopoints(p1, p2;
        source_planet=:Earth, target_planet=:SphericalEarth)

Return the rotation that maps ECEF vectors expressed in the local east-north-up
frame at `p1` to ECEF vectors in the local east-north-up frame at `p2`.

The points are reconstructed with the requested ellipsoids because `GeoPoint`
does not store which ellipsoid was used to create it.
"""
function rotation_between_geopoints(
    p1::GeoPoint,
    p2::GeoPoint;
    source_planet::Symbol=:Earth,
    target_planet::Symbol=:Earth,
)
    source_origin = GeoPoint(
        p1.lat, p1.lon;
        alt=p1.alt,
        ell=planet_ellipsoid(source_planet),
    )
    target_origin = GeoPoint(
        p2.lat, p2.lon;
        alt=p2.alt,
        ell=planet_ellipsoid(target_planet),
    )

    source_basis = SMatrix{3,3,Float64}(hcat(enuBasis(source_origin)...))
    target_basis = SMatrix{3,3,Float64}(hcat(enuBasis(target_origin)...))

    return target_basis * source_basis'
end

function theta_phi_shift(R, p1::GeoPoint)
    _, _, up1 = enuBasis(p1)
    up2 = normalize(R * up1)

    θ1 = acosd(clamp(up1[3], -1.0, 1.0))
    φ1 = atan(up1[2], up1[1]) * 180 / π

    θ2 = acosd(clamp(up2[3], -1.0, 1.0))
    φ2 = atan(up2[2], up2[1]) * 180 / π

    Δθ = θ2 - θ1
    Δφ = mod(φ2 - φ1 + 180, 360) - 180

    # StagYY angular nodes are expressed in radians. Keep the degree-valued
    # diagnostics as well because GeoPoint latitude/longitude use degrees.
    return (
        θshift=deg2rad(Δθ),
        ϕshift=deg2rad(Δφ),
        Δθ,
        Δφ,
        θ1,
        φ1,
        θ2,
        φ2,
    )
end


function realToVirtualEarthCoords(
    points::AbstractArray{<:GeoPoint},
    p1::GeoPoint,
    p2::GeoPoint,
    ;
    source_planet::Symbol=:Earth,
    target_planet::Symbol=:SphericalEarth,
    plane::Symbol=:xz,
)
    convertedGrids, rotationAngles = transform_geopoints(
        points,
        p1,
        p2;
        source_planet=source_planet,
        target_planet=target_planet,
    )

    effectiveCoords = effective_cartesian_coordinates(
        points,
        p1,
        p2;
        source_planet=source_planet,
        target_planet=target_planet,
        plane=plane,
    )

    return (convertedGrids=convertedGrids, rotationAngles=rotationAngles, effectiveCoords=effectiveCoords)
end
"""
    transform_geopoints(points, p1, p2;
        source_planet=:Earth, target_planet=:SphericalEarth)

Rigidly move an array of `GeoPoint`s from the local frame centred at `p1` to
the local frame centred at `p2`. Local east-north-up offsets are preserved,
and the returned points use the target planet's ellipsoid. The output has the
same shape as `points`.
"""
function transform_geopoints(
    points::AbstractArray{<:GeoPoint},
    p1::GeoPoint,
    p2::GeoPoint;
    source_planet::Symbol=:Earth,
    target_planet::Symbol=:Earth,
)
    source_ellipsoid = planet_ellipsoid(source_planet)
    target_ellipsoid = planet_ellipsoid(target_planet)
    source_origin = GeoPoint(p1.lat, p1.lon; alt=p1.alt, ell=source_ellipsoid)
    target_origin = GeoPoint(p2.lat, p2.lon; alt=p2.alt, ell=target_ellipsoid)
    rotation = rotation_between_geopoints(
        source_origin,
        target_origin;
        source_planet=source_planet,
        target_planet=target_planet,
    )

    angles = theta_phi_shift(rotation,source_origin)

    return map(points) do point
        source_point = GeoPoint(
            point.lat, point.lon;
            alt=point.alt,
            ell=source_ellipsoid,
        )
        target_ecef = target_origin.ecef + rotation * (source_point.ecef - source_origin.ecef)
        GeoPoint(target_ecef; ell=target_ellipsoid)
    end,
    angles
end

"""
    effective_cartesian_coordinates(points; plane=:xz)

Return point-wise Cartesian query coordinates derived from a transformed
`GeoPoint` array. These coordinates may be non-separable after an
ellipsoid-to-sphere transformation and are intended for sampling a spherical
field; they do not replace the regular `allGridsInCartesian` coordinates.

Supported planes are `:xy`, `:xz`, and `:yz`. The returned `X` and `Y` arrays
have the same shape as `points`.
"""
function effective_cartesian_coordinates(
    points::AbstractArray{<:GeoPoint};
    plane::Symbol=:xz,
)
    components = if plane === :xy
        (1, 2)
    elseif plane === :xz
        (1, 3)
    elseif plane === :yz
        (2, 3)
    else
        throw(ArgumentError("plane must be :xy, :xz, or :yz; got $plane"))
    end

    i, j = components
    X = map(point -> point.ecef[i], points)
    Y = map(point -> point.ecef[j], points)
    return (; X, Y)
end

"""
    effective_cartesian_coordinates(points, p1, p2;
        source_planet=:Earth, target_planet=:SphericalEarth, plane=:xz)

Construct point-wise lookup coordinates by radially deforming the source
ellipsoid into the target ellipsoid and rotating the result from `p1` to `p2`.
For an ellipsoid-to-sphere conversion, the fractional radius relative to the
source surface is preserved, so the centre maps to the centre and the complete
ellipsoidal surface maps to the spherical surface.
"""
function effective_cartesian_coordinates(
    points::AbstractArray{<:GeoPoint},
    p1::GeoPoint,
    p2::GeoPoint;
    source_planet::Symbol=:Earth,
    target_planet::Symbol=:SphericalEarth,
    plane::Symbol=:xz,
)
    source_ellipsoid = planet_ellipsoid(source_planet)
    target_ellipsoid = planet_ellipsoid(target_planet)
    rotation = rotation_between_geopoints(
        p1, p2;
        source_planet=source_planet,
        target_planet=target_planet,
    )

    transformed = map(points) do point
        source_point = GeoPoint(
            point.lat, point.lon;
            alt=point.alt,
            ell=source_ellipsoid,
        )
        radius = norm(source_point.ecef)
        if iszero(radius)
            SVector(0.0, 0.0, 0.0)
        else
            direction = source_point.ecef / radius
            source_surface_radius = inv(sqrt(
                (direction[1]^2 + direction[2]^2) / source_ellipsoid.a^2 +
                direction[3]^2 / source_ellipsoid.b^2
            ))
            target_surface_radius = inv(sqrt(
                (direction[1]^2 + direction[2]^2) / target_ellipsoid.a^2 +
                direction[3]^2 / target_ellipsoid.b^2
            ))
            rotation * direction *
                (radius / source_surface_radius) * target_surface_radius
        end
    end

    components = plane === :xy ? (1, 2) :
                 plane === :xz ? (1, 3) :
                 plane === :yz ? (2, 3) :
                 throw(ArgumentError("plane must be :xy, :xz, or :yz; got $plane"))
    i, j = components
    X = map(point -> point[i], transformed)
    Y = map(point -> point[j], transformed)
    return (; X, Y)
end
function localENU(p::GeoPoint, p0::GeoPoint)
    east, north, up = enuBasis(p0)
    d = p.ecef - p0.ecef

    xEast  = dot(d, east)
    yNorth = dot(d, north)
    zUp    = dot(d, up)

    return xEast, yNorth, zUp
end

function GeoPoints_to_local(points::AbstractVector{GeoPoint}, boxGrids3D)
    pOrigin = boxGrids3D.pOriginECEF 
    R = boxGrids3D.rotationMatrix

    return [
        p_ECEF_to_local(p.ecef, pOrigin, R)
        for p in points
    ]
end

"""
    PathProfile2D

A field profile sampled along one 2D source-to-detector ray.

`values` contains the field average in each segment, `baselines` contains the
corresponding segment lengths in kilometres, and all coordinates are in
metres. The ordering is always source → detector.
"""
struct PathProfile2D
    cosθ::Float64
    values::Vector{Float64}
    baselines::Vector{Float64}
    source::SVector{2,Float64}
    detector::SVector{2,Float64}
    segment_midpoints::Vector{SVector{2,Float64}}
end

function _positive_distance_to_box(
    point::SVector{2,Float64},
    direction::SVector{2,Float64},
    xlimits,
    ylimits,
)
    distances = Float64[]
    for (coordinate, component, limits) in zip(point, direction, (xlimits, ylimits))
        if component > 0
            push!(distances, (limits[2] - coordinate) / component)
        elseif component < 0
            push!(distances, (limits[1] - coordinate) / component)
        end
    end
    positive = filter(>(0.0), distances)
    isempty(positive) &&
        throw(ArgumentError("ray from detector does not enter the coordinate box"))
    return minimum(positive)
end

function _surface_distance(
    inside_interpolator,
    detector::SVector{2,Float64},
    direction::SVector{2,Float64},
    maximum_distance::Float64,
    search_step::Float64,
    tolerance::Float64,
)
    inside_interpolator(detector...) >= 0.5 ||
        throw(ArgumentError("detector must be inside material_mask"))

    lower = 0.0
    upper = min(search_step, maximum_distance)
    while upper < maximum_distance &&
          inside_interpolator((detector + upper * direction)...) >= 0.5
        lower = upper
        upper = min(upper + search_step, maximum_distance)
    end

    inside_interpolator((detector + upper * direction)...) < 0.5 ||
        throw(ArgumentError(
            "ray did not leave material_mask before reaching the interpolation grid boundary",
        ))

    while upper - lower > tolerance
        middle = (lower + upper) / 2
        if inside_interpolator((detector + middle * direction)...) >= 0.5
            lower = middle
        else
            upper = middle
        end
    end
    return (lower + upper) / 2
end

"""
    creationPaths(field, x, z, detectorLocal, cosθgrid;
        Δbaseline, material_mask, modifiedLongitude2Disk=0.0,
        center=(0.0, 0.0), boundary_search_step=nothing,
        boundary_tolerance=1.0, path_constructor=nothing)

Sample a regular 2D Cartesian field along incoming neutrino rays. `x`, `z`,
`detectorLocal`, `center`, and `Δbaseline` use metres. Segment baselines are
returned in kilometres for neutrino-oscillation routines.

The ray for `cosθ=-1` crosses the diameter. `cosθ=0` follows the
counter-clockwise tangent used by the original `creationPaths`. A detector
constructed with longitude 180° is the mirror of longitude 0° and therefore
selects the opposite half-disk.

The topographic source is found from the first transition out of
`material_mask` along the detector-to-source ray. Field averages use
three-point Gauss-Legendre quadrature. Profiles are stored in the physically
important source → detector order.

When `path_constructor` is supplied, it is called as
`path_constructor(profile.values, profile.baselines)` for every angle. This
lets callers construct `Neurthino.Path` objects without coupling GeoPoints to
the neutrino module.
"""
function creationPaths(
    field::AbstractMatrix,
    x_coordinates::AbstractVector,
    z_coordinates::AbstractVector,
    detectorLocal,
    cosθgrid;
    Δbaseline::Real,
    material_mask::AbstractMatrix,
    modifiedLongitude2Disk::Real=0.0,
    center=(0.0, 0.0),
    boundary_search_step=nothing,
    boundary_tolerance::Real=1.0,
    path_constructor=nothing,
)
    x = Float64.(collect(x_coordinates))
    z = Float64.(collect(z_coordinates))
    size(field) == (length(x), length(z)) ||
        throw(DimensionMismatch(
            "field size $(size(field)) must equal (length(x), length(z)) = " *
            "$((length(x), length(z)))",
        ))
    size(material_mask) == size(field) ||
        throw(DimensionMismatch("material_mask must have the same size as field"))
    all(diff(x) .> 0) || throw(ArgumentError("x coordinates must be strictly increasing"))
    all(diff(z) .> 0) || throw(ArgumentError("z coordinates must be strictly increasing"))
    Δbaseline > 0 || throw(ArgumentError("Δbaseline must be positive"))
    boundary_tolerance > 0 ||
        throw(ArgumentError("boundary_tolerance must be positive"))

    longitude = mod(Float64(modifiedLongitude2Disk), 360.0)
    (isapprox(longitude, 0.0; atol=1e-8) ||
     isapprox(longitude, 180.0; atol=1e-8)) ||
        throw(ArgumentError("modifiedLongitude2Disk must be 0° or 180°"))

    detector = SVector{2,Float64}(detectorLocal)
    origin = SVector{2,Float64}(center)
    radial = detector - origin
    norm(radial) > 0 || throw(ArgumentError("detectorLocal cannot equal center"))

    expected_x_sign = isapprox(longitude, 0.0; atol=1e-8) ? 1 : -1
    radial_x_tolerance = sqrt(eps(Float64)) * norm(radial)
    if abs(radial[1]) > radial_x_tolerance &&
       sign(radial[1]) != expected_x_sign
        throw(ArgumentError(
            "detectorLocal is inconsistent with modifiedLongitude2Disk=" *
            "$(modifiedLongitude2Disk)°",
        ))
    end

    radial_unit = radial / norm(radial)
    tangent_unit = SVector(-radial_unit[2], radial_unit[1])
    cosines = Float64.(collect(cosθgrid))
    all(-1.0 .<= cosines .<= 0.0) ||
        throw(ArgumentError("cosθgrid values must lie in [-1, 0]"))

    field_interpolator = LinearInterpolation(
        (x, z),
        Float64.(field);
        extrapolation_bc=Flat(),
    )
    inside_interpolator = LinearInterpolation(
        (x, z),
        Float64.(material_mask);
        extrapolation_bc=0.0,
    )

    grid_step = min(minimum(diff(x)), minimum(diff(z)))
    search_step = isnothing(boundary_search_step) ?
                  grid_step / 2 :
                  Float64(boundary_search_step)
    search_step > 0 ||
        throw(ArgumentError("boundary_search_step must be positive"))

    profiles = Vector{PathProfile2D}(undef, length(cosines))
    Threads.@threads for angle_index in eachindex(cosines)
        cosine = cosines[angle_index]
        direction = cosine * radial_unit +
                    sqrt(max(0.0, 1 - cosine^2)) * tangent_unit
        maximum_distance = _positive_distance_to_box(
            detector,
            direction,
            extrema(x),
            extrema(z),
        )
        path_length = _surface_distance(
            inside_interpolator,
            detector,
            direction,
            maximum_distance,
            search_step,
            Float64(boundary_tolerance),
        )
        source = detector + path_length * direction

        number_segments = max(1, ceil(Int, path_length / Δbaseline))
        edges = collect(range(0.0, path_length; length=number_segments + 1))
        baselines = diff(edges) .* 1e-3
        values = Vector{Float64}(undef, number_segments)
        midpoints = Vector{SVector{2,Float64}}(undef, number_segments)

        # Three-point Gauss-Legendre rule, normalized for a segment average.
        quadrature_nodes = (-sqrt(3 / 5), 0.0, sqrt(3 / 5))
        quadrature_weights = (5 / 18, 4 / 9, 5 / 18)
        source_to_detector = -direction
        for segment_index in eachindex(values)
            left = edges[segment_index]
            right = edges[segment_index + 1]
            middle = (left + right) / 2
            half_width = (right - left) / 2
            midpoints[segment_index] = source + middle * source_to_detector
            values[segment_index] = sum(
                weight * field_interpolator(
                    (source + (middle + half_width * node) * source_to_detector)...,
                )
                for (node, weight) in zip(quadrature_nodes, quadrature_weights)
            )
        end

        profiles[angle_index] = PathProfile2D(
            cosine,
            values,
            baselines,
            source,
            detector,
            midpoints,
        )
    end

    paths = isnothing(path_constructor) ?
            nothing :
            [path_constructor(profile.values, profile.baselines) for profile in profiles]
    return (;
        cosθgrid=cosines,
        profiles,
        paths,
        detector,
        center=origin,
        modifiedLongitude2Disk=longitude,
    )
end

function creationPaths(
    field::AbstractMatrix,
    x_coordinates::AbstractVector,
    z_coordinates::AbstractVector,
    detectorLocal,
    binning::NamedTuple;
    kwargs...,
)
    hasproperty(binning, :cosθgrid) ||
        throw(ArgumentError("binning must contain cosθgrid"))
    return creationPaths(
        field,
        x_coordinates,
        z_coordinates,
        detectorLocal,
        binning.cosθgrid;
        kwargs...,
    )
end
