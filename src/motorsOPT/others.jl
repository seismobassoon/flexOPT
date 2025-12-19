
timeDimensionString="t" 
# if the user does not want to use "t" for the time marching scheme, it should be changed
# and this "t" should be the last element in coordinates

function findCartesianDependency(expression,Ndimension,∂)
    expressionDependency=ones(Int,Ndimension)
    for iDimension in 1:Ndimension
        if typeof(expand_derivatives(∂[iDimension](expression))==0) <: Bool 
            if expand_derivatives(∂[iDimension](expression))==0
                expressionDependency[iDimension] = 0
            end
        end
    end
    return expressionDependency
end

function makeMixPartials(orders,coordinates;field=identity)
    # this function will give a matrix of mixed partial derivative operators 
    # coordinates should be an array of Symbolics variables 
    # orders is a matrix 


    Ndimension=length(coordinates)
    if length(orders)!== Ndimension
        @error "the highest orders array has not the same dimension as that of the coordinates"
    end

    ∂ = []
    for iDim in 1:Ndimension
        ∂ = push!(∂,Differential(coordinates[iDim]))
    end

    ∇ = Array{Any,Ndimension}(undef, Tuple(orders))
    R=CartesianIndices(∇)
        
    ∇ .= field
    for I in R
        for iDim in 1:Ndimension
            ∇[I] = (∂[iDim]^(I[iDim]-1))(∇[I])
        end
    end
   
    return ∇
end

function varMmaker(maxPointsUsed,coordinates,vars) 
    # this will make an array of material coeffs for with a local Cartesian grid (max points used for a node)
    Ndimension = length(coordinates)

    R = CartesianIndices(Tuple(maxPointsUsed))

    varM=Array{Any,2}(undef,length(vars),length(R))
   
    for iVar in eachindex(vars)


        newstring=split(string(vars[iVar]),"(")[1]
     
        
        CartesianDependency=findCartesianDependency(vars[iVar],Ndimension)
       
        smallVarM=Symbolics.variables(Symbol(newstring),1:length(R))
        for j in R
            linearJ=LinearIndices(R)[j]
            realJ=(car2vec(j).-1).*CartesianDependency .+1 # if there is no dependence on a direction, it should get the same name
            linearRealJ=LinearIndices(R)[CartesianIndex(realJ...)]
            smallVarM[linearJ]=smallVarM[linearRealJ]
        end
        varM[iVar,:]=smallVarM
    end
    
    return varM
end

function PDECoefFinder(orders,coordinates,expr,field,vars)
    # PDECoefFinder cannot detect the material partials × material partials for the moment!! 
    # I know how to do it, but eq. 40 should be then more generalised (kind of the product of partials of different materials)

    # maxPolynomialOrderMaterial is also a chelou thing, that I need to work on more systematically
    # like the powers of partials should also be included but here search for Rm[1], yeah, that's what I am doing


    Ndimension = length(coordinates)
    alpha=[]
    
    maxPolynomialOrderMaterial = 2*(maximum(orders)-1)
    ∇=makeMixPartials(orders,coordinates;field=field)
    R=CartesianIndices(∇)
    expr=mySimplify(expr)


    for i in R
        term_searched = ∇[i]

        tmpCoeff = myCoeff(expr,term_searched)
        if tmpCoeff !== 0
            isTmpCoeffAConstant=true
            for iVar in eachindex(vars)
                
                ∇m=makeMixPartials(orders,coordinates;field=vars[iVar]) # material partials
                Rm=CartesianIndices(∇m)
                for j in Rm
                    term_material_searched = ∇m[j]
                    tmpCoeffMaterial = myCoeff(tmpCoeff,term_material_searched)
                    
                    if tmpCoeffMaterial !==0
                        isTmpCoeffAConstant=false
                        isOKtoinclude =true
                        # This is to avoid partials of other material 
                        for jVar in eachindex(vars)
                            if jVar !== iVar
                                ∇n = makeMixPartials(orders,coordinates;field=vars[jVar]) 
                                for jj in Rm
                                    if jj !== Rm[1]
                                        term_material_searched_plus = ∇n[jj]
                                        differentialCoeff = myCoeff(tmpCoeffMaterial,term_material_searched_plus)
                                        if differentialCoeff !==0
                                            isOKtoinclude = false
                                        end
                                    end
                                end
                            end
                        end
                        if isOKtoinclude
                            #tmpCoeffMaterial=substitute(tmpCoeffMaterial,mapping)
                            specificMaterialTerm=tmpCoeffMaterial*vars[iVar]
                            tmpAlphaIJ = (node=specificMaterialTerm,nᶜ=j,n=i) # the famous n prime and n in the equation 56 or 40
                            alpha = push!(alpha,tmpAlphaIJ)
                        end
                    end
                end
                for matPower in 2:maxPolynomialOrderMaterial
                    tmpCoeffMaterial = myCoeff(tmpCoeff,vars[iVar]^matPower)
                    #tmpCoeffMaterial=substitute(tmpCoeffMaterial,mapping)
                    if tmpCoeffMaterial !==0
                        specificMaterialTerm=tmpCoeffMaterial*vars[iVar]^matPower
                        tmpAlphaIJ = (node=specificMaterialTerm,nᶜ=Rm[1],n=i) # the famous n prime and n in the equation 56 or 40
                        alpha = push!(alpha,tmpAlphaIJ)
                    end
                end
            end
            if isTmpCoeffAConstant
                specificMaterialTerm = tmpCoeff
                tmpAlphaIJ=(node=specificMaterialTerm,nᶜ=R[1],n=i)
                alpha = push!(alpha,tmpAlphaIJ)
            end
        end
    end
    alpha=unique(alpha)

    return alpha # varM: iVar and linearised cartesian indices
end 

function TaylorCoefInversion(coefInversionDict::Dict)

    # this function will give Float64 array (not Any)

    # the user might want to have a look on illposedTaylorCoefficientsInversion_legend, which is deprecated as of 10/06/2025.


    # based on the equation 27 (of the version 10/06/2025 FD2025 : \psi_{;\mu,\nu}^{(l)}[\nu+\mu]=\sum_\eta C_{\mu+\eta;\mu,\nu}^{(l)} \psi[\nu+\mu+\eta]), we need to perform this inversion anyways for all the point \mu inside L(\nu) (the concerned points for \nu)


    # be careful that pointsIndices is now a 1D array of integer vectors!!

    @unpack coordinates,multiOrdersIndices,pointsIndices, Δ, WorderBspline, modifiedμ = coefInversionDict

    Ndimension=length(coordinates)

    if length(Δ) !== Ndimension
        @error "the numerical delta increment has not the same dimension!"
    end

    #@show vcat(deep_flatten(pointsIndices))
    #@show pointsIndices = vec2car(vcat(deep_flatten(pointsIndices)),Ndimension)
    #pointsIndices=to_cartesian_list(pointsIndices,Ndimension)
    pointsIndices = vec(CartesianIndex.(Tuple.(pointsIndices)))
    μpointsIndices = pointsIndices # which can be changed 
    
    numberOfEtas = length(pointsIndices)
    numberOfLs   = length(multiOrdersIndices)

    numberOfMus = length(μpointsIndices) # which can be \mathbb{Z}/2 or something else ...

    CˡηGlobal = Array{Float64,3}(undef,numberOfEtas,numberOfLs,numberOfMus)

    # this is the C^{(l)}_{\mu+\eta; μ, \nu}

    for μ_oneD in axes(CˡηGlobal,3)
        #@show typeof(multiOrdersIndices),typeof(pointsIndices)
        CˡηGlobal[:,:,μ_oneD]=TaylorCoefInversion(numberOfLs,numberOfEtas,multiOrdersIndices,pointsIndices,μpointsIndices,Δ,μ_oneD,WorderBspline,modifiedμ)
    end 

    return @strdict(CˡηGlobal)

end

function TaylorCoefInversion(numberOfLs,numberOfEtas,multiOrdersIndices,pointsIndices,μpointsIndices,Δ,μ_oneD,WorderBspline,modifiedμ)


    # the old version is : illposedTaylorCoefficientsInversionSingleCentre

    # in fact, available points depend on the position of μ (=k here), we need to 'mute' some points
    # with Y_μ
    
    # for this pointsIndices are filtered for every μ
    
    Ndimension = length(WorderBspline)

    tmpPointsIndices = []
    linearIndicesUsed = []

    #modifiedμ_vector = Array{Float64,1}(undef,Ndimension)

  

    @show modifiedμ_vector = Float64.(car2vec(μpointsIndices[μ_oneD]))


   
    #for η_μ in pointsIndices
    for i in eachindex(pointsIndices)
        η_μ = pointsIndices[i]
        iSayWeSayGo = 1
        for iCoord in 1:Ndimension # Ndimension
            #@show 1,μ_oneD,μpointsIndices[μ_oneD][iCoord],WorderBspline[iCoord]+1, iCoord
            #tmp1=Num2Float64(modifiedμ[iCoord][1,μ,WorderBspline[iCoord]+1])
            #tmp2=Num2Float64(modifiedμ[iCoord][2,μ,WorderBspline[iCoord]+1])
            #tmp3=Num2Float64(modifiedμ[iCoord][3,μ,WorderBspline[iCoord]+1])
            
            tmp1=Num2Float64(safeget(modifiedμ[iCoord],1,μpointsIndices[μ_oneD][iCoord],WorderBspline[iCoord]+1))
            tmp2=Num2Float64(safeget(modifiedμ[iCoord],2,μpointsIndices[μ_oneD][iCoord],WorderBspline[iCoord]+1))
            tmp3=Num2Float64(safeget(modifiedμ[iCoord],3,μpointsIndices[μ_oneD][iCoord],WorderBspline[iCoord]+1))
            #@show tmp1, tmp2, tmp3

            #modifiedμ_vector[iCoord] = tmp3
            if WorderBspline[iCoord] === -1 # this will use Y everywhere (for ν+μ = ν)
                iSayWeSayGo *= 1
            elseif  tmp1 <=η_μ[iCoord] <= tmp2
                iSayWeSayGo *= 1
            else
                iSayWeSayGo *= 0
            end
        end
        if iSayWeSayGo === 1
            tmpPointsIndices=push!(tmpPointsIndices,η_μ)
            linearIndicesUsed=push!(linearIndicesUsed,i)
        end
    end

    tmpNumberOfEtas = length(tmpPointsIndices)
    tmpTaylorExpansionCoeffs = Array{Any,2}(undef,numberOfLs,tmpNumberOfEtas)


    for iAvailablePoint in eachindex(tmpPointsIndices)
        η_μ = tmpPointsIndices[iAvailablePoint]
        #η = tmpPointsIndices[i]-pointsIndices[μ]
        #η = tmpPointsIndices[i] - modifiedμ_vector
        η = Float64.(car2vec(η_μ)) .- modifiedμ_vector 
        distances= η .* Δ
        for j in multiOrdersIndices
            linearJ = LinearIndices(multiOrdersIndices)[j]
            orders = car2vec(j).-1
            numerator = prod(distances .^orders)
            denominator=prod(factorial.(orders))
            tmpTaylorCoeffs = numerator/denominator
            tmpTaylorExpansionCoeffs[linearJ,iAvailablePoint]=tmpTaylorCoeffs 

        end
    end

    # here we do the famous inversion (ttttttt) even though this code is essentially a forward problem
    
    aa=transpose(tmpTaylorExpansionCoeffs)*tmpTaylorExpansionCoeffs
    aa=Num2Float64.(aa)
    @show typeof(aa),aa,size(aa)
    invaa= myInv(aa)
    tmpCˡηlocal=invaa*transpose(tmpTaylorExpansionCoeffs)


    Cˡηlocal = Array{Any,2}(undef,numberOfEtas,numberOfLs)

    Cˡηlocal .= 0

    for j in eachindex(tmpPointsIndices)
        Cˡηlocal[linearIndicesUsed[j],:] = tmpCˡηlocal[j,:]
    end

    return Cˡηlocal
end


#
#
#
#
#
# below are the obsolete functions

function illposedTaylorCoefficientsInversion(coordinates,multiOrdersIndices,multiPointsIndices;testOnlyCentre=true,Δ=nothing,timeMarching=false)

    # this function is deprecated, no more testOnlyCentre/timeMarching options are allowed (nor Symbolics Δ)

    # here we propose a big Taylor expansion matrix with Δcoordinates, symbolically (when Δ=nothing or Δ as a symbolic array)
    #   or numerically otherwise


    Ndimension=length(coordinates)

    if Δ === nothing
       Δ = Symbolics.variables(:Δ,1:Ndimension)
    else
        if length(Δ) !== Ndimension
            @error "the numerical delta increment has not the same dimension!"
        end
    end

    numberOfEtas = length(multiPointsIndices)
    numberOfLs   = length(multiOrdersIndices)

    #TaylorExpansionCoeffs=Array{Any,Ndimension*2}(undef, multiOrdersIndices[end],multiPointsIndices[end])
    #TaylorExpansionCoeffs=Array{Any,Ndimension*2}(undef, Tuple(vcat(collect(Tuple(multiOrdersIndices[end])),collect(Tuple(multiPointsIndices[end])))))
    
    CˡηGlobal = Array{Any,3}(undef,numberOfEtas,numberOfLs,numberOfEtas)
   

    tmpVecForMiddlePoint = (car2vec(multiPointsIndices[end]).-1 ).÷2 .+1 # only valid for testOnlyCentre
    midTimeCoord=nothing
    if timeMarching
        midTimeCoord=car2vec(multiPointsIndices[end])[end]-1
        tmpVecForMiddlePoint[end]=midTimeCoord
    end
    midK=vec2car(tmpVecForMiddlePoint)

    for k in multiPointsIndices
        linearK = LinearIndices(multiPointsIndices)[k]
        
        if !testOnlyCentre || k === midK || (timeMarching && car2vec(k)[end] === midTimeCoord && !testOnlyCentre) # because we cannot predict more than one futures
            CˡηGlobal[:,:,linearK]=illposedTaylorCoefficientsInversionSingleCentre(numberOfLs,numberOfEtas,multiOrdersIndices,multiPointsIndices,Δ,k)
        end
    end 

    if testOnlyCentre
        CˡηCentre = CˡηGlobal[:,:,LinearIndices(multiPointsIndices)[midK]]
        CˡηGlobal = nothing
        return CˡηCentre,Δ,multiOrdersIndices
    else
        #@show CˡηGlobal
        return CˡηGlobal,Δ,multiOrdersIndices
    end
end

function illposedTaylorCoefficientsInversionSingleCentre(numberOfLs,numberOfEtas,multiOrdersIndices,multiPointsIndices,Δ,k)

    # this function is deprecated as of 10/06/2025

    # this should be completely numerical
    TaylorExpansionCoeffs = Array{Num,2}(undef,numberOfLs,numberOfEtas)
    for i in multiPointsIndices
        linearI = LinearIndices(multiPointsIndices)[i]
        η = car2vec(i-k)
        distances= η .* Δ
        for j in multiOrdersIndices
            linearJ = LinearIndices(multiOrdersIndices)[j]
            orders = car2vec(j).-1
            numerator = prod(distances .^orders)
            denominator=prod(factorial.(orders))
            tmpTaylorCoeffs = numerator/denominator
            TaylorExpansionCoeffs[linearJ,linearI]=tmpTaylorCoeffs 

        end
    end
    # here we do the famous inversion (ttttttt) even though this code is essentially a forward problem
    
    aa=transpose(TaylorExpansionCoeffs)*TaylorExpansionCoeffs
    invaa= myInv(aa)
    Cˡηlocal=invaa*transpose(TaylorExpansionCoeffs)
    return Cˡηlocal
end



function spaceCoordinatesConversionfunctions(absorbingBoundaries, NdimensionMinusTime)
    offset_model = vec2car(absorbingBoundaries[1, 1:NdimensionMinusTime])
    #offset_empty = vec2car(spacePointsUsed)

    model2whole(a::CartesianIndex) = a + offset_model
    whole2model(a::CartesianIndex) = a - offset_model
    #whole2empty(a::CartesianIndex) = a + offset_empty
    #empty2whole(a::CartesianIndex) = a - offset_empty
    #model2empty(a::CartesianIndex) = whole2empty(model2whole(a))
    #empty2model(a::CartesianIndex) = whole2model(empty2whole(a))
    return(; model2whole, whole2model)
    #return (; model2whole, whole2model, whole2empty, empty2whole, model2empty, empty2model)
end

function BouncingCoordinates(a::CartesianIndex,PointsUsed)
    #
    # this will bounce the boundary inside the PointsUsed vector
    #
    # i.e. get the nearby coordinates inside the domain to fake the continuity

    if length(a) !== length(PointsUsed)
        @error "cannot bound this CartesianIndex due to the dimension mismatch"
    end
    avector=car2vec(a)
    for iCoord in eachindex(avector)
        if avector[iCoord] < 1
            avector[iCoord] = 1
        elseif avector[iCoord] > PointsUsed[iCoord]
            avector[iCoord] = PointsUsed[iCoord]
        end
    end
    a=vec2car(avector)
    return a
end