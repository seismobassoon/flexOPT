# Analytical computation of WYYKK integral (look at demo1DKernelIntegral.ipynb)

using .commonBatchs, UnPack, Symbolics


function WYYKKIntegralNumerical(params;ImakeReport=true)
    
    @unpack orderBspline1D, YorderBspline1Dμᶜ, YorderBspline1Dμ, μᶜs, μs, maxNode, ν, lᶜ_nᶜ_max, l_n_max, Δ = params

    paramsForSymbolic = @strdict orderBspline1D YorderBspline1Dμᶜ YorderBspline1Dμ μᶜs μs maxNode ν lᶜ_nᶜ_max l_n_max ImakeReport

    output = myProduceOrLoad(WYYKKIntegralPureSymbolic,paramsForSymbolic,"WYYKKIntegralSymbolic")

    WYYKK_integral = output["WYYKK_integral"]

    coefWYYKK = Array{Float64, 5}(undef,l_n_max+1,lᶜ_nᶜ_max+1,length(μs),length(μᶜs),length(ν))

    nodes = WYYKK_integral.nodes
    numericalNodes = Δ .* nodes
    numNodes = WYYKK_integral.numberNodes
    x = WYYKK_integral.variables[1]
    Δx = WYYKK_integral.variables[2]

    for iν ∈ eachindex(ν), iμᶜ ∈ eachindex(μᶜs), iμ ∈ eachindex(μs), lᶜ_nᶜ ∈ 0:lᶜ_nᶜ_max, l_n ∈ 0:l_n_max
        l_n_slot=l_n+1
        lᶜ_nᶜ_slot = lᶜ_nᶜ+1
        tmpAntiDerivative=WYYKK_integral.data[:,1,l_n_slot,lᶜ_nᶜ_slot,iμ,iμᶜ,iν]
        tmpCoef = 0.0
        for ι in 1:numNodes-1
            xLeft = numericalNodes[ι]
            xRight = numericalNodes[ι+1]
            expr = tmpAntiDerivative[ι]
            rightValue = Symbolics.value(Symbolics.substitute(expr, Dict(x => xRight, Δx => Δ)))
            leftValue = Symbolics.value(Symbolics.substitute(expr, Dict(x => xLeft, Δx => Δ)))
            tmpCoef += rightValue-leftValue
        end
        coefWYYKK[l_n_slot,lᶜ_nᶜ_slot,iμ,iμᶜ,iν] = tmpCoef
    end
    return coefWYYKK
end

function WYYKKIntegralPureSymbolic(params::Dict)
    # Δ should be strictly Float64

    # orders: -1 -> indicator function, 0 -> box car, >=1 -> B-spline

    # this computes the analytical value of the 1D integral between B-spline fns and weighted Taylor kernels
    # \int dx Bspline Y_μᶜ Y_μ  K_{lᶜ-nᶜ}(y-y_μᶜ) K_{l-n}(y-y_μ)

    # unlike the previous integralBsplineTaylorKernels1D, it computes for a specific ν
    # Cˡη;μ are computed for a specific geometry, so even though the boundary condition reduces
    # the number of available points, each Taylor expansion for K_{l-n}(y-y_μ) should be Ok

    @unpack orderBspline1D,YorderBspline1Dμᶜ,YorderBspline1Dμ,μᶜs,μs,maxNode,ν,lᶜ_nᶜ_max,l_n_max,ImakeReport = params

    nodesFromOne = collect(1:1:maxNode) # ∈ Z like [1,2,3], an array of integers collect(1:1:N) (nothing else!!)

    allNodes = unique(sort(vcat(
        Float64.(nodesFromOne),
        Float64.(ν),
        Float64.(μs),
        Float64.(μᶜs),
    )))

    to_indices(xs, master) = searchsortedfirst.(Ref(master), Float64.(xs))

    idx_nodesFromOne = to_indices(nodesFromOne, allNodes)
    idx_ν            = to_indices(ν, allNodes)
    idx_μs           = to_indices(μs, allNodes)
    idx_μᶜs          = to_indices(μᶜs, allNodes)


    # for B-spline

    paramsBSν  = (maximumOrder=orderBspline1D, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_refPoints = idx_nodesFromOne, idx_selectedPoints = idx_ν)
    paramsBSμᶜ = (maximumOrder=YorderBspline1Dμᶜ, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_refPoints = idx_μs, idx_selectedPoints = idx_μᶜs)
    paramsBSμ  = (maximumOrder=YorderBspline1Dμ, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_refPoints = idx_μs, idx_selectedPoints = idx_μs)
    # idx_nodesNum : an ordinary consecutive integer increment from 1 (the numerical nodes with Δy)
    # idx_refPoints_original : supporting nodes to construct the Bspline family
    # idx_selectedPoints  : the node addresses that user needs to take, should be a subset of idx_refPoints_original


    # for Taylor Expansions

    paramsTEμᶜ = (maxL_MINUS_N=lᶜ_nᶜ_max, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_selectedPoints = idx_μᶜs)
    paramsTEμ  = (maxL_MINUS_N=l_n_max, allNodes = allNodes, idx_nodesNum = idx_nodesFromOne, idx_selectedPoints = idx_μs)

    # Computing Wν,Yμᶜ,Yμ 

    Wν = constructBsplineFamily(paramsBSν);
    Yμᶜ = constructBsplineFamily(paramsBSμᶜ);
    Yμ = constructBsplineFamily(paramsBSμ);

    # Computing Kμᶜ, Kμ

    Kμᶜ=constructTaylorExpansions(paramsTEμᶜ);
    Kμ=constructTaylorExpansions(paramsTEμ);

    # Computing WYYKK

    WYYKK=CompactSymbolicFunctions(Yμ.b.nodes,1;auxDims=(l_n_max+1, lᶜ_nᶜ_max+1, length(μs),length(μᶜs),length(ν)),variables=Yμ.b.variables)

    # note that W and Ys are computed on the fly without DrWatson (for the moment) 
    # therefore, boxcar functions should be properly called
    # however, we need to be careful with further developments around this 

    YorderSlot = maximum((YorderBspline1D + 1,1)) # for Yorder = -1 : just boxcar
    orderSlot = maximum((orderBspline1D + 1,1)) # for order = -1 : just boxcar
    derivSlot = 1 # no derivatives

    for iν ∈ eachindex(ν), iμᶜ ∈ eachindex(μᶜs), iμ ∈ eachindex(μs), lᶜ_nᶜ ∈ 0:lᶜ_nᶜ_max, l_n ∈ 0:l_n_max
        l_n_slot=l_n+1
        lᶜ_nᶜ_slot = lᶜ_nᶜ+1
        WYYKK.data[:,1,l_n_slot, lᶜ_nᶜ_slot, iμ, iμᶜ,iν] = mySimplify(Wν.b.data[:,iν,derivSlot,orderSlot].*Yμᶜ.b.data[:,iμᶜ,derivSlot,YorderSlot].*Yμ.b.data[:,iμ,derivSlot,YorderSlot].*Kμᶜ.k.data[:,iμᶜ,lᶜ_nᶜ_slot].*Kμ.k.data[:,iμ,l_n_slot])
    end

    WYYKK_integral = integrate(WYYKK,x)

    reportDir = nothing

    if ImakeReport
        
        reportDir = joinpath(pwd(), "WYYKKIntegralPureSymbolic_report")
        mkpath(reportDir)

        # Page 1: B-spline families
        save_bspline_report_plot(
            joinpath(reportDir, "01_Wnu.pdf"),
            Wν.b;
            derivOrder=0,
            order=max(orderBspline1D, 0),
            N=80,
            title="Wν",
        )

        save_bspline_report_plot(
            joinpath(reportDir, "02_Ymuc.pdf"),
            Yμᶜ.b;
            derivOrder=0,
            order=max(YorderBspline1Dμᶜ, 0),
            N=80,
            title="Yμᶜ",
        )

        save_bspline_report_plot(
            joinpath(reportDir, "03_Ymu.pdf"),
            Yμ.b;
            derivOrder=0,
            order=max(YorderBspline1Dμ, 0),
            N=80,
            title="Yμ",
        )

        # Kμᶜ pages
        for iμᶜ in eachindex(μᶜs), lᶜ_nᶜ in 0:lᶜ_nᶜ_max
            lᶜ_nᶜ_slot = lᶜ_nᶜ + 1
            save_report_plot(
                joinpath(reportDir, "10_Kmuc_i$(iμᶜ)_l$(lᶜ_nᶜ).pdf"),
                Kμᶜ.k,
                (lᶜ_nᶜ_slot,);
                N=80,
                title="Kμᶜ, iμᶜ=$iμᶜ, lᶜ-nᶜ=$lᶜ_nᶜ",
                ylabel="Taylor kernel",
            )
        end

        # Kμ pages
        for iμ in eachindex(μs), l_n in 0:l_n_max
            l_n_slot = l_n + 1
            save_report_plot(
                joinpath(reportDir, "20_Kmu_i$(iμ)_l$(l_n).pdf"),
                Kμ.k,
                (l_n_slot,);
                N=80,
                title="Kμ, iμ=$iμ, l-n=$l_n",
                ylabel="Taylor kernel",
            )
        end

        # WYYKK pages
        for iν in eachindex(ν), iμᶜ in eachindex(μᶜs), iμ in eachindex(μs), lᶜ_nᶜ in 0:lᶜ_nᶜ_max, l_n in 0:l_n_max
            l_n_slot = l_n + 1
            lᶜ_nᶜ_slot = lᶜ_nᶜ + 1

            save_report_plot(
                joinpath(reportDir, "30_WYYKK_inu$(iν)_imuc$(iμᶜ)_imu$(iμ)_lc$(lᶜ_nᶜ)_l$(l_n).pdf"),
                WYYKK,
                (l_n_slot, lᶜ_nᶜ_slot, iμ, iμᶜ, iν);
                N=80,
                title="WYYKK, iν=$iν, iμᶜ=$iμᶜ, iμ=$iμ, lᶜ-nᶜ=$lᶜ_nᶜ, l-n=$l_n",
                ylabel="integrand",
            )
        end

        # WYYKK_integral pages
        for iν in eachindex(ν), iμᶜ in eachindex(μᶜs), iμ in eachindex(μs), lᶜ_nᶜ in 0:lᶜ_nᶜ_max, l_n in 0:l_n_max
            l_n_slot = l_n + 1
            lᶜ_nᶜ_slot = lᶜ_nᶜ + 1

            save_report_plot(
                joinpath(reportDir, "40_WYYKK_integral_inu$(iν)_imuc$(iμᶜ)_imu$(iμ)_lc$(lᶜ_nᶜ)_l$(l_n).pdf"),
                WYYKK_integral,
                (l_n_slot, lᶜ_nᶜ_slot, iμ, iμᶜ, iν);
                N=80,
                title="∫WYYKK, iν=$iν, iμᶜ=$iμᶜ, iμ=$iμ, lᶜ-nᶜ=$lᶜ_nᶜ, l-n=$l_n",
                ylabel="integral",
            )
        end
    end

    
    return @strdict(WYYKK_integral=WYYKK_integral, reportDir=reportDir)

end


function save_report_plot(
    filepath::AbstractString,
    csf::CompactSymbolicFunctions,
    slots::NTuple;
    N=80,
    Δxval=1.0,
    title="",
    ylabel="value",
)
    fig = plotCompactSymbolicFunctions(csf, slots; N=N, Δxval=Δxval, ylabel=ylabel)
    ax = content(fig[1, 1])
    ax.title = title
    save(filepath, fig)
    return filepath
end

function save_bspline_report_plot(
    filepath::AbstractString,
    csf::CompactSymbolicFunctions;
    derivOrder=0,
    order=0,
    N=80,
    Δxval=1.0,
    title="",
)
    slots = (derivOrder + 1, order + 1)
    ylabel = derivOrder == 0 ? "B-spline value" : "Derivative order $derivOrder"
    return save_report_plot(filepath, csf, slots; N=N, Δxval=Δxval, title=title, ylabel=ylabel)
end

