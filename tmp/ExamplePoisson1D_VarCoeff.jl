using Plots, Printf # need to install Plots package 
import LinearAlgebra: norm

𝑢(𝑥)    =  cos(𝑥) # manufactured solution
∇𝑢(𝑥)   = -sin(𝑥) 
Δ𝑢(𝑥)   = -cos(𝑥)
β(𝑥)    =  cos(𝑥)
∇β(𝑥)   = -sin(𝑥)
# β(𝑥)    =  1.0
# ∇β(𝑥)   =  0.0
𝑞(𝑥)    = -β(𝑥)*∇𝑢(𝑥) 
∇𝑞(𝑥)   = -(∇β(𝑥)*∇𝑢(𝑥) + β(𝑥)*Δ𝑢(𝑥))

@views function main(ncx, method) 

    L   = 1.0 # Length

    # Parameters
    x   = (min=-L/2, max=L/2)
    Δx  = (x.max - x.min) / ncx
    xce = LinRange(x.min-Δx/2., x.max+Δx/2., ncx+2)
    xve = LinRange(x.min-Δx, x.max+Δx, ncx+3)

    # Allocations
    u      = zeros(ncx+2)
    u_ex   = zeros(ncx+2)
    b      = zeros(ncx+2)
    Δb     = zeros(ncx+2)
    k      = zeros(ncx+3)
    dudx   =  ones(ncx+3)
    q      =  ones(ncx+3)
    q_ex   =  ones(ncx+3)
    dqdx   =  ones(ncx+2)
    F_dudx = zeros(ncx+1) # residual dudx
    ∂a∂τ   = zeros(ncx+1)
    F_dqdx = zeros(ncx-0) # residual dqdx
    ∂b∂τ   = zeros(ncx-0)
    F      = zeros(ncx-0) # residual Poisson
    ∂u∂τ   = zeros(ncx-0) 

    qb     =  ones(ncx+3)
    Δdudx  =  ones(ncx+3)
    dudx   =  ones(ncx+3)

    # Initial condition: Evaluate exact initial solution
    for i in eachindex(u)
        u[i]    =  𝑢(xce[i])
        u_ex[i] =  𝑢(xce[i])
        b[i]    =-∇𝑞(xce[i])
        dqdx[i] = ∇𝑞(xce[i])
    end
    for i in eachindex(q)
        dudx[i] = ∇𝑢(xve[i])
        q_ex[i] =  𝑞(xve[i])
        q[i]    =  𝑞(xve[i])
        k[i]    = β.(xve[i])
    end

    # Pseudo-transient integration
    Δτ = Δx^2/maximum(k)/2.1 
    θ  = 0.038/Δx/1.6 

    if method===:Spotz         # Discretisation of Spotz et al., 1995
        for iter=1:10000
            # BC
            # u[1]   = 2* 𝑢(-L/2)   - u[2]       # this is only second order - do not use
            # u[end] = 2* 𝑢( L/2)   - u[end-1]
            # b[1]   = 2*∇𝑞(-L/2,β) - b[2]       # this is OK as the Δb is anyway O(2)
            # b[end] = 2*∇𝑞( L/2,β) - b[end-1]
            # Flux
            dudx[2:end-1]  .= -1
            q[2:end-1]  .= -k[2:end-1].*diff(u)/Δx
            qb[2:end-1] .= -k[2:end-1].*diff(b)/Δx
            # Laplacian of source term
            Δb[2:end-1] .= (b[1:end-2] .+ b[3:end-0] .- 2*b[2:end-1])/Δx^2
            Δdudx[2:end-1] .= -(dudx[1:end-2] .+ dudx[3:end-0] .- 2*dudx[2:end-1])/Δx^2
            # Δb[2:end-1] .= -diff(qb[2:end-1])/Δx
            # ΔqE[2:end-1] .= (b[1:end-2] .+ b[3:end-0] .- 2*b[2:end-1])/Δx^2
            # Residual
            c2 = Δx^2/12 .*(k[3:end-1].-k[2:end-2])/Δx .* 0.5.*(Δdudx[3:end-1] .+ Δdudx[2:end-2])
            # c2 = Δx/12 .*(k[3:end-1].*Δdudx[3:end-1] .- k[2:end-2].*Δdudx[2:end-2]) 
            F          .= -diff(q[2:end-1])/Δx - b[2:end-1] .- Δx^2/12 .* Δb[2:end-1] #.- c2
            # Check
            nF          = norm(F)/sqrt(length(F))
            nF < 1e-12 ? break : nothing
            mod(iter, 100)==0 ? @printf("Iter %05d --- abs. res. = %1.4e\n", iter, nF) : nothing 
            # Rate update
            ∂u∂τ       .= (1-θ).*∂u∂τ .+ F 
            # Variable update
            u[2:end-1] .+= Δτ*∂u∂τ
        end

    elseif method===:Abide # Discretisation of Abide et al., 2020
        α = 1/22
        a = 12/11
        # α = 0.
        # a = 1.0    
        dudx[2:end-1] .= 0
        dqdx[2:end-1] .= 0
        for iter=1:100000
            # BC
            # # ---> nothing, I assume we know the functions at ghost nodes locations, as above
            # # Residuals
            # F_dudx     .= -(α*dudx[1:end-2] .+ dudx[2:end-1] .+ α*dudx[3:end-0] .- a*(u[2:end]   .- u[1:end-1])/Δx)
            # F_dqdx     .= -(α*dqdx[1:end-2] .+ dqdx[2:end-1] .+ α*dqdx[3:end-0] .- a*(q[3:end-1] .- q[2:end-2])/Δx)
            # F          .= -dqdx[2:end-1]    .- b[2:end-1]
            # q[2:end-1] .= -k[2:end-1].*dudx[2:end-1]
            # # Check
            # nF          = norm(F)/sqrt(length(F))
            # nF < 1e-11 ? break : nothing
            # mod(iter, 100)==0 ? @printf("Iter %05d --- abs. res. = %1.4e %1.4e %1.4e\n", iter, nF, norm(F_dudx), norm(F_dqdx)) : nothing 
            # # Rate update
            # ∂a∂τ       .= (1-θ*0.4).*∂a∂τ .+ F_dudx
            # ∂b∂τ       .= (1-θ*0.4).*∂b∂τ .+ F_dqdx
            # ∂u∂τ       .= (1-θ*0.025).*∂u∂τ .+ F 
            # # Variable update
            # dudx[2:end-1] .+= ∂a∂τ*1.7
            # dqdx[2:end-1] .+= ∂b∂τ*1.7
            # u[2:end-1]    .+= ∂u∂τ*1e-5/6

            #######################################################

            # BC
            # ---> nothing, I assume we know the functions at ghost nodes locations, as above
            # Residuals
            F_dudx     .= -(α*dudx[1:end-2] .+ dudx[2:end-1] .+ α*dudx[3:end-0] .- a*(u[2:end]   .- u[1:end-1])/Δx)
            # Rate update
            ∂a∂τ       .= (1-θ*0.95).*∂a∂τ .+ F_dudx
            # Variable update
            dudx[2:end-1] .+= ∂a∂τ*1.
            q[2:end-1] .= -k[2:end-1].*dudx[2:end-1]

            # Residuals
            F_dqdx     .= -(α*dqdx[1:end-2] .+ dqdx[2:end-1] .+ α*dqdx[3:end-0] .- a*(q[3:end-1] .- q[2:end-2])/Δx)
            # Rate update
            ∂b∂τ       .= (1-θ*0.95).*∂b∂τ .+ F_dqdx
            # Variable update
            dqdx[2:end-1] .+= ∂b∂τ*1.

            # Residuals
            F          .= -dqdx[2:end-1]    .- b[2:end-1]
            # Rate update
            ∂u∂τ       .= (1-θ*0.01).*∂u∂τ .+ F 
            # Variable update
            u[2:end-1]    .+= ∂u∂τ*0.8e-5

            # Check
            nF          = norm(F)/sqrt(length(F))
            nF < 1e-11 ? break : nothing
            mod(iter, 100)==0 ? @printf("Iter %05d --- abs. res. = %1.4e %1.4e %1.4e\n", iter, nF, norm(F_dudx), norm(F_dqdx)) : nothing 
           
        end
        @show Δx
    end

    # Discretisation error
    err   = norm(u_ex .- u)/sqrt(length(u))

    # Plot
    p1 = plot(title = @sprintf("u error = %1.2e", err), xlabel="x", ylabel="u")
    plot!(xve, k, label="coefficient")
    p2 = plot(title = @sprintf("u error = %1.2e", err), xlabel="x", ylabel="u")
    plot!(xce, u_ex, label="exact")
    scatter!(xce, u, label="num")
    display(plot(p1, p2, layout=(2,1)))

    return err
end

# @printf("Convergence for Spotz (1996):\n")
# @show sqrt(main(10, :Spotz)/main(20, :Spotz))
# @show sqrt(main(20, :Spotz)/main(40, :Spotz))

@printf("Convergence for Abide (2020):\n")
@show sqrt(main(10, :Abide)/main(20, :Abide))
@show sqrt(main(20, :Abide)/main(40, :Abide))
