using SparseArrays,Symbolics #SparseDiffTools

include("../src/batchUseful.jl")
include("../src/batchEnzyme.jl")

function buildNumericalFunctions(costfunctions, symbUnknownField, symbKnownField, symbKnownForce)
    # Collect all known inputs (constants) and unknown inputs (variables)
    knownInputs = vcat(reduce(vcat, reduce(vcat, symbKnownField;init=[]);init=[]), reduce(vcat, symbKnownForce))
    unknownInputs = reduce(vcat, symbUnknownField)
    all_inputs = vcat(unknownInputs, knownInputs)
    
    residual_func = Array{Function, 1}(undef, length(costfunctions))
    #@show costfunctions[length(costfunctions)÷2],symbKnownForce
    #@show costfunctions[44],costfunctions[45],costfunctions[46]
    #@show costfunctions[54],costfunctions[55],costfunctions[56]
    #@show costfunctions[64],costfunctions[65],costfunctions[66]
    for i in eachindex(costfunctions)
        # Create the symbolic function
        residual_func_expr = build_function(costfunctions[i], all_inputs; expression = Val{false})
        
        # Evaluate the symbolic function and store it as a numerical function
        residual_func[i] = eval(residual_func_expr)
    end
    #@show residual_func[120]
    return residual_func
end


function Residual_OLD!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)

    mapping = Dict()
    
    for k in eachindex(knownField)
        for j in eachindex(knownField[k])
            #mapping[symbKnownField[k][j]] = knownField[k][j]
        end
    end

    for j in eachindex(unknownField)
        #mapping[symbUnknownField[j]] = unknownField[j]
    end

    for j in eachindex(knownForce)
        #mapping[symbKnownForce[j]] = knownForce[j]
    end

    for i in eachindex(F)
        F[i] = substitute(costfunctions[i],mapping)
    end
    return
end

function makeInputsForNumericalFunctions(unknownField,knownField,knownForce)
    knownInputs = vcat(reduce(vcat,reduce(vcat,knownField; init=[]); init=[] ),reduce(vcat,reduce(vcat,knownForce)))
    unknownInputs = reduce(vcat,unknownField)
    #all_inputs = vcat(unknownInputs,knownInputs)
    #@show length(knownInputs),length(unknownInputs)
    return unknownInputs,knownInputs
end

function Residual!(F,f,unknownInputs,knownInputs)
    all_inputs=vcat(unknownInputs,knownInputs)
    for i in eachindex(f)
        F[i]=f[i]((all_inputs))
    end
end
function temporaryConstantResidualFunction(f::Vector{F}, U::AbstractVector, knownInputs::AbstractVector;boundaryConditionForced=false) where {F <: Function}
    N_U = length(U)
    N_K = length(knownInputs)
    function f_specific!(Fout::AbstractVector{T}, unknownInputs::AbstractVector{T}) where {T}
        all_inputs = Vector{T}(undef, N_U + N_K)
        @inbounds begin
            all_inputs[1:N_U] .= unknownInputs
            for i in 1:N_K
                all_inputs[N_U + i] = T(knownInputs[i])
            end
            for i in eachindex(f)
                Fout[i] = f[i](all_inputs)
            end
            if boundaryConditionForced
                
                Fout[1]=unknownInputs[1]-1.0
                Fout[end]=unknownInputs[end]-1.0
            end
        end
        return nothing
    end
    return f_specific!
end


function sparseColouring(f,unknownField,knownField,knownForce)
    # this function is colouring the matrix
    knownField .= 0.0
    knownForce .= 0.0
    unknownField .=0.0
    U,knownInputs = makeInputsForNumericalFunctions(unknownField,knownField,knownForce)

    nCostfunctions = length(f)
    nU = length(U)
    U = rand(nU)
    F=rand(nCostfunctions)
    Res_closed_look! = (F,U) -> Residual!(F,f,U,knownInputs)
    #Res_closed_look! = temporaryConstantResidualFunction(f,U,knownInputs)
    sparsity    = Symbolics.jacobian_sparsity(Res_closed_look!,F,U)
    J           = Float64.(sparse(sparsity))
    V=rand(nU)
    cache=ForwardColorJacCache(Res_closed_look!, V)
    #J = (sparse(sparsity))
    #colors      = matrix_colors(J)
    #coloring_cache = ForwardColorJacCache(Res_closed_look!, J, U; sparsity = J)
    return J,cache
end


function timeStepOptimisation!(f,unknownField,knownField,knownForce,J,cache,NpointsSpace,NField;nIteration=10,smallNumber =1.e-8,boundaryConditionForced=false)


    #danger !!! 
    boundaryConditionForced = true
    #danger!!!
    # DANGER!!!



    nEq = length(f)    
    # normalisation by the number of equations
    normalisation = 1.0/nEq
    r1 = 1.0
    unknownField .= 0.0
    U,knownInputs = makeInputsForNumericalFunctions(unknownField,knownField,knownForce)
    #@show maximum(knownField)

    δU = U
    F=zeros(nEq)
    f_specific! = temporaryConstantResidualFunction(f,U,knownInputs;boundaryConditionForced=boundaryConditionForced)
    for iter in 1:nIteration
        
        #Residual!(F,costfunctions,symbUnknownField,unknownField,symbKnownField,knownField,symbKnownForce,knownForce)
        #Res_closed! = (F,U) -> Residual!(F,f,U,knownInputs)
        #Res_closed!(F,U)

        f_specific!(F,U)

      
        r = norm(F)*normalisation
        @show iter, r
        if iter==1 r1 = r; end
        if r === 0.0 break end
        if r/r1 < smallNumber break end
      
        # Jacobian assembly
    
        @time forwarddiff_color_jacobian!(J, f_specific!, U, cache)
        #@time handMadeJacobianComputation!(J, f_specific!, F, U, rows, cols)

        # Solve
        #@time factor = lu(J)  # Or try `ldlt`, `cholesky`, or `qr` depending on J's properties
        #invJac=inv(factor)
        invJac = myInv(J)
        @show round.(F,sigdigits=2)

        @time δU = - invJac * F
        #@time δU   .= .-J\F
        @show round.(δU,sigdigits=2)

        α = 1.0
        U    .+= α .* δU
        
    end
    #@show maximum(U)
    
    unknownField .= reshape(U,NpointsSpace,NField)
    
    return 
end


function handMadeJacobianComputation!(J, f_specific!, F, U, rows, cols)
    # Define the residual wrapper function to give to Enzyme
    #residual_wrapper(F::Vector{Float64}, U::Vector{Float64}) = f_specific(F, U)

    # perturbation 
    ΔU = maximum((1.e-8,maximum(U)))

    # Loop through the non-zero entries
    for k in eachindex(rows)
        i = rows[k]
        j = cols[k]

        # Perturbation vector
        dU = zeros(length(U))
        dU[j] = ΔU

        # dF: output for directional derivative
        #dF = zeros(length(F))
        
        # Call Enzyme's forward-mode autodiff on the wrapper function
        #Enzyme.autodiff(Enzyme.Forward,residual_wrapper,Const, Duplicated(F, dF), Duplicated(U, dU))
        newU = U .+ dU
        newF = similar(F)
        f_specific!(newF, newU)
        dF_i = newF[i]-F[i]

        # Store entry in the sparse Jacobian
        J[i, j] = dF_i/dU[j]
    end
end







# detection of backend
_has_cuda, _has_metal = detect_backend()

