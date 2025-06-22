using Optim, LinearAlgebra

"""
    optimize_LBFGS(g_init::Vector{<:Real}, init::Init; cost = QCFL(), ∇_tol::Real=1e-16, maxiter::Integer = 1000, show_trace::Bool=true, print_result::Bool = true)

Minimize cost function using the [LBFGS-optimizer from `Optim.jl`](https://julianlsolvers.github.io/Optim.jl/stable/algo/lbfgs/) and return 
the optimal parameters and the cost function value at the minimum.

# Required arguments
- `g_init::Vector{<:Real}`: initial parameters.
- `init::Init`: struct containing settings, variational Ansatz and buffers.

# Keyword arguments 
- `cost = QCFL(1)`: cost function to be minimized, default is QCFL with ``T_\\text{max} =1``.
- `∇_tol::Real=1e-16`: maximum gradient norm in order to stop minimizing.
- `maxiter::Integer = 1000`: maximum number of iterations in order to stop minimizing.
- `show_trace::Bool=true`: true for showing the trace of the minizing procedure, false otherwise.
- `print_result::Bool = true` true to print optimal parameters and result, false otherwise.
"""
function optimize_LBFGS(g_init::Vector{<:Real}, cost; ∇_tol::AbstractFloat=1e-16, maxiter::Integer = 200, show_trace::Bool=true, print_result::Bool = true)

    cost_function, g1, numBlocks = cost

    _∇_tol = Float64(∇_tol)
    _g_init = Float64.(g_init)

    if isnan(g1)
        @assert length(_g_init) == numBlocks "You entered $(length(_g_init)) parameters but $(numBlocks) blocks. 
        The amount of parameters and blocks need to be equal!"
        optimize_free(_g_init, _∇_tol, maxiter, show_trace, print_result, cost_function)
    else 
        num_free_params = length(g_init)
        @assert num_free_params+1 == numBlocks "You entered $(num_free_params+1) parameters (from which is one fixed to $(g1)) but $(numBlocks) blocks.
        The amount of parameters and blocks need to be equal!"
        
        #if num_free_params +1 != length(init.buff.C_G) || num_free_params +1 != length(init.buff.C_G_result)
        #    pop!(init.buff.C_G)
        #    pop!(init.buff.C_G_result)
        #end 

        optimize_fixed(_g_init, g1, _∇_tol, maxiter, show_trace, print_result, cost_function)
    end
end

function optimize_free(cost::AbstractCostFunction, g_init::Vector{<:AbstractFloat}; ∇_tol::AbstractFloat=1e-16, maxiter::Integer=1000, show_trace::Bool=true, print_result::Bool= true)
    
    result = optimize(g -> cost(g), (G,g) -> gradient!(G, cost, g), g_init, LBFGS(), Optim.Options(g_tol = ∇_tol,
                                                                store_trace = false,
                                                                show_trace = show_trace,
                                                                show_warnings = true, iterations = maxiter))
    g_opt = Optim.minimizer(result)                                                                
    if print_result == true                                                   
        println(result)
        println(g_opt)
    end 
    return g_opt, cost(g_opt) 
end



function optimize_fixed(g_init::Vector{<:AbstractFloat},  g1::Float64, ∇_tol::AbstractFloat, maxiter::Integer, show_trace::Bool, print_result::Bool, cost_function::Function)

    result = optimize(Optim.only_fg!((F, G, g) -> cost_function(F, G, vcat(g1,g))), g_init, LBFGS(), Optim.Options(g_tol = ∇_tol,
                                                                    store_trace = false,
                                                                    show_trace = show_trace,
                                                                    show_warnings = true, iterations = maxiter))

    g_opt = Optim.minimizer(result)             
    if print_result == true                                                   
        println(result)
        println(g_opt)
    end 
    return g_opt, cost_function(1., nothing,  vcat(g1,g_opt))
end