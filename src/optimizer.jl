using Optim, LinearAlgebra

export optimize

"""
    optimize(cost::AbstractCostFunction, g_init::Vector{Real}; ∇_tol::Real=1e-16, maxiter::Integer=1000, show_trace::Bool=true, print_result::Bool= true)

Minimize cost function using the [LBFGS-optimizer from `Optim.jl`](https://julianlsolvers.github.io/Optim.jl/stable/algo/lbfgs/) and return 
the optimal parameters and the cost function value at the minimum.

# Required arguments
- `cost`: cost function.
- `g_init`: initial parameters.

# Keyword arguments 
- `∇_tol`: minimum infinity norm of the gradient in order to stop the minimization.
- `maxiter`: maximum number of iterations in order to stop minimizing.
- `show_trace`: true for showing the trace of the minizing procedure, false otherwise.
- `print_result` true to print optimal parameters and result, false otherwise.
"""
function optimize(cost::AbstractCostFunction, g_init::Vector{<:Real}; ∇_tol::Real=1e-16, maxiter::Integer=1000, show_trace::Bool=true, print_result::Bool= true)
    
    free_indices = get_free_indices(cost)

    @assert length(g_init) == length(free_indices) "You provided $(length(g_init)) initial parameters but there are $(length(free_indices)) free parameters!"  

    result = Optim.optimize(Optim.only_fg!((F, G, g) ->  fg!(F, G, cost, g)), Float64.(g_init), LBFGS(), Optim.Options(g_tol = ∇_tol,
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
