using Optim, LinearAlgebra

"""
    optimize_LBFGS(g_init::Vector{<:Real}, init::Init; cost::Symbol = :QCFL, integration_method = tanh_sinh(), g1::Real=NaN, gtol::AbstractFloat=1e-16, maxiter::Integer = 200, show_trace::Bool=true, print_result::Bool = true)

Minimize cost function using the LBFGS-optimizer from `Optim.jl`.

# Required arguments
- `g_init::Vector{<:AbstractFloat}`: initial parameters.
- `init::Init`: struct containing settings, variational Ansatz, buffers and integration table for tanh-sinh integration.

# Keyword arguments 
- `g1::AbstractFloat=NaN`: provide number if the first parameter should be fixed.
- `gtol::AbstractFloat=1e-12`: maximum gradient norm in order to stop minimizing.
- `maxiter::Integer = 200`: maximum number of iterations in order to stop minimizing.
- `show_trace::Bool=true`: true for showing the trace of the minizing procedure, false otherwise.
- `print_result::Bool = true` true to print optimal parameters, false otherwise 
"""
function optimize_LBFGS(g_init::Vector{<:Real}, init::Init; cost = QCFL(), gtol::AbstractFloat=1e-16, maxiter::Integer = 200, show_trace::Bool=true, print_result::Bool = true)

   
    cost_function, g1 = cost
    
    _g_init = Float64.(g_init)

    if isnan(g1)
        @assert length(_g_init) == length(init.blocks) "You entered $(length(_g_init)) parameters but $(length(init.blocks)) blocks. 
        The amount of parameters and blocks need to be equal!"
        optimize_free(_g_init, init, gtol, maxiter, show_trace, print_result, cost_function)
    else 
        num_free_params = length(g_init)
        @assert num_free_params+1 == length(init.blocks) "You entered $(num_free_params+1) parameters (from which is one fixed to $(g1)) but $(length(init.blocks)) blocks.
        The amount of parameters and blocks need to be equal!"
        
        if num_free_params +1 != length(init.buff.C_G) || num_free_params +1 != length(init.buff.C_G_result)
            pop!(init.buff.C_G)
            pop!(init.buff.C_G_result)
        end 

        optimize_fixed(_g_init, init, g1, gtol, maxiter, show_trace, print_result, cost_function)
    end
end

function optimize_free(g_init::Vector{<:AbstractFloat}, init::Init, gtol::AbstractFloat, maxiter::Integer, show_trace::Bool, print_result::Bool, cost_function::Function)
    
    result = optimize(Optim.only_fg!((F, G, g) -> cost_function(F, G, g, init)), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                store_trace = false,
                                                                show_trace = show_trace,
                                                                show_warnings = true, iterations = maxiter))
    g_opt = Optim.minimizer(result)                                                                
    if print_result == true                                                   
        println(result)
        println(g_opt)
    end 
    return g_opt, cost_function(1., nothing,  g_opt, init) 
end



function optimize_fixed(g_init::Vector{<:AbstractFloat}, init::Init, g1::Float64, gtol::AbstractFloat, maxiter::Integer, show_trace::Bool, print_result::Bool, cost_function::Function)

    result = optimize(Optim.only_fg!((F, G, g) -> cost_function(F, G, vcat(g1,g), init)), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = false,
                                                                    show_trace = show_trace,
                                                                    show_warnings = true, iterations = maxiter))

    g_opt = Optim.minimizer(result)             
    if print_result == true                                                   
        println(result)
        println(g_opt)
    end 
    return g_opt, cost_function(1., nothing,  vcat(g1,g_opt), init)
end