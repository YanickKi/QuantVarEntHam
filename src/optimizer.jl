using Optim, LinearAlgebra

"""
    optimize_LBFGS(g_init::Vector{<:AbstractFloat}, init::Init; g1::AbstractFloat=NaN, gtol::AbstractFloat=1e-12, maxiter::Integer = 200, show_trace::Bool=true, print_result::Bool = true)

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
function optimize_LBFGS(g_init::Vector{<:Real}, init::Init; methods = tanh_sinh(), g1::Real=NaN, gtol::AbstractFloat=1e-16, maxiter::Integer = 200, show_trace::Bool=true, print_result::Bool = true)

    _g_init = Float64.(g_init)
    _g_1 = Float64(g1)

    if isnan(g1)
        @assert length(_g_init) == length(init.blks.blocks) "You entered $(length(_g_init)) parameters but $(length(init.blks.blocks)) blocks. 
        The amount of parameters and blocks need to be equal!"
        optimize_free(_g_init, init, gtol, maxiter, show_trace, print_result, methods)
    else 
        @assert length(_g_init)+1 == length(init.blks.blocks) "You entered $(length(_g_init)+1) parameters (from which is one fixed to $(_g_1)) but $(length(init.blks.blocks)) blocks.
        The amount of parameters and blocks need to be equal!"
        optimize_fixed(_g_init, init, _g_1, gtol, maxiter, show_trace, print_result, methods)
    end
end


function optimize_free(g_init::Vector{<:AbstractFloat}, init::Init, gtol::AbstractFloat, maxiter::Integer, show_trace::Bool, print_result::Bool, methods)
    
    result = optimize(Optim.only_fg!((F, G, g) -> cost_grad!(F, g, G , init, methods)), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                store_trace = false,
                                                                show_trace = show_trace,
                                                                show_warnings = true, iterations = maxiter))
    g_opt = Optim.minimizer(result)                                                                
    if print_result == true                                                   
        println(result)
        println(g_opt)
    end 
    return g_opt, cost_grad!(1.,  g_opt, nothing, init, methods) 
end



function optimize_fixed(g_init::Vector{<:AbstractFloat}, init::Init, g1::Float64, gtol::AbstractFloat, maxiter::Integer, show_trace::Bool, print_result::Bool, methods)

    result = optimize(Optim.only_fg!((F, G, g) -> cost_grad_fixed!(F, vcat(g1,g), G , init, methods)), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = false,
                                                                    show_trace = show_trace,
                                                                    show_warnings = true, iterations = maxiter))

    g_opt = Optim.minimizer(result)             
    if print_result == true                                                   
        println(result)
        println(g_opt)
    end 
    return g_opt, cost_grad_fixed!(1.,  vcat(g1,g_opt), nothing, init, methods)
end

function optimize_relativeEntropy(g_init::Vector{<:AbstractFloat}, init::Init; gtol::AbstractFloat=1e-12, maxiter::Integer = 200, show_trace::Bool=true, print_result::Bool = true)

    function rel_entropy(g::Vector{<:AbstractFloat}, init::Init)
        @unpack ρ_A  = init.set
        get_H_A!(g, init) 
        S1 = tr(ρ_A.state * init.buff.H_A)
        S2 = log(tr(exp(-init.buff.H_A)))
        return S1 + S2
    end 

    function grad_relentropy(G, g::Vector{<:AbstractFloat}, init::Init)
        @unpack ρ_A  = init.set
        get_H_A!(g, init) 
        expHA = exp(- init.buff.H_A)
        for i in eachindex(g)
            dS1 = tr(ρ_A.state*init.blks.matrices[i])
            dS2 = 1/(tr(expHA)) * tr(-expHA*init.blks.matrices[i])
            G[i] = dS1+dS2
        end 
    end 
    result = optimize(g ->  rel_entropy(g, init), (G, g) ->  grad_relentropy(G, g, init) ,g_init, BFGS(), Optim.Options(g_tol = gtol,
        store_trace = false,
        show_trace = show_trace,
        show_warnings = true, iterations = maxiter))

    g_opt = Optim.minimizer(result)                                                                
    if print_result == true                                                   
        println(result)
        println(g_opt)
    end 
    return g_opt, rel_entropy(g_opt, init)

end 


function comm_opt_fixed(g_init::Vector{<:AbstractFloat}, init::Init, g1::AbstractFloat, gtol::AbstractFloat, maxiter = 100)
    result = optimize(g -> comm_cost(vcat(g1,g),init), g_init, LBFGS(), Optim.Options(g_tol = gtol,
    store_trace = false,
    show_trace = true,
    show_warnings = true, iterations = maxiter))

    g_opt = Optim.minimizer(result)

    println(result)
    println(g_opt)
    return g_opt, comm_cost(vcat(g1,g_opt),init)
end 

function comm_opt(g_init::Vector{<:AbstractFloat}, init::Init; g1::AbstractFloat=NaN, gtol::AbstractFloat = 1e-12)
    result = optimize(g -> comm_cost(g,init), g_init, LBFGS(), Optim.Options(g_tol = gtol,
    store_trace = false,
    show_trace = true,
    show_warnings = true, iterations = 10000))

    g_opt = Optim.minimizer(result)

    println(result)
    println(g_opt)
    return g_opt, comm_cost(g_opt,init)
end 


# FREE TMAX

function optimize_T_max_g1_fixed_custrgad(g_init::Vector{<:AbstractFloat}, init::Init, g1::AbstractFloat; gtol::AbstractFloat = 1e-12, maxiter::Integer = 200)
    
    result = optimize(Optim.only_fg!((F, G, g) -> cost_grad_Tmax_g1_fixed!(F, vcat(g[1],g1,g[2:end]), G , init)), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = false,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    g_opt = Optim.minimizer(result)                                                                
    println(result)
    println(g_opt)
    return g_opt, cost_grad_Tmax_g1_fixed!(1.,  vcat(g_opt[1],g1,g_opt[2:end]), nothing, init)
end

function optimize_midpoint(g_init::Vector{<:AbstractFloat}, init::Init; gtol::AbstractFloat = 1e-12, maxiter::Integer=200)
    
    result = optimize(Optim.only_fg!((F, G, g) -> cost_grad_midpoint!(F, g, G , init)), g_init, BFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = false,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    g_opt = Optim.minimizer(result)                                                                
    println(result)
    println(g_opt)
    return g_opt, cost_grad_midpoint!(1.,  g_opt, nothing, init)
end
