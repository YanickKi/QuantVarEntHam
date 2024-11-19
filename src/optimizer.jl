using Optim

function optimize_LBFGS(g_init::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var; g1::AbstractFloat=NaN, gtol::AbstractFloat=1e-12
    , maxiter::Integer = 200, multithreading::Bool = false)
    @unpack T_max, N, N_A = set
    println("g_init: ", g_init)
    println("N: ", N)
    println("N_A: ", N_A)
    println("T_max: ", T_max)
    if isnan(g1)
        @assert length(g_init) == length(blks.blocks) "You entered $(length(g_init)) parameters but $(length(blks.blocks)) blocks. 
        The amount of parameters and blocks need to be equal!"
        AD_grad(g_init, set, blks, gtol, maxiter, multithreading)
    else 
        @assert length(g_init)+1 == length(blks.blocks) "You entered $(length(g_init)+1) parameters (from which is one fixed to $(g1)) but $(length(blks.blocks)) blocks.
        The amount of parameters and blocks need to be equal!"
        AD_grad_fixed(g_init, set, blks, g1, gtol, maxiter, multithreading)
    end
end


function AD_grad(g_init::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, gtol::AbstractFloat, maxiter::Integer, multithreading::Bool)
    
    result = optimize(Optim.only_fg!((F, G, g) -> fg!(F, G, g, set, blks)) ,g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = false,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    println(result)
    println(Optim.minimizer(result))
    return Optim.minimizer(result), cost_for_grad(Optim.minimizer(result), set, blks)
end

function AD_grad_fixed(g_init::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, g1::Float64, gtol::AbstractFloat, maxiter::Integer, multithreading::Bool)
    
    result = optimize(Optim.only_fg!((F, G, g) -> fg_fixed!(F, G, g, set, blks, g1)), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = false,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    println(result)
    println(Optim.minimizer(result))
    return Optim.minimizer(result), cost_for_grad(vcat(g1, Optim.minimizer(result)), set, blks)
end


function fg!(F::AbstractFloat, G::Vector{<:AbstractFloat}, g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    C, ∇ =  withgradient(g -> cost_for_grad(g, set, blks), g)
    if G !== nothing
        copyto!(G, ∇[1])
    end
    if F !== nothing
      return C
    end
end

function fg_fixed!(F::AbstractFloat, G::Vector{<:AbstractFloat}, g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var,  g1::AbstractFloat)
    println(set)
    C, ∇ =  withgradient(g -> cost_for_grad(vcat(g1,g), set, blks), g)
    if G !== nothing
        copyto!(G, ∇[1])
    end
    if F !== nothing
      return C
    end
end


#####################################################################
#                                                                   #
#           DOWN BELOW HERE GRAD AND C EVALUATED SEPERATELY         #
#                                                                   #
#####################################################################

#=
function AD_grad_old(g_init::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, gtol::AbstractFloat, maxiter::Integer, multithreading::Bool)
    result = optimize(g -> cost(g, set, blks, multithreading), (G, g)-> grad(G, g, set, blks),g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = true,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    println(result)
    println(Optim.minimizer(result))
    return Optim.minimizer(result), cost(Optim.minimizer(result), set, blks, multithreading)
end

function AD_grad_fixed_old(g_init::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, g1::Float64, gtol::AbstractFloat, maxiter::Integer, multithreading::Bool)
    result = optimize(g -> cost(vcat(g1, g), set, blks, multithreading), (G, g)-> grad_fixed(G, g, set, blks, g1),g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = true,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    println(result)
    println(Optim.minimizer(result))
    return Optim.minimizer(result), cost(vcat(g1, Optim.minimizer(result)), set, blks, multithreading)
end

function grad(G::Vector{Float64}, g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    G .= gradient(g -> cost_for_grad(g, set, blks), g)[1]
end 

function grad_fixed(G::Vector{Float64}, g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, g1::AbstractFloat)
    G .= gradient(g -> cost_for_grad(vcat(g1, g), set, blks), g)[1]
end 

=#



#####################################################
#                                                   #
#           DOWN BELOW HERE FOR FINITE DIFF         #
#                                                   #
#####################################################



#=


function _old_optimize_LBFGS(g_init::Vector{<:AbstractFloat}, set::Settings, H_A_Var::Function; g1::AbstractFloat=NaN, grad_parallel::Bool = false, gtol::AbstractFloat=1e-12
    , maxiter::Integer = 200)
    @unpack T_max, N, N_A = set
    println("g_init: ", g_init)
    println("N: ", N)
    println("N_A: ", N_A)
    println("T_max: ", T_max)

    if Threads.nthreads() > 1
        println('\n')
       @info "Multithreaded integration activated"
       println('\n')
    end

    if grad_parallel == true 
        if isnan(g1)
            parallel(g_init, set, H_A_Var, gtol, maxiter)
        else 
            parallel_fixed(g_init, g1, set, H_A_Var, gtol, maxiter)
        end
    else
        if isnan(g1)
            serial(g_init, set, H_A_Var, gtol, maxiter)
        else 
            serial_fixed(g_init, g1, set, H_A_Var, gtol, maxiter)
        end
    end 
end

function serial(g_init::Vector{<:AbstractFloat}, set::Settings, H_A_Var::Function, gtol::AbstractFloat=1e-12
    , maxiter::Integer = 200 )
    result = optimize(g -> cost(g, set, H_A_Var), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = true,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    println(result)
    println(Optim.minimizer(result))
    return Optim.minimizer(result), cost(Optim.minimizer(result), set, H_A_Var)
end



function serial_fixed(g_init::Vector{<:AbstractFloat}, g1::AbstractFloat, set::Settings, H_A_Var::Function, gtol::AbstractFloat=1e-12
    , maxiter::Integer = 200 )

    result = optimize(g -> cost(vcat(g1,g), set, H_A_Var), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = true,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    println(result)
    println(Optim.minimizer(result))
    return Optim.minimizer(result), cost(vcat(g1, Optim.minimizer(result)), set, H_A_Var)
end


function parallel(g_init::Vector{<:AbstractFloat}, set::Settings, H_A_Var::Function, gtol::AbstractFloat=1e-12
    , maxiter::Integer = 200 )

    h::Float64 = eps(Float64)^(1/3)
    Gshared::SharedArray{Float64} = zeros(length(g_init))
    
    result = optimize(g -> cost(g, set, H_A_Var), (G, g)-> g!(G, g, set, h, Gshared, H_A_Var), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = true,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    println(result)
    println(Optim.minimizer(result))
    return Optim.minimizer(result), cost(vcat(g1, Optim.minimizer(result)), set, H_A_Var)
end


function parallel_fixed(g_init::Vector{<:AbstractFloat}, g1::AbstractFloat, set::Settings, H_A_Var::Function, gtol::AbstractFloat=1e-12
    , maxiter::Integer = 200 )


    h::Float64 = eps(Float64)^(1/3)
    Gshared::SharedArray{Float64} = zeros(length(g_init))
    
    result = optimize(g -> cost(vcat(g1,g), set, H_A_Var), (G, g) -> g_fixed!(G, g, g1, set, h, Gshared, H_A_Var), g_init, LBFGS(), Optim.Options(g_tol = gtol,
                                                                    store_trace = true,
                                                                    show_trace = true,
                                                                    show_warnings = true, iterations = maxiter))

    println(result)
    println(Optim.minimizer(result))
    return Optim.minimizer(result), cost(vcat(g1, Optim.minimizer(result)), set, H_A_Var)
end


function add_inOneDirec(a::Vector{<:AbstractFloat}, index::Integer, h::Float64)
        b::Vector{Float64} = copy(a)
        b[index]+=h
    return b 
end

function finite_diff(g::Vector{<:AbstractFloat}, index::Integer, set::Settings, h::Float64, blks::H_A_Var)
    return (cost(add_inOneDirec(g, index, h), set, blks) - cost(add_inOneDirec(g, index, -h), set, blks))/(2*h)
end

function g!(G::Vector{<:AbstractFloat}, g::Vector{<:AbstractFloat}, set::Settings, h::Float64, Gshared::SharedArray{<:AbstractFloat}, blks::H_A_Var)
    Gshared = pmap(index -> finite_diff(g, index, set, h, H_A_Var), eachindex(g))
    G .= Gshared     
end 

function g_fixed!(G::Vector{<:AbstractFloat}, g::Vector{<:AbstractFloat}, g1::AbstractFloat, set::Settings, h::Float64, blks::H_A_Var)
    G = map(index -> finite_diff(vcat(g1, g), index+1, set, h, blks), eachindex(g))
    return G 
end

#####################################################
#                                                   #
#           DOWN BELOW HERE MULTITHREADED           #
#                                                   #
#####################################################

function g_fixed_backup!(G::Vector{Float64}, g::Vector{Float64}, g1::AbstractFloat, set::Settings, h::Float64, G_padded::Vector{Float64}, H_A_Var::Function)
    
    Threads.@threads for index in eachindex(g)
        G_padded[index*8-7] = finite_diff(vcat(g1, g), index+1, set, h, H_A_Var)
    end 
    for index in eachindex(g)
        G[index] = G_padded[index*8-7]
    end
end 

=#