using LinearAlgebra

function deviation(j::Int64, rhoA_copy::DensityMatrix,  set::Settings)::Float64
    @unpack observables, meas0 = set
    return (expect(observables[j], rhoA_copy) - meas0[j])
end 

function integrand(set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64})
    @unpack N_A, rhoA, observables, meas0, T_max, mtrxObs = set
    U = exp(-1im*t*H_A)
    rhoEvolved = U*rhoA.state*U'
    im_res = @inbounds sum((tr(mtrxObs[j]*rhoEvolved) - meas0[j])^2 for j in eachindex(observables))
    return real(im_res)
end


function cost(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    @unpack observables, T_max, q = set 
    H_A = @inbounds sum(g[i]*blks.matrices[i] for i in eachindex(g))
    return tanh_sinh(t -> integrand(set, t, H_A),0., T_max, q)[1]/(length(observables)*T_max)
end 


#################################################################################################
#                                                                                               #
#        old cost function using ThreadedDoubleExponentialFormula                               #
#                                                                                               #   
#################################################################################################

#=

function cost_for_grad(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    @unpack observables, T_max = set 
    H_A = @inbounds sum(g[i]*blks.matrices[i] for i in eachindex(g))
    return quadde(t -> integrand_for_grad(set, t, H_A),0,T_max, multithreading = false)[1]/(length(observables)*T_max)
end 

=#

#################################################################################################
#                                                                                               #
#        old cost function for seperate evaluation of the gradient and cost function            #
#                                                                                               #   
#################################################################################################
#=
function integrand(set::Settings, t::AbstractFloat, blks::H_A_Var)
    @unpack N_A, rhoA, observables, meas0, T_max= set
    rhoA_copy= copy(rhoA)
    blks.U.dt = t
    c = 0.0
    rhoA_copy |> blks.U
    for j in eachindex(observables)
        c+= deviation(j, rhoA_copy, set)^2
    end
    return c
end

function cost(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, multithreading::Bool)
    @unpack observables, T_max = set 
    H_A = sum(g[i]*blks.blocks[i] for i in eachindex(g))
    blks.U.H = H_A
    return quadde(t -> integrand(set, t, blks),0,T_max, multithreading = multithreading)[1]/(length(observables)*T_max)
end

=#


#####################################################################
#                                                                   #
#        worse integration techniques                               #
#                                                                   #
#####################################################################

#=


function cost_with_QuadGK(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, multithreading::Bool = false)
    @unpack observables, T_max = set 
    @unpack observables, T_max = set 
    H_A = sum(g[i]*blks.blocks[i] for i in eachindex(g))
    blks.U.H = H_A
    return quadgk(t -> integrand(set, t, blks), 0, T_max)[1]/(length(observables) * T_max)
end

=#

#=



function cost_with_QuadGK(g::Vector{<:AbstractFloat}, set::Settings,get_HA::Function = H_A_BW)
    @unpack observables, T_max = set 

    f(x) = integrand(g, set, x, get_HA)
    function f!(y::Vector{Float64}, x::Vector{Float64})
           n::Int = Threads.nthreads()
           Threads.@threads for i in 1:n
                y[i:n:end] .= f.(@view(x[i:n:end]))
           end
       end
    return quadgk(BatchIntegrand{<:AbstractFloat}(f!), 0, T_max, rtol = 1e-5)[1]/(length(observables) * T_max)
end


function cost_with_midpointrule(g::Vector{<:AbstractFloat}, set::Settings, get_HA::Function = H_A_BW)
    @unpack dt, N_T, N_A, rhoA, observables, meas0, T_max = set
    C::Float64 = 0.0
    rhoA_copy::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}} = copy(rhoA)
    H_A::Add{2} = get_HA(g, set)
    U_halve = time_evolve(cache(H_A), dt/2, check_hermicity=false)
    U = time_evolve(cache(H_A), dt, check_hermicity=false)
    Yao.apply!(rhoA_copy, U_halve)
    for j in eachindex(observables)
        C += deviation(j, rhoA_copy, set)^2
    end 
    for _ in 1:N_T-1
        Yao.apply!(rhoA_copy, U)
        for j in eachindex(observables)
            C += deviation(j, rhoA_copy, set)^2
        end
    end
    return C/(T_max*length(observables))*dt
end

#function wadj()
#    g = [1.,2.,3.]
#    blks = [put(3, i=>Z) for i in 1:3]
#    g'*blks
#end 
#
#function wmpadreduce()
#    g = [1.,2.,3.]
#    blks = [put(3, i=>Z) for i in 1:3]
#    mapreduce()
#end 

=#