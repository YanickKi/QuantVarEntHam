using LinearAlgebra
using ChainRules, ChainRulesCore


function get_H_A!(g::Vector{<:AbstractFloat}, init::Init)
    fill!(init.buff.H_A, 0)
    @fastmath @inbounds @simd for i in eachindex(g)
        init.buff.H_A .+= g[i].*init.blks.matrices[i]
    end 
end 

function cost_grad!(F::Union{AbstractFloat, Nothing}, g::Vector{<:AbstractFloat}, G::Union{Vector{<:AbstractFloat}, Nothing}, init::Init)
    @unpack observables, T_max, atol, rtol = init.set
    get_H_A!(g, init)
    if G !== nothing
        init.buff.C_G_result .= (tanh_sinh(t -> integrand(t, init),0., T_max, init.q, atol = atol, rtol = rtol)/(length(observables)*T_max))
        G[:] .= @view  init.buff.C_G_result[2:end]
    end
    if F !== nothing 
        if G === nothing
            return tanh_sinh(t -> integrand_onlycost(t, init),0., T_max, init.q, atol = atol, rtol = rtol)/(length(observables)*T_max)
        else 
            return init.buff.C_G_result[1]
        end 
    end
end 

function integrand(t::AbstractFloat, init::Init)
    @unpack ρ_A, meas0, mtrxObs, T_max = init.set
    buff = init.buff
    buff.C_G[1] = 0.
    mul!(buff.H_A_forexp, -1im*t, buff.H_A)
    U, pullback =  rrule(exp, buff.H_A_forexp)
    mul!(buff.ρ_A_right, ρ_A.state, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)

    @fastmath @inbounds @simd for i in eachindex(mtrxObs)
        mul!(buff.evolobs[i], mtrxObs[i], buff.ρ_A_evolved)
        buff.dev[i] = real(tr(buff.evolobs[i])) - meas0[i]
        buff.sumobs .+= buff.dev[i].*mtrxObs[i]
        buff.C_G[1] += buff.dev[i]^2
    end 

    mul!(buff.dAforpb, buff.ρ_A_right, buff.sumobs)

    adj = pullback(buff.dAforpb')[2]
    
    @fastmath @inbounds @simd for i in eachindex(buff.G_buffer)
        mul!(buff.adjtimesBlockMatrices[i], adj', init.blks.matrices[i])
    end     

    @fastmath @inbounds @simd for i in eachindex(buff.G_buffer)
        buff.G_buffer[i] = 4*t*imag(tr(buff.adjtimesBlockMatrices[i]))
    end 

    fill!(buff.sumobs, 0)
    buff.C_G[2:end] .= @view buff.G_buffer[:]
    return buff.C_G
end


function integrand_onlycost(t::AbstractFloat,  init::Init)
    @unpack ρ_A, meas0, mtrxObs = init.set
    buff = init.buff 
    mul!(buff.H_A_forexp, -t, buff.H_A)
    U = cis(buff.H_A_forexp)
    mul!(buff.ρ_A_right, ρ_A.state, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)
    c::Float64 = 0.
    for i in eachindex(mtrxObs)
        @inbounds mul!(buff.evolobs[i], mtrxObs[i], buff.ρ_A_evolved)
        @inbounds buff.dev[i] = real(tr(buff.evolobs[i]))- meas0[i]
        c += buff.dev[i]^2
    end 
    return c
end

#=
function cost!(F::Union{AbstractFloat, nothing}, g::Vector{<:AbstractFloat}, G::Union{Vector{<:AbstractFloat}, nothing}, init::Init)
    @unpack observables, T_max = init.set 
    H_A = @inbounds sum(g[i].*init.blks.matrices[i] for i in eachindex(g))
    if isnothing(G)
        return tanh_sinh(t -> integrand_onlycost(t, H_A, init),0., T_max, init.q)/(length(observables)*T_max)
    end 
    init.buff.C_G_result .= (tanh_sinh(t -> integrand(t, H_A, init),0., T_max, init.q)/(length(observables)*T_max))
    G[:] .= @view  init.buff.C_G_result[2:end]
    return  init.buff.C_G_result[1]
end 
=#

#=
function deviation(j::Int64, ρ_A_copy::DensityMatrix,  init::Init)::Float64
    @unpack observables, meas0 = init.set
    return (expect(observables[j], ρ_A_copy) - meas0[j])
end 

function cost_forpb_scalarintegration(g::Vector{<:AbstractFloat}, G::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, buff::Buffers)
    @unpack observables, T_max, q = set 
    H_A = @inbounds sum(g[i]*blks.matrices[i] for i in eachindex(g))
    fill!(G, 0)
    return tanh_sinh(t -> integrand_forpb(G, set, t, H_A, blks, te),0., T_max, q)[1]/(length(observables)*T_max)
end 

function integrand_forpb_scalar_integration(G::Vector{<:AbstractFloat}, set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64}, blks::H_A_Var, buff::Buffers)
    @unpack N_A, ρ_A, observables, meas0, mtrxObs, T_max = set
    d = 2^N_A 
    U, pullback = rrule(exp, -1im*t*H_A)
    mul!(buff.ρ_A_right, ρ_A.state, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)

    for i in eachindex(mtrxObs)
        @inbounds mul!(buff.evolobs[i], mtrxObs[i], buff.ρ_A_evolved)
        @inbounds buff.dev[i] = real(tr(buff.evolobs[i]))- meas0[i]
        @inbounds buff.sumobs += mtrxObs[i] * buff.dev[i] 
    end 
    
    adj = pullback(adjoint(ρ_A.state*U'*buff.sumobs))[2]

    for i in eachindex(G)
        @inbounds G[i] += 4*t*imag(tr(adj'*blks.matrices[i]))/(T_max * length(observables))
    end 
    fill!(buff.sumobs, 0)
    return sum(buff.dev.^2)
end


function integrand_forZygote(set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64})
    @unpack ρ_A, observables, meas0, mtrxObs = set

    U = exp(-1im*t*H_A)
    ρ_A_evolved = U*ρ_A.state*U'
    im_res = @inbounds sum((tr(mtrxObs[j]*ρ_A_evolved) - meas0[j])^2 for j in eachindex(observables))
    return real(im_res)
end


function cost_forZygote(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    @unpack observables, T_max, q = set 
    H_A = @inbounds sum(g[i]*blks.matrices[i] for i in eachindex(g))
    return tanh_sinh(t -> integrand_forZygote(set, t, H_A),0., T_max, q)[1]/(length(observables)*T_max)
end 

function cost_freeT(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    @unpack observables, q = set 
    H_A = @inbounds sum(g[i]*blks.matrices[i-1] for i in 2:lastindex(g))
    inte =  tanh_sinh(t -> integrand(set, t, H_A),0., g[1], q)[1]/(length(observables)*g[1])
    return inte
end 


function costgrad(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    @unpack observables, T_max, q = set 

    G = zeros(2^N_A, 2^N_A)

    H_A = @inbounds sum(g[i]*blks.matrices[i] for i in eachindex(g))
    return tanh_sinh(t -> integrand(set, t, H_A),0., T_max, q)[1]/(length(observables)*T_max)
end 



function cost_with_midpointrule(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, te::tes)
    @unpack N_A, ρ_A, observables, meas0, T_max, mtrxObs = set
    dt = 0.005
    N_T = T_max ÷ dt
    C::Float64 = 0.0
    H_A = @inbounds sum(g[i]*blks.matrices[i] for i in eachindex(g))
    copyto!(buff.Uhalve, exp(-1im*dt/2*H_A))
    mul!(buff.U, buff.Uhalve, buff.Uhalve)
    mul!(buff.ρ_A_right, ρ_A.state, buff.Uhalve')
    mul!(buff.ρ_A_evolved, buff.Uhalve, buff.ρ_A_right)
    C += real(sum((tr(mtrxObs[j]*buff.ρ_A_evolved) - meas0[j])^2 for j in eachindex(observables)))
    for _ in 1:N_T-1
        mul!(buff.ρ_A_right, buff.ρ_A_evolved, buff.U')
        mul!(buff.ρ_A_evolved, buff.U, buff.ρ_A_right)
        for i in 1:4
            mul!(buff.evolobs[i], mtrxObs[i], buff.ρ_A_evolved)
        end
        C += real(sum((tr(buff.evolobs[j]) - meas0[j])^2 for j in eachindex(observables)))
    end
    return C/(T_max*length(observables))*dt
end



function mul_integrand(set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64}, buff::Buffers)
    @unpack N_A, ρ_A, observables, meas0, T_max, mtrxObs = set
    @time copyto!(buff.U, exp(-1im*t*H_A))
    copyto!(buff.U, exp(-1im*t*H_A))
    mul!(buff.ρ_A_right, ρ_A.state, buff.U')
    mul!(buff.ρ_A_evolved, buff.U, buff.ρ_A_right)
    for i in 1:4
        mul!(buff.evolobs[i], mtrxObs[i], buff.ρ_A_evolved)
    end
    im_res = @inbounds sum((tr(buff.evolobs[j]) - meas0[j])^2 for j in eachindex(observables))
    return real(im_res)
end


function mul_cost(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, buff::Buffers)
    @unpack observables, T_max, q = set 
    H_A = @inbounds sum(g[i]*blks.matrices[i] for i in eachindex(g))
    return tanh_sinh(t -> mul_integrand(set, t, H_A, buff),0., T_max, q)[1]/(length(observables)*T_max)
end 

using BenchmarkTools
function inegrand_costgrad(G::Vector{<:AbstractFloat}, set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64}, blks::H_A_Var, buff::Buffers)
    @unpack N_A, mtrxObs, ρ_A, meas0 = set
    d = 2^N_A 
    copyto!(buff.U, cis(-t*H_A))
    #@btime cis(-$t*Hermitian($H_A))
    mul!(buff.ρ_A_right, ρ_A.state, buff.U')
    mul!(buff.ρ_A_evolved, buff.U, buff.ρ_A_right)

    for i in eachindex(mtrxObs)
        @inbounds mul!(buff.evolobs[i], mtrxObs[i], buff.ρ_A_evolved)
        @inbounds buff.dev[i] = real(tr(buff.evolobs[i]))- meas0[i]
        @inbounds buff.sumobs += mtrxObs[i] * buff.dev[i] 
    end 
    mul!(buff.Minput, buff.ρ_A_right, buff.sumobs)
    copyto!(buff.M, [H_A buff.Minput; 0I H_A])
    #@btime cis(-$t*Hermitian($buff.M))
    copyto!(buff.expm, cis(-t*buff.M))
    copyto!(buff.divU, @view buff.expm[1:d, end-d+1:end])

    for i in eachindex(G)
        @inbounds G[i] = 4*real(tr(buff.divU*blks.matrices[i])) 
    end 
    fill!(buff.sumobs, 0)
        ##dev[i] =  tr(buff.evolobs[j]) - meas0[j]
        #mul!(buff.multdivrho[i], blks.matrices[i], buff.ρ_A_right)
        #mul!(buff.divevolvedrho[i], buff.divU[i], buff.multdivrho[i])
        #mul!(buff.evolobs[i], mtrxObs[i], buff.divevolvedrho[i])
        #G[i] = 4*real(tr(buff.evolobs[i])) 

end 

function mul_cost(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, buff::Buffers)
    @unpack observables, T_max, q = set 
    H_A = @inbounds sum(g[i]*blks.matrices[i] for i in eachindex(g))
    return tanh_sinh(t -> mul_integrand(set, t, H_A, buff),0., T_max, q)[1]/(length(observables)*T_max)
end 

=#
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
    @unpack N_A, ρ_A, observables, meas0, T_max= set
    ρ_A_copy= copy(ρ_A)
    blks.U.dt = t
    c = 0.0
    ρ_A_copy |> blks.U
    for j in eachindex(observables)
        c+= deviation(j, ρ_A_copy, set)^2
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
    @unpack dt, N_T, N_A, ρ_A, observables, meas0, T_max = set
    C::Float64 = 0.0
    ρ_A_copy::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}} = copy(ρ_A)
    H_A::Add{2} = get_HA(g, set)
    U_halve = time_evolve(cache(H_A), dt/2, check_hermicity=false)
    U = time_evolve(cache(H_A), dt, check_hermicity=false)
    Yao.apply!(ρ_A_copy, U_halve)
    for j in eachindex(observables)
        C += deviation(j, ρ_A_copy, set)^2
    end 
    for _ in 1:N_T-1
        Yao.apply!(ρ_A_copy, U)
        for j in eachindex(observables)
            C += deviation(j, ρ_A_copy, set)^2
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