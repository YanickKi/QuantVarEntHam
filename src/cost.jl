using LinearAlgebra
using QuadGK


function get_H_A!(g::Vector{<:AbstractFloat}, init::Init)
    init.buff.H_A .= g[1].*init.blocks[1]

    @fastmath @inbounds @simd for i in 2:length(g)
        init.buff.H_A .+= g[i].*init.blocks[i]
    end 
end 

function QCFL(;integration_method = tanh_sinh(), g1::Real=NaN)
    return (F, G, g, init) -> fg!(F, G, g, init, integration_method, fix_first_index = !isnan(g1)), Float64(g1)
end 

function fg!(F::Union{AbstractFloat, Nothing}, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, init::Init, integration_method; fix_first_index::Bool = false)
    @unpack observables, T_max = init.set
    get_H_A!(g, init)
    init.buff.H_A .*= -1im
    method_scalar, method_vector, need_buffer = integration_method
    if G !== nothing
        if need_buffer == true
            method_vector(init.buff.C_G_result, t -> integrand_gradient(t, init, fix_first_index), T_max, init.buff.Σ)
        else 
            method_vector(init.buff.C_G_result, t -> integrand_gradient(t, init, fix_first_index), T_max)
        end 
        init.buff.C_G_result ./= length(observables)*T_max
        G[:] .= @view init.buff.C_G_result[2:end]
    end
    if F !== nothing 
        if G === nothing
            In = method_scalar(t -> integrand_cost(t, init), T_max)/(length(observables)*T_max)
            return In
        else 
            return init.buff.C_G_result[1]
        end 
    end
end 


function integrand_gradient(t::AbstractFloat, init::Init, fix_first_index::Bool)
    @unpack ρ_A, meas0, observables = init.set
    buff = init.buff
    
    buff.H_A_forexp .=  t .* buff.H_A 
    pullback =  own_rrule(exp, buff.H_A_forexp, init.exp_buf)
    
    evolve(init)  
    
    δ = dev(init, 1)
    buff.sumobs .= δ.*observables[1]
    buff.C_G[1] = δ^2

    @fastmath @inbounds @simd for i in 2:lastindex(observables)
        δ = dev(init, i)
        buff.sumobs .+= δ.*observables[i]
        buff.C_G[1] += δ^2
    end 

    mul!(buff.dAforpb, buff.ρ_A_right, buff.sumobs)

    pullback(buff.dAforpb)

    @fastmath @inbounds @simd for i in 1+fix_first_index:lastindex(buff.C_G)-1+fix_first_index
        buff.C_G[i+1-fix_first_index] = gradient_component(init, i, t)
    end    
    
    return buff.C_G
end

function integrand_cost(t::AbstractFloat,  init::Init)
    @unpack ρ_A, meas0, observables = init.set 

    buff = init.buff

    buff.H_A_forexp .=  t .* buff.H_A 
    exp_bufered!(buff.H_A_forexp, init.exp_buf)
    evolve(init)

    c = 0.
    @fastmath @inbounds @simd for i in eachindex(observables)
        c += dev(init, i)^2
    end 
    return c
end

function evolve(init::Init)
    @unpack ρ_A = init.set
    buff = init.buff
    U = init.exp_buf.X
    mul!(buff.ρ_A_right, ρ_A, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)
end 

function dev(init::Init, index_observable::Integer)
    @unpack meas0, observables = init.set
    @inbounds mul!(init.buff.evolob, observables[index_observable], init.buff.ρ_A_evolved)
    @inbounds δ = real(tr(init.buff.evolob))- meas0[index_observable]
    return δ
end 

function gradient_component(init::Init, index_parameter, t::Real)
    buff = init.buff
    frech = init.exp_buf.∂X
    mul!(buff.frechTimesBlock, frech, init.blocks[index_parameter])
    return 4*t*imag(tr(buff.frechTimesBlock))
end


#########################################################################################
#                                                                                       #
#                   Commutator cost function                                            #
#                                                                                       #
#########################################################################################

function commutator()
    return (F, G, g, init) -> comm_fg!(F, G, g, init), NaN 
end 

function comm_fg!(F::Union{AbstractFloat, Nothing}, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, init::Init)
 
    get_H_A!(g, init)
    C = 0.
    if G !== nothing
        C = comm_grad!(G, eachindex(g), init)
    end
    if F !== nothing 
        if G === nothing
            C = comm_cost(init)
            return C
        else 
            return C
        end 
    end
end 


function comm_cost(init::Init)

    comm = init.buff.dAforpb
    H_A = init.buff.H_A
    ρ_A = init.set.ρ_A

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1)
    
    return norm(comm)/(2*norm(H_A)*norm(ρ_A))
end 

function comm_grad!(G::Vector{<:AbstractFloat}, indices::AbstractUnitRange{Int64}, init::Init)

    comm = init.buff.dAforpb
    H_A = init.buff.H_A
    ρ_A = init.set.ρ_A    
    temp = init.buff.H_A_forexp

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1)


    Λ_comm = norm(comm)
    Λ_H_A = norm(H_A)
    Λ_ρ_A = norm(ρ_A)

    for index in indices
        mul!(mul!(temp, ρ_A, init.blocks[index]), init.blocks[index], ρ_A, 1, -1)
        ∂Λ_comm = 1/Λ_comm * real(sum(had!(temp, temp, comm)))
        ∂Λ_H_A = 1/Λ_H_A * real(sum(had!(temp, H_A, init.blocks[index])))
        G[index] = 1/(2*Λ_ρ_A*Λ_H_A^2)*(∂Λ_comm * Λ_H_A - Λ_comm * ∂Λ_H_A) 
    end 
    return Λ_comm/(2*Λ_H_A * Λ_ρ_A) 
end 

function had!(buf::AbstractMatrix, A::AbstractMatrix,B::AbstractMatrix)
    n = size(A)[1]
    for j in 1:n
       for i in 1:n
         @inbounds buf[i,j] = A[i,j] *B[i,j]
       end
    end
    return buf
end


#########################################################################################
#                                                                                       #
#                   Relative entropy                                                    #
#                                                                                       #
#########################################################################################

function relative_entropy()
    return (F, G, g, init) -> relative_entropy_fg!(F, G, g, init), NaN
end 

function relative_entropy_fg!(F::Union{AbstractFloat, Nothing}, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, init::Init)
 
    get_H_A!(g, init)
    C = 0
    if G !== nothing
        C = relative_entropy_grad!(G, length(g), init)
    end
    if F !== nothing 
        if G === nothing
            C = relative_entropy_cost(init)
            return C
        else 
            return C
        end 
    end
end 

function relative_entropy_cost(init::Init)
    @unpack ρ_A  = init.set
    buff = init.buff
    
    buff.H_A_forexp .= -1 .*init.buff.H_A
    exp_bufered!(buff.H_A_forexp, init.exp_buf)
    exp_H_A = init.exp_buf.X

    S1 = tr(mul!(buff.dAforpb, ρ_A, buff.H_A))
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 

function relative_entropy_grad!(G::Union{AbstractVector, Nothing}, num_params::Integer, init::Init)
    @unpack ρ_A  = init.set
    buff = init.buff
    buff.H_A_forexp .= -1 .*buff.H_A

    exp_bufered!(buff.H_A_forexp, init.exp_buf)
    exp_H_A = init.exp_buf.X
    γ = -1/(tr(exp_H_A))

    for i in 1:num_params
        ∂S1 = tr(mul!(buff.dAforpb, ρ_A, init.blocks[i]))
        ∂S2 = γ * tr(mul!(buff.dAforpb, exp_H_A, init.blocks[i]))
        G[i] = ∂S1+∂S2
    end

    S1 = tr(mul!(buff.dAforpb, ρ_A, buff.H_A))
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 

#########################################################################################
#                                                                                       #
#                   INTEGRAND AFTER TANH-SINH TRANSFORMATION                            #
#                                                                                       #
#########################################################################################

function tanh_sinh_integrand(g::Vector{<:AbstractFloat}, u1::Real, u2::Real, num_points::Int, init::Init)
    @unpack T_max, observables = init.set
    u1_ = Float64(u1)
    u2_ = Float64(u2)
    u = range(start = u1_, stop = u2_, length = num_points)
    sinhu = sinh.(u)
    Φ = tanh.(sinhu .* π/2)
    ϕ′ = (cosh.(u) .*π/2) ./ (cosh.(sinhu .* π/2)).^2
    c = Vector{Float64}(undef, length(u))
    get_H_A!(g, init)
    for i in eachindex(u) 
        c[i] = integrand_onlycost(T_max/2*(Φ[i]+1), init) * ϕ′[i]
    end 

    return c.*T_max/(2*length(observables))
end 



#########################################################################################
#                                                                                       #
#                                COUNT INTEGRAND EVALUTAIONS                            #
#                                                                                       #
#########################################################################################

function cost_count(g::Vector{<:AbstractFloat}, init::Init; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Int = 12, h0::Float64 = 1.)
    @unpack observables, T_max = init.set
    q = integration_tables(maxlevel = maxlevel, h0 = h0)
    count = [0]
    get_H_A!(g, init)
    init.buff.H_A .*= -1im 
    I = tanh_sinh(t -> integrand_onlycost_count(t, init, count),0., T_max, q, atol = atol, rtol = rtol)/(length(observables)*T_max)
    return I, count[1]
end 

function integrand_count(t::AbstractFloat,  init::Init, count::Vector{Int})
    @unpack ρ_A, meas0, observables = init.set
    count[1] += 1
    return integrand_onlycost(t, init)
end


#########################################################################################
#                                                                                       #
#                                QUADGK                                                 #
#                                                                                       #
#########################################################################################

function quadgk_count!(g::Vector{<:AbstractFloat}, init::Init)
    @unpack observables, T_max, atol, rtol = init.set
    get_H_A!(g, init)
    init.buff.H_A .*= -1im
    I, E, count =  quadgk_count(t -> integrand_onlycost(t, init),0., T_max)
    return I/(length(observables)*T_max), count 
end 


#=

#########################################################
#                                                       #
#           TORIC                                       #
#                                                       #
#########################################################

function cg_toric(F::Union{AbstractFloat, Nothing}, g::Vector{<:Real}, G::Union{Vector{<:AbstractFloat}, Nothing}, init::Init, integration_method)
    @unpack observables, T_max = init.set
    get_H_A!(g, init)

    decomp = eigen(init.buff.H_A)

    method_scalar, method_vector = integration_method
    if G !== nothing
        method_vector(t -> integrand(t, init, decomp), T_max, init.buff.C_G_result, init.buff.Σ)
        init.buff.C_G_result ./= length(observables)*T_max
        G[:] .= @view init.buff.C_G_result[2:end]
    end
    if F !== nothing 
        if G === nothing
            return method_scalar(t -> integrand_onlycost(t, init, decomp), T_max)/(length(observables)*T_max)
        else 
            return init.buff.C_G_result[1]
        end 
    end
end 

function integrand(t::AbstractFloat, init::Init, decomp)
    @unpack ρ_A, meas0, observables = init.set
    buff = init.buff
    buff.C_G[1] = 0.

    U = decomp.vectors * Diagonal(cis.(-t * decomp.values)) * decomp.vectors'
    mul!(buff.ρ_A_right, ρ_A, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)

    @fastmath @inbounds @simd for i in eachindex(observables)
        mul!(buff.evolobs[i], observables[i], buff.ρ_A_evolved)
        buff.dev[i] = real(tr(buff.evolobs[i])) - meas0[i]
        buff.sumobs .+= buff.dev[i].*observables[i]
        buff.C_G[1] += buff.dev[i]^2
    end 


    mul!(buff.dAforpb, buff.ρ_A_evolved, buff.sumobs) 

    @fastmath @inbounds @simd for i in eachindex(buff.G_buffer)
        mul!(buff.frechTimesBlock, buff.dAforpb, init.blocks[i])
    end     

    @fastmath @inbounds @simd for i in eachindex(buff.G_buffer)
        buff.G_buffer[i] = 4*t*imag(tr(buff.frechTimesBlock))
    end 

    fill!(buff.sumobs, 0)
    buff.C_G[2:end] .= @view buff.G_buffer[:]
    return buff.C_G
end

function integrand_onlycost(t::AbstractFloat,  init::Init, decomp)
    @unpack ρ_A, meas0, observables = init.set
    buff = init.buff 
    mul!(buff.H_A_forexp, -t, buff.H_A)
    U = decomp.vectors * Diagonal(cis.(-t * decomp.values)) * decomp.vectors'
    mul!(buff.ρ_A_right, ρ_A, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)
    c = 0.
    @fastmath @inbounds @simd for i in eachindex(observables)
        @inbounds mul!(buff.evolobs[i], observables[i], buff.ρ_A_evolved)
        @inbounds buff.dev[i] = real(tr(buff.evolobs[i]))- meas0[i]
        c += buff.dev[i]^2
    end 
    return c
end


=#


#=

                SAFE WITH RRULE FROM CHAINRULES

function cost_grad!(F::Union{AbstractFloat, Nothing}, g::Vector{<:Real}, G::Union{Vector{<:AbstractFloat}, Nothing}, init::Init, integration_method)
    @unpack observables, T_max = init.set
    get_H_A!(g, init)
    method_scalar, method_vector = integration_method
    if G !== nothing
        method_vector(t -> integrand(t, init), T_max, init.buff.C_G_result, init.buff.Σ)
        init.buff.C_G_result ./= length(observables)*T_max
        G[:] .= @view init.buff.C_G_result[2:end]
    end
    if F !== nothing 
        if G === nothing
            return method_scalar(t -> integrand_onlycost(t, init), T_max)/(length(observables)*T_max)
        else 
            return init.buff.C_G_result[1]
        end 
    end
end 


function integrand(t::AbstractFloat, init::Init)
    @unpack ρ_A, meas0, observables = init.set
    buff = init.buff
    buff.C_G[1] = 0.
    mul!(buff.H_A_forexp, -1im*t, buff.H_A)
    U, pullback =  rrule(exp, buff.H_A_forexp)
    mul!(buff.ρ_A_right, ρ_A, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)

    @fastmath @inbounds @simd for i in eachindex(observables)
        mul!(buff.evolobs[i], observables[i], buff.ρ_A_evolved)
        buff.dev[i] = real(tr(buff.evolobs[i])) - meas0[i]
        buff.sumobs .+= buff.dev[i].*observables[i]
        buff.C_G[1] += buff.dev[i]^2
    end 

    mul!(buff.dAforpb, buff.ρ_A_right, buff.sumobs)

    adj = pullback(buff.dAforpb')[2]
    
    @fastmath @inbounds @simd for i in eachindex(buff.G_buffer)
        mul!(buff.frechTimesBlock, adj', init.blocks[i])
    end     

    @fastmath @inbounds @simd for i in eachindex(buff.G_buffer)
        buff.G_buffer[i] = 4*t*imag(tr(buff.frechTimesBlock))
    end 

    fill!(buff.sumobs, 0)
    buff.C_G[2:end] .= @view buff.G_buffer[:]
    return buff.C_G
end
=#

#=

#########################################################################################
#                                                                                       #
#                                VARYING TMAX                                           #
#                                                                                       #
#########################################################################################

function cost_T_max(g::Vector{<:AbstractFloat}, init::Init)
    @unpack observables, atol, rtol = init.set
    get_H_A!(vcat(g[2:end]), init)

    return tanh_sinh(t -> integrand_onlycost(t, init),0., g[1], init.q, atol = atol, rtol = rtol)/(length(observables)*g[1])
end 


function cost_grad_Tmax_g1_fixed!(F::Union{AbstractFloat, Nothing}, g::Vector{<:AbstractFloat}, G::Union{Vector{<:AbstractFloat}, Nothing}, init::Init)
    @unpack observables, atol, rtol = init.set
    get_H_A!(g[2:end], init)
    if G !== nothing
        init.buff.C_G_result[1:end-1] .= (tanh_sinh(t -> integrand_Tmax_g1_fixed(t, init),0., g[1], init.q, atol = atol, rtol = rtol)./(length(observables)*g[1]))
        G[2:end] .= @view  init.buff.C_G_result[2:end-1]
        G[1] = 1/g[1] * (integrand_onlycost(g[1], init)/length(observables) - init.buff.C_G_result[1])
    end
    if F !== nothing 
        if G === nothing
            return tanh_sinh(t -> integrand_onlycost(t, init),0., g[1], init.q, atol = atol, rtol = rtol)/(length(observables)*g[1])
        else 
            return init.buff.C_G_result[1]
        end 
    end
end 

function integrand_Tmax_g1_fixed(t::AbstractFloat, init::Init)
    @unpack ρ_A, meas0, observables = init.set
    
    numBlocks = length(init.blocks)

    buff = init.buff
    buff.C_G[1] = 0.
    mul!(buff.H_A_forexp, -1im*t, buff.H_A)
    U, pullback =  rrule(exp, buff.H_A_forexp)
    mul!(buff.ρ_A_right, ρ_A, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)

    @fastmath @inbounds @simd for i in eachindex(observables)
        mul!(buff.evolobs[i], observables[i], buff.ρ_A_evolved)
        buff.dev[i] = real(tr(buff.evolobs[i])) - meas0[i]
        buff.sumobs .+= buff.dev[i].*observables[i]
        buff.C_G[1] += buff.dev[i]^2
    end 

    mul!(buff.dAforpb, buff.ρ_A_right, buff.sumobs)

    adj = pullback(buff.dAforpb')[2]
    
    @fastmath @inbounds @simd for i in 2:numBlocks
        mul!(buff.frechTimesBlock, adj', init.blocks[i])
    end     

    @fastmath @inbounds @simd for i in 2:numBlocks
        buff.G_buffer[i] = 4*t*imag(tr(buff.frechTimesBlock))
    end 

    fill!(buff.sumobs, 0)
    buff.C_G[2:end-1] .= @view buff.G_buffer[2:end]
    return buff.C_G[1:end-1]
end


=#
#=
function cost!(F::Union{AbstractFloat, nothing}, g::Vector{<:AbstractFloat}, G::Union{Vector{<:AbstractFloat}, nothing}, init::Init)
    @unpack observables, T_max = init.set 
    H_A = @inbounds sum(g[i].*init.blocks[i] for i in eachindex(g))
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
    H_A = @inbounds sum(g[i]*blocks[i] for i in eachindex(g))
    fill!(G, 0)
    return tanh_sinh(t -> integrand_forpb(G, set, t, H_A, blks, te),0., T_max, q)[1]/(length(observables)*T_max)
end 

function integrand_forpb_scalar_integration(G::Vector{<:AbstractFloat}, set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64}, blks::H_A_Var, buff::Buffers)
    @unpack N_A, ρ_A, observables, meas0, observables, T_max = set
    d = 2^N_A 
    U, pullback = rrule(exp, -1im*t*H_A)
    mul!(buff.ρ_A_right, ρ_A, U')
    mul!(buff.ρ_A_evolved, U, buff.ρ_A_right)

    for i in eachindex(observables)
        @inbounds mul!(buff.evolobs[i], observables[i], buff.ρ_A_evolved)
        @inbounds buff.dev[i] = real(tr(buff.evolobs[i]))- meas0[i]
        @inbounds buff.sumobs += observables[i] * buff.dev[i] 
    end 
    
    adj = pullback(adjoint(ρ_A*U'*buff.sumobs))[2]

    for i in eachindex(G)
        @inbounds G[i] += 4*t*imag(tr(adj'*blocks[i]))/(T_max * length(observables))
    end 
    fill!(buff.sumobs, 0)
    return sum(buff.dev.^2)
end


function integrand_forZygote(set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64})
    @unpack ρ_A, observables, meas0, observables = set

    U = exp(-1im*t*H_A)
    ρ_A_evolved = U*ρ_A*U'
    im_res = @inbounds sum((tr(observables[j]*ρ_A_evolved) - meas0[j])^2 for j in eachindex(observables))
    return real(im_res)
end


function cost_forZygote(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    @unpack observables, T_max, q = set 
    H_A = @inbounds sum(g[i]*blocks[i] for i in eachindex(g))
    return tanh_sinh(t -> integrand_forZygote(set, t, H_A),0., T_max, q)[1]/(length(observables)*T_max)
end 

function cost_freeT(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    @unpack observables, q = set 
    H_A = @inbounds sum(g[i]*blocks[i-1] for i in 2:lastindex(g))
    inte =  tanh_sinh(t -> integrand(set, t, H_A),0., g[1], q)[1]/(length(observables)*g[1])
    return inte
end 


function costgrad(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var)
    @unpack observables, T_max, q = set 

    G = zeros(2^N_A, 2^N_A)

    H_A = @inbounds sum(g[i]*blocks[i] for i in eachindex(g))
    return tanh_sinh(t -> integrand(set, t, H_A),0., T_max, q)[1]/(length(observables)*T_max)
end 



function cost_with_midpointrule(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, te::tes)
    @unpack N_A, ρ_A, observables, meas0, T_max, observables = set
    dt = 0.005
    N_T = T_max ÷ dt
    C::Float64 = 0.0
    H_A = @inbounds sum(g[i]*blocks[i] for i in eachindex(g))
    copyto!(buff.Uhalve, exp(-1im*dt/2*H_A))
    mul!(buff.U, buff.Uhalve, buff.Uhalve)
    mul!(buff.ρ_A_right, ρ_A, buff.Uhalve')
    mul!(buff.ρ_A_evolved, buff.Uhalve, buff.ρ_A_right)
    C += real(sum((tr(observables[j]*buff.ρ_A_evolved) - meas0[j])^2 for j in eachindex(observables)))
    for _ in 1:N_T-1
        mul!(buff.ρ_A_right, buff.ρ_A_evolved, buff.U')
        mul!(buff.ρ_A_evolved, buff.U, buff.ρ_A_right)
        for i in 1:4
            mul!(buff.evolobs[i], observables[i], buff.ρ_A_evolved)
        end
        C += real(sum((tr(buff.evolobs[j]) - meas0[j])^2 for j in eachindex(observables)))
    end
    return C/(T_max*length(observables))*dt
end



function mul_integrand(set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64}, buff::Buffers)
    @unpack N_A, ρ_A, observables, meas0, T_max, observables = set
    @time copyto!(buff.U, exp(-1im*t*H_A))
    copyto!(buff.U, exp(-1im*t*H_A))
    mul!(buff.ρ_A_right, ρ_A, buff.U')
    mul!(buff.ρ_A_evolved, buff.U, buff.ρ_A_right)
    for i in 1:4
        mul!(buff.evolobs[i], observables[i], buff.ρ_A_evolved)
    end
    im_res = @inbounds sum((tr(buff.evolobs[j]) - meas0[j])^2 for j in eachindex(observables))
    return real(im_res)
end


function mul_cost(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, buff::Buffers)
    @unpack observables, T_max, q = set 
    H_A = @inbounds sum(g[i]*blocks[i] for i in eachindex(g))
    return tanh_sinh(t -> mul_integrand(set, t, H_A, buff),0., T_max, q)[1]/(length(observables)*T_max)
end 

using BenchmarkTools
function inegrand_costgrad(G::Vector{<:AbstractFloat}, set::Settings, t::AbstractFloat, H_A::Matrix{ComplexF64}, blks::H_A_Var, buff::Buffers)
    @unpack N_A, observables, ρ_A, meas0 = set
    d = 2^N_A 
    copyto!(buff.U, cis(-t*H_A))
    #@btime cis(-$t*Hermitian($H_A))
    mul!(buff.ρ_A_right, ρ_A, buff.U')
    mul!(buff.ρ_A_evolved, buff.U, buff.ρ_A_right)

    for i in eachindex(observables)
        @inbounds mul!(buff.evolobs[i], observables[i], buff.ρ_A_evolved)
        @inbounds buff.dev[i] = real(tr(buff.evolobs[i]))- meas0[i]
        @inbounds buff.sumobs += observables[i] * buff.dev[i] 
    end 
    mul!(buff.Minput, buff.ρ_A_right, buff.sumobs)
    copyto!(buff.M, [H_A buff.Minput; 0I H_A])
    #@btime cis(-$t*Hermitian($buff.M))
    copyto!(buff.expm, cis(-t*buff.M))
    copyto!(buff.divU, @view buff.expm[1:d, end-d+1:end])

    for i in eachindex(G)
        @inbounds G[i] = 4*real(tr(buff.divU*blocks[i])) 
    end 
    fill!(buff.sumobs, 0)
        ##dev[i] =  tr(buff.evolobs[j]) - meas0[j]
        #mul!(buff.multdivrho[i], blocks[i], buff.ρ_A_right)
        #mul!(buff.divevolvedrho[i], buff.divU[i], buff.multdivrho[i])
        #mul!(buff.evolobs[i], observables[i], buff.divevolvedrho[i])
        #G[i] = 4*real(tr(buff.evolobs[i])) 

end 

function mul_cost(g::Vector{<:AbstractFloat}, set::Settings, blks::H_A_Var, buff::Buffers)
    @unpack observables, T_max, q = set 
    H_A = @inbounds sum(g[i]*blocks[i] for i in eachindex(g))
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
    H_A = @inbounds sum(g[i]*blocks[i] for i in eachindex(g))
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
    H_A = sum(g[i]*blocks[i] for i in eachindex(g))
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
    H_A = sum(g[i]*blocks[i] for i in eachindex(g))
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