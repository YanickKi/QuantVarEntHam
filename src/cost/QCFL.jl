struct QCFL{M, I, O} <: AbstractCostFunction
    model::M 
    blocks::Vector{Matrix{ComplexF64}}
    integrator::I
    observables::Vector{O} 
    T_max::Float64
    meas0::Vector{Float64}
    buff::QCFL_buffer
end
 
function QCFL(model::AbstractModel, blocks::Vector{<:AbstractMatrix}, T_max::Real; integrator::Union{Nothing,AbstractIntegrator} = nothing, observables::Union{Nothing, Vector{<:AbstractMatrix}} = nothing)
    
    observables = something(observables, [repeat(model.N_A, Z, (i,i+1), S = model.S) for i in 1:model.N_A-1])
    integrator = something(integrator, Tanh_sinh(length(blocks)+1))
    meas0 = [expect(observable, model.ρ_A) for observable in observables]
    return QCFL(
        model, 
        blocks, 
        integrator,
        observables,
        Float64(T_max),
        meas0,
        make_QCFL_buffer(model, length(blocks), observables)
    )
end 

function shorten_buffers(c::QCFL, fixed_indices::AbstractVector)
    if typeof(c.integrator) == Tanh_sinh
        for _ in eachindex(fixed_indices)
            pop!(c.integrator.buffer)
        end 
    end 
    for _ in eachindex(fixed_indices)
        pop!(c.buff.qcfl_buff.C_G)
        pop!(c.buff.qcfl_buff.C_G_result)
    end 

end 

function (c::QCFL)(g::Vector{<:Real})
    get_H_A!(c, g)
    return c.integrator.scalar_integrate(t -> integrand_cost(c, t), c.T_max)/(length(c.observables)*c.T_max)
end 

function _gradient!(c::QCFL, G, g::Vector{<:Real}, free_indices)
    get_H_A!(c, g)
    c.integrator.vector_integrate(c.buff.qcfl_buff.C_G_result, t -> integrand_gradient(c, t, free_indices), c.T_max)
    c.buff.qcfl_buff.C_G_result ./= length(c.observables)*c.T_max
    G[:] .= @view c.buff.qcfl_buff.C_G_result[2:end]
    return c.buff.qcfl_buff.C_G_result[1]
end


function pre_computations_gradient(c::QCFL, t::AbstractFloat)
    ρ_A = c.model.ρ_A 
    observables = c.observables
    meas0 = c.meas0 
    buff = c.buff

    qcfl_buff = buff.qcfl_buff

    qcfl_buff.H_A_forexp .=  -1im*t .* buff.H_A 
    pullback =  own_rrule(exp, qcfl_buff.H_A_forexp, buff.exp_buff)
    
    evolve(ρ_A, buff.exp_buff.X, qcfl_buff)  
    
    δ = dev(observables[1], meas0[1], qcfl_buff)
    qcfl_buff.sumobs .= δ.*observables[1]
    qcfl_buff.C_G[1] = δ^2

    @fastmath @inbounds @simd for i in 2:lastindex(observables)
        δ = dev(observables[i], meas0[i], qcfl_buff)
        qcfl_buff.sumobs .+= δ.*observables[i]
        qcfl_buff.C_G[1] += δ^2
    end 

    pullback(mul!(qcfl_buff.dAforpb, qcfl_buff.ρ_A_right, qcfl_buff.sumobs))

end 

function integrand_gradient(c::QCFL, t::AbstractFloat, free_indices)

    pre_computations_gradient(c, t)

    @fastmath @inbounds @simd for i in eachindex(free_indices)
        c.buff.qcfl_buff.C_G[i+1] = gradient_component(c.buff, c.blocks[free_indices[i]], t)
    end    
    
    return c.buff.qcfl_buff.C_G
end


function integrand_cost(c::QCFL, t::AbstractFloat)
    
    ρ_A = c.model.ρ_A 
    observables = c.observables
    meas0 = c.meas0 
    buff = c.buff

    qcfl_buff = buff.qcfl_buff

    qcfl_buff.H_A_forexp .=  -1im*t .* buff.H_A 
    exp_only_buffered!(qcfl_buff.H_A_forexp, buff.exp_buff)

    evolve(ρ_A, buff.exp_buff.X, qcfl_buff)  

    c = 0.
    @fastmath @inbounds @simd for i in eachindex(observables)
        c += dev(observables[i], meas0[i], qcfl_buff)^2
    end 
    return c
end

function evolve(ρ_A::AbstractMatrix, U::AbstractMatrix, qcfl_buff)
    mul!(qcfl_buff.ρ_A_evolved, U, mul!(qcfl_buff.ρ_A_right, ρ_A, U'))
end 

function dev(observable::AbstractMatrix, expect_at_0::AbstractFloat, qcfl_buff)
    @inbounds mul!(qcfl_buff.evolob, observable, qcfl_buff.ρ_A_evolved)
    @inbounds δ = real(tr(qcfl_buff.evolob))- expect_at_0
    return δ
end 

function gradient_component(buff, block::AbstractMatrix, t::Real)
    qcfl_buff = buff.qcfl_buff
    frech = buff.exp_buff.∂X
    mul!(qcfl_buff.frechTimesBlock, frech, block)
    return 4*t*imag(tr(qcfl_buff.frechTimesBlock))
end