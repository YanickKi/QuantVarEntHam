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

function (c::QCFL)(g::Vector{<:Real})
    get_H_A!(c.buff.qcfl_buff.H_A, g, c.blocks)
    return c.integrator.scalar_integrate(t -> integrand_cost(c, g, t), c.T_max)/(length(c.observables)*c.T_max)
end 

function gradient!(G, c::QCFL, g::Vector{<:Real})
    get_H_A!(c.buff.qcfl_buff.H_A, g, c.blocks)
    c.integrator.vector_integrate(c.buff.qcfl_buff.C_G_result, t -> integrand_gradient(c, g, t), c.T_max)
    c.buff.qcfl_buff.C_G_result ./= length(c.observables)*c.T_max
    G[:] .= @view c.buff.qcfl_buff.C_G_result[2:end]
end


function integrand_gradient(c::QCFL, g::Vector{<:Real}, t::AbstractFloat)
    
    ρ_A = c.model.ρ_A 
    observables = c.observables
    blocks = c.blocks 
    meas0 = c.meas0 
    buff = c.buff

    fix_first_index = false 

    get_H_A!(buff.qcfl_buff.H_A, g, blocks)

    qcfl_buff = buff.qcfl_buff

    qcfl_buff.H_A_forexp .=  -1im*t .* qcfl_buff.H_A 
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

    @fastmath @inbounds @simd for i in 1+fix_first_index:lastindex(qcfl_buff.C_G)-1+fix_first_index
        qcfl_buff.C_G[i+1-fix_first_index] = gradient_component(buff, blocks[i], t)
    end    
    
    return qcfl_buff.C_G
end

function integrand_cost(c::QCFL, g::Vector{<:Real}, t::AbstractFloat)
    
    ρ_A = c.model.ρ_A 
    observables = c.observables
    blocks = c.blocks 
    meas0 = c.meas0 
    buff = c.buff

    get_H_A!(buff.qcfl_buff.H_A, g, blocks)
    qcfl_buff = buff.qcfl_buff

    qcfl_buff.H_A_forexp .=  -1im*t .* qcfl_buff.H_A 
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