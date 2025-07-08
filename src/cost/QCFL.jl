struct QCFLOnlyBuffer{O<:AbstractMatrix}
    C_G::Vector{Float64}
    C_G_result::Vector{Float64}
    ρ_A_evolved::Matrix{ComplexF64}
    ρ_A_right::Matrix{ComplexF64}
    direction_frech::Matrix{ComplexF64}
    H_A_forexp::Matrix{ComplexF64}
    sumobs::O
end 

function QCFLOnlyBuffer(d::Integer, numBlocks::Integer, observables::Vector{<:AbstractMatrix}) 
    return QCFLOnlyBuffer(
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks+1),
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d),
    zeros(ComplexF64, d, d),
    zero(observables[1]),
    )
end


"""
    QCFLBuffer{O<:AbstractMatrix}

Struct containing the buffers for the [`QCFL`](@ref).
"""
struct QCFLBuffer{O<:AbstractMatrix}
    qcfl_buff::QCFLOnlyBuffer{O}
    exp_buff::ExpFrechBuffer{ComplexF64}
    H_A::Matrix{ComplexF64}
end


"""
    QCFLBuffer(model::AbstractModel, blocks::Vector{<:AbstractMatrix}, observables::Vector{<:AbstractMatrix})

Outer constructor for [`QCFLBuffer`](@ref) given a `model`, `blocks` and `observables`.
"""
function QCFLBuffer(model::AbstractModel, blocks::Vector{<:AbstractMatrix}, observables::Vector{<:AbstractMatrix})
    numBlocks = length(blocks) 
    d = size(model.ρ_A)[1] # Hilbert space dimension
    return QCFLBuffer(
        QCFLOnlyBuffer(d, numBlocks, observables), 
        ExpFrechBuffer(ComplexF64, d), 
        zeros(ComplexF64, d, d) 
    )
end 

"""
    QCFL{M<:AbstractModel, O<:AbstractMatrix} <: AbstractFreeCostFunction

The cost function from the QCFL as an object defined as 

```math
\\mathcal{C}(\\vec{g}) = \\frac{1}{T_\\text{max}} \\int_0^{T_\\text{max}} \\frac{1}{N_\\text{O}} \\sum_{j=1}^{N_\\text{O}} \\left ( \\langle \\mathcal{O}_j^\\text{A} \\rangle_t - \\langle \\mathcal{O}_j^\\text{A} \\rangle_0 \\right ) dt
```
with the expectation value of the time evolved observables
```math 
\\langle \\mathcal{O}_j^\\text{A} \\rangle_t =  \\text{Tr}_{\\text{A}} \\left [ \\mathcal{O}_j^\\text{A} e^{- i H_\\text{A}^\\text{Var}(\\vec{g})t} \\rho_\\text{A}  e^{i H_\\text{A}^\\text{Var}(\\vec{g})t} \\right ].
```
``N_\\text{O}`` denotes the number of observables.
# Fields 
- `model::M`: model 
- `blocks::Vector{Matrix{ComplexF64}}`: blocks of the variational Ansatz
- `integrator::Integrator`: object containing scalar and vector integration (NOTE: INTERNAL ONLY)
- `observables::Vector{O}`: monitored observables
- `T_max::Float64`: maximum integration time aka how long the system is sampled 
- `meas0::Vector{Float64}`: expectation values at time `t=0` 
- `buff::QCFLBuffer{O}`: see [`QCFLBuffer`](@ref)

# Gradient 

The gradient is given by 

```math
    \\partial_{g_k} \\mathcal{C}(\\vec{g}) = \\frac{4}{T_\\text{max} N_\\text{O}} \\int_0^{T_\\text{max}} t \\, \\text{Im} \\left \\{ \\text{Tr}_{\\text{A}} 
    \\left [ \\mathcal{L}_{e^X} ( - i t  H_\\text{A}^\\text{Var}, \\rho_\\text{A} U_\\text{A}^\\dagger \\Xi_\\text{A}) h_\\text{k} \\right ] \\right   \\} dt ,
```
where ``\\Xi_\\text{A} = \\sum_j \\langle \\mathcal{O}_j^\\text{A} \\rangle_t - \\langle \\mathcal{O}_j^\\text{A} \\rangle_0 \\mathcal{O}_j^\\text{A}``
and ``\\mathcal{L}_{e^X} ( - i t H_\\text{A}^\\text{Var}, \\rho_\\text{A} U_\\text{A}^\\dagger \\Xi_\\text{A})``
denotes the Frechet derivative of the matrix exponential at ``- i t H_\\text{A}^\\text{Var}``
in the direction of ``\\rho_\\text{A} U_\\text{A}^\\dagger \\Xi_\\text{A}``.
"""
struct QCFL{M<:AbstractModel, O<:AbstractMatrix, S<:AbstractScalarIntegrator, V<:AbstractVectorIntegrator} <: AbstractFreeCostFunction
    model::M 
    blocks::Vector{Matrix{ComplexF64}}
    integrator::Integrator{S,V}
    observables::Vector{O} 
    T_max::Float64
    meas0::Vector{Float64}
    buff::QCFLBuffer{O}
end
 

"""
    QCFL(model::AbstractModel, blocks::Vector{<:AbstractMatrix}, T_max::Real; integrator::Union{Nothing,AbstractIntegrator} = nothing, observables::Union{Nothing, Vector{<:AbstractMatrix}} = nothing,
    buffer::Union{Nothing, QCFLBuffer} = nothing) 

Outer constructor for [`QCFL`](@ref) s.t. the correct buffers will be automatically constructed for 
the given `model`, `blocks` and `observables`.

Provides default values for observables and the integrator.
The default values for observables is `` \\{ Z_i Z_{i+1} | 1 \\leq i < N_\\text{A} \\}``.
The tanh-sinh quadrature is set as the default integrator with its default values (see [`TanhSinh(; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Integer=12, h0::Real=1)`](@ref)).

# Required arguments

- `model`: model 
- `blocks`: blocks for the variational Ansatz
- `T_max`: maximum integration time aka how long the system is sampled

# Keyword arguments
- `integrator`: settings for integration method (see [`AbstractIntegrator`](@ref))
- `observables`: observables
- `buffer`: see [`QCFLBuffer`](@ref) 
"""
function QCFL(model::AbstractModel, blocks::Vector{<:AbstractMatrix}, T_max::Real; integrator::Union{Nothing,AbstractIntegrator} = nothing, observables::Union{Nothing, Vector{<:AbstractMatrix}} = nothing,
    buffer::Union{Nothing, QCFLBuffer} = nothing)
    
    observables = something(observables, [Diagonal(Matrix(repeat(model.N_A, Z, (i,i+1), S = model.S))) for i in 1:model.N_A-1])
    integrator = something(integrator, TanhSinh())
    buffer = something(buffer, QCFLBuffer(model, blocks, observables))
    meas0 = [expect(observable, model.ρ_A) for observable in observables]


    complete_integrator = make_integrator(blocks, integrator)
    
    return QCFL(
        model, 
        blocks, 
        complete_integrator,
        observables,
        Float64(T_max),
        meas0,
        buffer
    )
end 

function shorten_buffers!(c::QCFL, how_often::Integer)
    shorten_buffer!(buffertrait(c.integrator.vector_integrate),c.integrator.vector_integrate, how_often)
    for _ in 1:how_often
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
    
    δ = dev(observables[1], meas0[1], qcfl_buff.ρ_A_evolved)
    qcfl_buff.sumobs .= δ.*observables[1]
    qcfl_buff.C_G[1] = δ^2

    @fastmath @inbounds @simd for i in 2:lastindex(observables)
        δ = dev(observables[i], meas0[i], qcfl_buff.ρ_A_evolved)
        qcfl_buff.sumobs .+= δ.*observables[i]
        qcfl_buff.C_G[1] += δ^2
    end 

    pullback(mul!(qcfl_buff.direction_frech, qcfl_buff.ρ_A_right, qcfl_buff.sumobs))

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
        c += dev(observables[i], meas0[i], qcfl_buff.ρ_A_evolved)^2
    end 
    return c
end

function evolve(ρ_A::AbstractMatrix, U::AbstractMatrix, qcfl_buff::QCFLOnlyBuffer)
    mul!(qcfl_buff.ρ_A_evolved, U, mul!(qcfl_buff.ρ_A_right, ρ_A, U'))
end 

function dev(observable::AbstractMatrix, expect_at_0::AbstractFloat, ρ_A_evolved::AbstractMatrix)
    return  real(dot(observable, ρ_A_evolved)) - expect_at_0
end 

function gradient_component(buff, block::AbstractMatrix, t::Real)
    frech = buff.exp_buff.∂X
    return 4*t*imag(dot(block,frech))
end