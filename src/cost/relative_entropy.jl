"""
    RelativeEntropyBuffer{S,N_A}
    RelativeEntropyBuffer(::AbstractModel{S,N_A}) where {S,N_A} 

Buffers for the [`RelativeEntropy`](@ref) for a given [`AbstractModel`](@ref) .
"""
struct RelativeEntropyBuffer{S,N_A}
    H_A::Matrix{ComplexF64}
    H_A_forexp::Matrix{ComplexF64}
    exp_buff::ExpBuffer{ComplexF64}
end 


function RelativeEntropyBuffer(::AbstractModel{S,N_A}) where {S,N_A} 

    d = Int((2*S+1)^N_A)
    return RelativeEntropyBuffer{S,N_A}(zeros(ComplexF64, d, d), zeros(ComplexF64, d, d), 
        ExpBuffer(ComplexF64, d))
end 

"""
    RelativeEntropy{M<:AbstractModel, A<:AbstractAnsatz} <: AbstractFreeCostFunction
    RelativeEntropy(model::AbstractModel{S,N_A}, ansatz::AbstractAnsatz{S,N_A}, buffer::Union{Nothing, RelativeEntropyBuffer{S,N_A}} = nothing) where {S,N_A}

Contains the model, blocks, observables and buffers for the relative entropy as cost function.

Existing buffers can be provided, and are constructed automatically otherwise. 

The cost function is given by 
```math 
\\mathcal{C}(\\vec{g}) = \\text{Tr}_{\\text{A}} [ρ_A H_\\text{A}^\\text{Var}(\\vec{g}) ] + \\ln  \\text{Tr}_\\text{A} [e^{-H_\\text{A}^\\text{Var}(\\vec{g})}]
```

# Gradient 

The gradient is given by 

```math
\\partial_{g_k} \\mathcal{C}(\\vec{g}) = \\text{Tr}_{\\text{A}} [  \\rho_\\text{A} h_k ] - \\frac{\\text{Tr}_{\\text{A}} [ e^{-H_\\text{A}^\\text{Var}(\\vec{g})} h_k ] }
{\\text{Tr}_{\\text{A}} [ e^{-H_\\text{A}^\\text{Var}(\\vec{g})}]}.
```
"""
mutable struct RelativeEntropy{M<:AbstractModel, A<:AbstractAnsatz} <: AbstractFreeCostFunction
    model::M
    blocks_mat::Vector{Matrix{ComplexF64}}
    ansatz::A
    buff::RelativeEntropyBuffer
    trace_exp_H_A::Float64 # cached variable to save intermediate calculations in gradient for cost
end 

function RelativeEntropy(model::AbstractModel{S,N_A}, ansatz::AbstractAnsatz{S,N_A}, buffer::Union{Nothing, RelativeEntropyBuffer{S,N_A}} = nothing) where {S,N_A}
    blocks_mat = Matrix.(mat.(ansatz.blocks))
    buffer = @something buffer RelativeEntropyBuffer(model)
    return RelativeEntropy(
        model, 
        blocks_mat,
        ansatz,
        buffer,
        0.
    )
end 


function (c::RelativeEntropy)(g)
    buff = c.buff
    ρ_A = c.model.ρ_A

    get_H_A!(c, g)
    buff.H_A_forexp  .= -1 .* buff.H_A
    exp_only_buffered!(buff.H_A_forexp, buff.exp_buff)
    exp_H_A = buff.exp_buff.X

    S1 = dot(ρ_A, buff.H_A)
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 


function _gradient!(c::RelativeEntropy, G, g::Vector{<:Real}, free_indices)
    
    get_H_A!(c, g)

    pre_computations_gradient(c)

    @fastmath @inbounds @simd for i in eachindex(free_indices)
        G[i] = gradient_component(c, free_indices[i])
    end

    return G 
end 

function pre_computations_gradient(c::RelativeEntropy)
    
    H_A_forexp = c.buff.H_A_forexp
    exp_buff = c.buff.exp_buff
   
    H_A_forexp .= -1 .* c.buff.H_A

    exp_only_buffered!(H_A_forexp, exp_buff)
    exp_H_A = exp_buff.X
    
    c.trace_exp_H_A = tr(exp_H_A)

end 


function gradient_component(c::RelativeEntropy, index::Integer)

    exp_H_A = c.buff.exp_buff.X # renaming for readability 

    γ = - 1/c.trace_exp_H_A

    ∂S1 = dot(c.model.ρ_A, c.blocks_mat[index])
    ∂S2 = γ * dot(exp_H_A, c.blocks_mat[index])
    return ∂S1+∂S2
end 


function cost_from_gradient_intermediates(c::RelativeEntropy)
    S1 = dot(c.model.ρ_A, c.buff.H_A)
    S2 = log(c.trace_exp_H_A)
    return S1+S2
end 

cost_from_gradient_intermediates(fc::FixedCost{C}) where {C<:RelativeEntropy} = unwrap(cost_from_gradient_intermediates, fc)


function fg!(F, G::Union{Vector{<:Real}, Nothing}, c::Union{C, FixedCost{C}}, g::Vector{<:Real}) where {C<:RelativeEntropy}
    if !isnothing(G)
        gradient!(G, c, g)
        if !isnothing(F)
            return cost_from_gradient_intermediates(c)
        end 
    end 
    if !isnothing(F)
        return c(g)
    end 
end 
