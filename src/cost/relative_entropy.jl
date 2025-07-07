"""
    RelativeEntropyBuffer

Struct containing the buffers for the [`RelativeEntropy`](@ref).
"""
struct RelativeEntropyBuffer
    H_A::Matrix{ComplexF64}
    H_A_forexp::Matrix{ComplexF64}
    exp_buff::ExpBuffer{ComplexF64}
end 

"""
    RelativeEntropyBuffer(model::AbstractModel)

Outer constructor for [`RelativeEntropyBuffer`](@ref) given a `model`. 
"""
function RelativeEntropyBuffer(model::AbstractModel)

    d = size(model.ρ_A)[1]
    return RelativeEntropyBuffer(zeros(ComplexF64, d, d), zeros(ComplexF64, d, d), 
        ExpBuffer(ComplexF64, d))
end 

"""
    RelativeEntropy{M<:AbstractModel} <: AbstractCostFunction

Relative entropy as a cost function as an object defined as 
```math 
\\mathcal{C}(\\vec{g}) = \\text{Tr} [ρ_A H_\\text{A}^\\text{Var}(\\vec{g}) ] + \\ln  \\text{Tr}_\\text{A} [e^{-H_\\text{A}^\\text{Var}(\\vec{g})}]
```
# Fields 

- `model::M`: model
- `blocks::Vector{Matrix{ComplexF64}}`: blocks for the variational Ansatz
- `buff::RelativeEntropyBuffer`: buffers (see [`RelativeEntropyBuffer`](@ref)) 
- `trace_exp_H_A::Float64`: cached variable to store an intermiadte result

# Gradient 

The gradient is given by 

```math
\\partial_{g_k} \\mathcal{C}(\\vec{g}) = \\text{Tr}_{\\text{A}} [  \\rho_\\text{A} h_k ] - \\frac{\\text{Tr}_{\\text{A}} [ e^{-H_\\text{A}^\\text{Var}(\\vec{g})} h_k ] }
{\\text{Tr}_{\\text{A}} [ e^{-H_\\text{A}^\\text{Var}(\\vec{g})}]}.
```
"""
mutable struct RelativeEntropy{M<:AbstractModel} <: AbstractCostFunction
    model::M
    blocks::Vector{Matrix{ComplexF64}}
    buff::RelativeEntropyBuffer
    trace_exp_H_A::Float64 # cached variable to save intermediate calculations in gradient for cost
end 

"""
    RelativeEntropy(model::AbstractModel, blocks::Vector{<:AbstractMatrix})

Outer constructor for [`RelativeEntropy`](@ref) s.t. the correct `buffers´ (see [`RelativeEntropyBuffer`](@ref)) will be automatically constructed
for a given `model` and `blocks`.
An already existing `buffer` can be provided.
"""
function RelativeEntropy(model::AbstractModel, blocks::Vector{<:AbstractMatrix})
    buffer = RelativeEntropyBuffer(model)
    return RelativeEntropy(
        model, 
        blocks,
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

    ∂S1 = dot(c.model.ρ_A, c.blocks[index])
    ∂S2 = γ * dot(exp_H_A, c.blocks[index])
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
