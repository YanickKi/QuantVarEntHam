mutable struct Relative_entropy{M} <: AbstractCostFunction
    model::M
    blocks::Vector{Matrix{ComplexF64}}
    buff::Relative_entropy_buffer
    trace_exp_H_A::Float64 # cached variable to save intermediate calculations in gradient for cost
end 

function Relative_entropy(model::AbstractModel, blocks::Vector{<:AbstractMatrix})
    return Relative_entropy(
        model, 
        blocks, 
        make_relative_entropy_buffer(model),
        0.
    )
end 

function (c::Relative_entropy)(g)
    buff = c.buff
    ρ_A = c.model.ρ_A

    get_H_A!(c, g)
    buff.H_A_forexp  .= -1 .* buff.H_A
    exp_only_buffered!(buff.H_A_forexp, buff.exp_buff)
    exp_H_A = buff.exp_buff.X

    S1 = tr(mul!(buff.ρ_A_times_H_A, ρ_A, buff.H_A))
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 


function _gradient!(c::Relative_entropy, G, g::Vector{<:Real}, free_indices)
    
    get_H_A!(c, g)

    pre_computations_gradient(c)

    @fastmath @inbounds @simd for i in eachindex(free_indices)
        G[i] = gradient_component(c, free_indices[i])
    end

    return G 
end 

function pre_computations_gradient(c::Relative_entropy)
    
    H_A_forexp = c.buff.H_A_forexp
    exp_buff = c.buff.exp_buff
   
    H_A_forexp .= -1 .* c.buff.H_A

    exp_only_buffered!(H_A_forexp, exp_buff)
    exp_H_A = exp_buff.X
    
    c.trace_exp_H_A = tr(exp_H_A)

end 


function gradient_component(c::Relative_entropy, index::Integer)

    ρ_A_times_block = c.buff.ρ_A_times_H_A  # renaming for readability 
    exp_H_A = c.buff.exp_buff.X             # renaming for readability 

    γ = - 1/c.trace_exp_H_A

    ∂S1 = tr(mul!(ρ_A_times_block, c.model.ρ_A, c.blocks[index]))
    ∂S2 = γ * tr(mul!(ρ_A_times_block, exp_H_A, c.blocks[index]))
    return ∂S1+∂S2
end 


function cost_from_gradient_intermediates(c::Relative_entropy)
    S1 = tr(mul!(c.buff.ρ_A_times_H_A, c.model.ρ_A, c.buff.H_A))
    S2 = log(c.trace_exp_H_A)
    return S1+S2
end 

cost_from_gradient_intermediates(fc::FixedCost{C}) where {C<:Relative_entropy} = unwrap(cost_from_gradient_intermediates, fc)


function fg!(F, G::Union{Vector{<:Real}, Nothing}, c::Union{C, FixedCost{C}}, g::Vector{<:Real}) where {C<:Relative_entropy}
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
