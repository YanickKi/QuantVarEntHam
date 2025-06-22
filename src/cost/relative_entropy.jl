struct Relative_entropy{M} <: AbstractCostFunction
    model::M
    blocks::Vector{Matrix{ComplexF64}}
    buff::Relative_entropy_buffer
end 

function Relative_entropy(model::AbstractModel, blocks::Vector{<:AbstractMatrix})
    return Relative_entropy(
        model, 
        blocks, 
        make_relative_entropy_buffer(model)
    )
end 

function (c::Relative_entropy)(g)
    buff = c.buff
    ρ_A = c.model.ρ_A

    get_H_A!(buff.H_A, g, c.blocks)
    buff.H_A_forexp  .= -1 .* buff.H_A
    exp_only_buffered!(buff.H_A_forexp, buff.exp_buff)
    exp_H_A = buff.exp_buff.X

    S1 = tr(mul!(buff.ρ_A_times_H_A, ρ_A, buff.H_A))
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 

function gradient!(G, c::Relative_entropy, g::Vector{<:Real})
    
    H_A = c.buff.H_A
    H_A_forexp = c.buff.H_A_forexp
    ρ_A_times_H_A = c.buff.ρ_A_times_H_A
    exp_buff = c.buff.exp_buff
    ρ_A = c.model.ρ_A

    get_H_A!(H_A, g, c.blocks)
   
    H_A_forexp .= -1 .*H_A


    ρ_A_times_block = ρ_A_times_H_A # just a renaming so it's more readable

    exp_only_buffered!(H_A_forexp, exp_buff)
    exp_H_A = exp_buff.X
    γ = -1/(tr(exp_H_A))

    @fastmath @inbounds @simd for i in eachindex(G)
        ∂S1 = tr(mul!(ρ_A_times_block, ρ_A, c.blocks[i]))
        ∂S2 = γ * tr(mul!(ρ_A_times_block, exp_H_A, c.blocks[i]))
        G[i] = ∂S1+∂S2
    end

    S1 = tr(mul!(ρ_A_times_block, ρ_A, H_A))
    S2 = log(tr(exp_H_A))
end 