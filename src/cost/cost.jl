abstract type AbstractCostFunction end 

function get_H_A!(H_A::AbstractMatrix, g::Vector{<:AbstractFloat}, blocks::Vector{<:AbstractMatrix})
    H_A .= g[1].*blocks[1]

    @fastmath @inbounds @simd for i in 2:length(g)
        H_A .+= g[i].*blocks[i]
    end 
    return H_A 
end 

include("QCFL.jl")
include("relative_entropy.jl")
include("commutator.jl")