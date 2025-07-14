unwrap(f::Function, fc::FixedCost, args...) = f(fc.c, args...)

get_free_indices(fc::FixedCost) = fc.free_indices
get_free_indices(c::AbstractFreeCostFunction) = eachindex(c.blocks)


function fg!(F, G::Union{Vector{<:Real}, Nothing}, c::AbstractCostFunction, g::Vector{<:Real})
    if !isnothing(G)
        C = gradient!(G, c, g)
        return C 
    end 
    if !isnothing(F)
        return c(g)
    end 
end 


function get_H_A!(c::AbstractCostFunction, g::Vector{<:AbstractFloat})
    c.buff.H_A .= g[1] .* c.blocks[1]

    @fastmath @inbounds @simd for i in 2:length(g)
        c.buff.H_A .+= g[i] .* c.blocks[i]
    end 
    return c.buff.H_A 
end 

get_H_A!(fc::FixedCost, g::Vector{<:AbstractFloat}) = get_H_A!(fc.c, fill_full_g(fc, g))