abstract type AbstractCostFunction end 


struct FixedCost{C} <: AbstractCostFunction
    c::C 
    fixed_indices::Vector{Int}
    fixed_values::Vector{Float64}
    free_indices::Vector{Int} 
    full_g::Vector{Float64}
end 

function FixedCost(c::AbstractCostFunction, fixed_indices::Vector{<:Integer}, fixed_values::Vector{<:Real})
    @assert typeof(c) != FixedCost "Cannot wrap a fixed cost into a new fixed cost!"
    c = deepcopy(c)
    full_g = zeros(length(c.blocks))
    
    full_g[fixed_indices] .= fixed_values
    free_indices = setdiff(eachindex(c.blocks), fixed_indices)
    return FixedCost(
        c, 
        fixed_indices,
        Float64.(fixed_values),
        free_indices,
        full_g
    ) 
end 


function get_H_A!(c::AbstractCostFunction, g::Vector{<:AbstractFloat})
    c.buff.H_A .= g[1].* c.blocks[1]

    @fastmath @inbounds @simd for i in 2:length(g)
        c.buff.H_A .+= g[i].* c.blocks[i]
    end 
    return c.buff.H_A 
end 

function get_H_A!(fc::FixedCost, g::Vector{<:AbstractFloat})

    fill_full_g(fc, g)

    get_H_A!(fc.c, fc.full_g)
end 

function (fc::FixedCost)(g::Vector{<:Real})
    fill_full_g(fc, g)
    return fc.c(fc.full_g)
end 

function fg!(F, G::Union{Vector{<:Real}, Nothing}, c::AbstractCostFunction, g::Vector{<:Real})
    if !isnothing(G)
        C = gradient!(G, c, g)
        return C 
    end 
    if !isnothing(F)
        return c(g)
    end 
end 

function fill_full_g(fc::FixedCost, g::Vector{Float64})
   fc.full_g[fc.free_indices] .= g
end 

include("QCFL.jl")
include("relative_entropy.jl")
include("commutator.jl")