"""
    FixedCost{C<:AbstractFreeCostFunction} <: AbstractCostFunction
    FixedCost(cost::AbstractFreeCostFunction, fixed_indices::Vector{<:Integer}, fixed_values::Vector{<:Real})
    FixedCost(cost::AbstractFreeCostFunction, fixed_indices::Int, fixed_values::Real)
    FixedCost(cost::AbstractFreeCostFunction, fixed_indices::NTuple{L,Int}, fixed_values::NTuple{L,Int}) where {L}

Wrapper for fixing parameters of a cost function object. 

It stores a given cost function without fixed parameters `c` of the type [`AbstractFreeCostFunction`](@ref) together with the indices (`fixed_indices`) and the values (`fixed_values`) of the fixed parameters.
"""
struct FixedCost{C<:AbstractFreeCostFunction} <: AbstractCostFunction
    c::C
    fixed_indices::Vector{Int}
    fixed_values::Vector{Float64}
    free_indices::Vector{Int}
    full_g::Vector{Float64}
end

function FixedCost(
    cost::AbstractFreeCostFunction,
    fixed_indices::Vector{<:Integer},
    fixed_values::Vector{<:Real},
)
    cost = deepcopy(cost)

    full_g = zeros(length(cost.ansatz.blocks))

    full_g[fixed_indices] .= fixed_values
    free_indices = setdiff(eachindex(cost.ansatz.blocks), fixed_indices)

    shorten_buffers!(cost, length(fixed_indices))

    return FixedCost(cost, fixed_indices, Float64.(fixed_values), free_indices, full_g)
end

function FixedCost(cost::AbstractFreeCostFunction, fixed_indices::Int, fixed_values::Real)
    FixedCost(cost, [fixed_indices], [fixed_values])
end

function FixedCost(
    cost::AbstractFreeCostFunction, fixed_indices::NTuple{L,Int}, fixed_values::NTuple{L,Int}
) where {L}
    FixedCost(cost, [fixed_indices...], [fixed_values...])
end

shorten_buffers!(::AbstractFreeCostFunction, ::Integer) = nothing

"""
    fill_full_g(fc::FixedCost, g::Vector{<:Real})

Return the complete parameter vector consisting of the fixed parameter in `fc` and the passed parameters `g`.
"""
fill_full_g(fc::FixedCost, g::Vector{<:Real}) = fill_full_g(fc, Float64.(g))

function fill_full_g(fc::FixedCost, g::Vector{Float64})
    fc.full_g[fc.free_indices] .= g
    return fc.full_g
end

function (fc::FixedCost)(g::Vector{<:Real})
    fill_full_g(fc, g)
    return fc.c(fc.full_g)
end

function _gradient!(fc::FixedCost, G::Vector{<:Real}, g::Vector{<:Real}, free_indices)
    unwrap(_gradient!, fc, G, fill_full_g(fc, g), free_indices)
end
