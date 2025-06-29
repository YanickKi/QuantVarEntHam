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

    shorten_buffers!(c,length(fixed_indices))

    return FixedCost(
        c, 
        fixed_indices,
        Float64.(fixed_values),
        free_indices,
        full_g
    ) 
end 


shorten_buffers!(::AbstractCostFunction, ::Integer) = nothing 

function fill_full_g(fc::FixedCost, g::Vector{Float64})
   fc.full_g[fc.free_indices] .= g
   return fc.full_g
end 


function (fc::FixedCost)(g::Vector{<:Real})
    fill_full_g(fc, g)
    return fc.c(fc.full_g)
end 

_gradient!(fc::FixedCost, G::Vector{<:Real}, g::Vector{<:Real}, free_indices) = unwrap(_gradient!, fc, G, fill_full_g(fc, g), free_indices)
