"""
    FixedCost{C <: AbstractFreeCostFunction} <: AbstractCostFunction

Wrapper for fixing parameters of a cost function object. 

# Fields 

- `c::C`: cost function 
- `fixed_indices::Vector{Int}`: the indices of the parameters which are to be fixed 
- `fixed_values::Vector{Float64}`: values of the parameters which are to be fixed 
- `free_indices::Vector{Int}`: free indices which are to be optimized 
- `full_g::Vector{Float64}`: a buffer to construct the full parameter vector with fixed and free parameters
"""
struct FixedCost{C <: AbstractFreeCostFunction} <: AbstractCostFunction
    c::C 
    fixed_indices::Vector{Int}
    fixed_values::Vector{Float64}
    free_indices::Vector{Int} 
    full_g::Vector{Float64}
end 


"""
    FixedCost(c::AbstractFreeCostFunction, fixed_indices::Vector{<:Integer}, fixed_values::Vector{<:Real})

Outer constructor for [`FixedCost`](@ref).

# Arguments

- `c`: cost function 
- `fixed_indices`: the indices of the parameters which are to be fixed 
- `fixed_values`: values of the parameters which are to be fixed 
"""
function FixedCost(c::AbstractFreeCostFunction, fixed_indices::Vector{<:Integer}, fixed_values::Vector{<:Real})
        
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


shorten_buffers!(::AbstractFreeCostFunction, ::Integer) = nothing 

"""
    fill_full_g(fc::FixedCost, g::Vector{Float64})

Return the complete parameter vector consisting of the fixed parameter in `fc` and the passed parameters `g`.
"""
function fill_full_g(fc::FixedCost, g::Vector{Float64})
   fc.full_g[fc.free_indices] .= g
   return fc.full_g
end 

function (fc::FixedCost)(g::Vector{<:Real})
    fill_full_g(fc, g)
    return fc.c(fc.full_g)
end 

_gradient!(fc::FixedCost, G::Vector{<:Real}, g::Vector{<:Real}, free_indices) = unwrap(_gradient!, fc, G, fill_full_g(fc, g), free_indices)
