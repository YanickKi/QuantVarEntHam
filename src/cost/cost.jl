abstract type AbstractCostFunction end 

export gradient!
export AbstractCostFunction, QCFL, Commutator, Relative_entropy, FixedCost
export QCFL_buffer, Commutator_buffer, Relative_entropy_buffer

"""
    gradient!

Computes the gradient of a given cost function `c` at point `g` and saves it in `G`.

# Arguments
-`G::Vector{<:Real}`: buffer to save gradient in 
-`c::AbstractCostFunction`: cost function
-`g::Vector{<:Real}`: parameters

"""
function gradient!(G::Vector{<:Real}, c::AbstractCostFunction, g::Vector{<:Real})
    
    free_indices = get_free_indices(c)
    
    return _gradient!(c,G,g, free_indices) 
end 


include("FixedCost.jl")
include("QCFL.jl")
include("relative_entropy.jl")
include("commutator.jl")
include("utils.jl")