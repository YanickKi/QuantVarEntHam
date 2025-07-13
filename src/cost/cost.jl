export gradient!
export AbstractCostFunction, AbstractFreeCostFunction
export QCFL, Commutator, RelativeEntropy, FixedCost
export QCFLBuffer, CommutatorBuffer, RelativeEntropyBuffer
export fill_full_g
"""
    AbstractCostFunction

Abstract type for cost functions including where parameters are fixed.
"""
abstract type AbstractCostFunction end 


"""
    AbstractFreeCostFunction

Abstract type for cost functions only where all parameters are free.
"""
abstract type AbstractFreeCostFunction <: AbstractCostFunction end



"""
    gradient!(G::Vector{<:Real}, c::AbstractCostFunction, g::Vector{<:Real})

Computes the gradient of a given cost function `c` at point `g` and saves it in `G`.
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
include("printing.jl")