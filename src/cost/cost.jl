export AbstractCostFunction, AbstractFreeCostFunction
export QCFL, Commutator, RelativeEntropy, FixedCost
export QCFLBuffer, CommutatorBuffer, RelativeEntropyBuffer
export fill_full_g
export gradient
export print_blocks, print_observables, print_model
export getblocks, getobservables, getmodel

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
    gradient(cost::AbstractCostFunction, g::Vector{<:Real})

Return the gradient of a given `cost` at `g`.
"""
function gradient(cost::AbstractCostFunction, g::Vector{<:Real})
    
    G = similar(g)

    free_indices = get_free_indices(cost)
    
    _gradient!(cost,G,g, free_indices) 
    
    return G 
end

function gradient!(G::Vector{<:Real}, cost::AbstractCostFunction, g::Vector{<:Real})
    
    free_indices = get_free_indices(cost)
    
    return _gradient!(cost,G,g, free_indices) 
end 

include("FixedCost.jl")
include("QCFL.jl")
include("relative_entropy.jl")
include("commutator.jl")
include("utils.jl")
include("printing.jl")
include("getter.jl")