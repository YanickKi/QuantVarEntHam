module QuantVarEntHam

using LinearAlgebra

function divisible_by_half(S)
    iszero(S % 1//2) || throw(ArgumentError("S must be divisible by two!"))
end

include("integration/integration.jl")
include("spinoperators/spinoperators.jl")
include("models/models.jl")
include("frechet_exp/frechet_exp.jl")
include("cost/cost.jl")
include("interface_cost_integration.jl")
include("optimizer.jl")
include("universal_ratios.jl")
end
