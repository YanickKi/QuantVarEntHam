module QuantVarEntHam

using LinearAlgebra


function check_S_N_type(S,N)
    isa(S, Rational{Int}) || throw(ArgumentError("wrong type; S must be of type Rational{Int$(Sys.WORD_SIZE)}"))
    isa(N, Int) || throw(ArgumentError("wrong type; N must be of type Int$(Sys.WORD_SIZE)"))
end 

include("integration/integration.jl")
include("spinoperators/spinoperators.jl")
include("models/models.jl")
include("frechet_exp/frechet_exp.jl")
include("cost/cost.jl")
include("interface_cost_integration.jl")
include("optimizer.jl")
#include("utils.jl")
end
