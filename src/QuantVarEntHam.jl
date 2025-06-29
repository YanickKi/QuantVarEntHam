module QuantVarEntHam

using Parameters, LinearAlgebra


include("integration/integration.jl")
include("hamiltonians/hamiltonians.jl")
include("frechet_exp/frechet_exp.jl")
include("cost/cost.jl")
include("interface_cost_integration.jl")
include("optimizer.jl")
#include("utils.jl")
end
