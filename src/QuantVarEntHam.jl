module QuantVarEntHam

export H_XXZ, H_TFIM, H_pollmann
export XXZ, TFIM, pollmann
export H_A_BW, H_A_not_BW
export get_rhoA
export tanh_sinh, midpoint
export optimize_LBFGS
#export universal_ratios
export QCFL, Commutator, Relative_entropy, tanh_sinh_integrand, cost_count, quadgk_count!
export make_QCFL_buffer, make_commutator_buffer, make_relative_entropy_buffer
using Parameters

include("integration/integration.jl")
include("hamiltonians/hamiltonians.jl")
include("frechet_exp/frechet_exp.jl")
include("buffer.jl")
include("cost/cost.jl")
include("optimizer.jl")
#include("utils.jl")
end
