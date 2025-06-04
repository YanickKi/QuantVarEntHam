module QuantVarEntHam

export Settings, Settings_XXZ, Settings_TFIM, Settings_pollmann
export H_XXZ, H_TFIM, H_pollmann
export XXZ, TFIM, pollmann
export H_A_BW, H_A_not_BW
export Init, initialize, get_rhoA
export tanh_sinh, midpoint
export optimize_LBFGS
export universal_ratios

using Parameters

include("integration/integration.jl")
include("hamiltonians/hamiltonians.jl")
include("frechet_exp/frechet_exp.jl")
include("initialize.jl")
include("cost.jl")
include("optimizer.jl")
include("utils.jl")
end
