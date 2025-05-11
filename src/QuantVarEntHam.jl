module QuantVarEntHam

export H_A_BW, H_A_not_BW
export optimize_LBFGS
export get_rhoA, initialize, XXZ, TFIM
export universal_ratios, print_H_A, plot_universal_ratios
export Settings, Settings_XXZ, Settings_TFIM
export H_XXZ, H_TFIM, H_pollmann
export Init, H_A_Var
export tanh_sinh 
export midpoint
export pollmann
using Parameters
using StaticArrays

include("integration/tanh-sinh.jl")
include("integration/mapsum.jl")
include("integration/tanh-sinh_vector.jl")
include("integration/mapsum_vector.jl")
include("integration/midpoint.jl")
include("hamiltonians/hamiltonians.jl")
include("frechet_exp/frechet_exp.jl")
include("initialize.jl")
include("cost.jl")
include("optimizer.jl")
include("utils.jl")
end
