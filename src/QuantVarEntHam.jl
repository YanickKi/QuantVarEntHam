module QuantVarEntHam

export Settings_XXZ, Settings_TFIM
export H_A_BW, H_A_not_BW
export utilBlocks
export optimize_LBFGS
export get_rhoA

using Yao, Parameters

include("hamiltonians/hamiltonians.jl")
include("cost.jl")
include("optimizer.jl")

end


#= todo 
-ThreadedDoubleExponentialFormula is always generating the same tables for the same Time interval. Since the table of evaluation points does not change throughout one optimization 
procedure, one could evaluate the table at the beginning of the simulation thus saving computation time 
-ThreadedDoubleExponentialFormula has some bugs for multithreading, sometimes the integral is not correctly evaluated (observed for higher amounts of threads ~ 8)
-some functions are not type stable in module QuantVarEntHam
=#