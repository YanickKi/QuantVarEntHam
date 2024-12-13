include("src/QuantVarEntHam.jl")
using .QuantVarEntHam


function whole()
    init = initialize(XXZ(10, 5, 1., 1.,dt = 0.5), H_A_BW)
    g = [1., 2., 3., 4., 5.]
    G = rand(length(g))
    QuantVarEntHam.optimize_midpoint(g, init, gtol=  1e-14, maxiter= 200)
end 

whole()