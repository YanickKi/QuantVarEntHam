using QuantVarEntHam
using BenchmarkTools

function main()
    init = initialize(TFIM(8, 4, 1, 2), H_A_not_BW)
    g = Float64.([1, 1,2,2,2,3,3])
    g, c = QuantVarEntHam.optimize_LBFGS(g, init, cost = commutator())
    println(g)
end 

main()