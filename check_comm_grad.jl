using QuantVarEntHam
using BenchmarkTools

function main()
    tfim = TFIM(8, 4, 1)
    cost1 = Relative_entropy(tfim, H_A_not_BW(tfim))
    fixed_cost = FixedCost(cost1, [1], [1.5592665717881316])
    ginit = Float64.([1,1,2,2,2,3,3])
    ginit_fixed = Float64.([1,2,2,2,3,3])
 
    #g, c = QuantVarEntHam.optimize_free(cost1, ginit)
    g, c = QuantVarEntHam.optimize_free(fixed_cost, ginit_fixed)
    ##println(g/g[1])
    ##g = Float64.([1, 1,2,2,2,3,3])
    ##g, c = optimize_LBFGS(ginit, QCFL(set, 2, H_A_not_BW, integration_method = midpoint()))
    #g, c = optimize_LBFGS(ginit, commutator(set, H_A_not_BW))
end 

main()