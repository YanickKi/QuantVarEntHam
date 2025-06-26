using QuantVarEntHam
using BenchmarkTools
function main()
    model = TFIM(20, 10, 1)
    cost1 = Relative_entropy(model, H_A_notBW(model))
    fixed_cost = FixedCost(cost1, [1], [1])
    ginit = Float64.([1,1,2,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9])
    ginit_fixed = Float64.([1,2,2,3,3])
    g, c = optimize(cost1, ginit)
    #g, c = QuantVarEntHam.optimize_free(fixed_cost, ginit_fixed)
    ##println(g/g[1])
    ##g = Float64.([1, 1,2,2,2,3,3])
    ##g, c = optimize_LBFGS(ginit, QCFL(set, 2, H_A_not_BW, integration_method = midpoint()))
    #g, c = optimize_LBFGS(ginit, commutator(set, H_A_not_BW))
end 

main()