using QuantVarEntHam
using BenchmarkTools
function main()
    model = TFIM(8, 4, 1)
    cost1 = QCFL(model, H_A_notBW(model),1, integrator = MidPoint(1e-2))
    fixed_cost = FixedCost(cost1, [1], [1])
    ginit = Float64.([1,1,2,2,2,3,3])
    ginit_fixed = Float64.([1,2,2,2,3,3])
    g, c = optimize(cost1, ginit)
    g, c = QuantVarEntHam.optimize(fixed_cost, ginit_fixed)
    ##println(g/g[1])
    ##g = Float64.([1, 1,2,2,2,3,3])
    ##g, c = optimize_LBFGS(ginit, QCFL(set, 2, H_A_not_BW, integration_method = midpoint()))
    #g, c = optimize_LBFGS(ginit, commutator(set, H_A_not_BW))
end 

main()