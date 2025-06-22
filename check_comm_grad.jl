using QuantVarEntHam
using BenchmarkTools

function main()
    tfim = TFIM(12, 6, 1)
    cost1 = QCFL(tfim, H_A_not_BW(tfim), 1)
    println(length(cost1.blocks))
    
    ginit = Float64.([1,1,2,2,2,3,3,4,4,5,5])
    println(length(ginit))
    c = cost1(ginit)
    println(c)
    g, c = QuantVarEntHam.optimize_free(cost1, ginit)
    ##println(g/g[1])
    ##g = Float64.([1, 1,2,2,2,3,3])
    ##g, c = optimize_LBFGS(ginit, QCFL(set, 2, H_A_not_BW, integration_method = midpoint()))
    #g, c = optimize_LBFGS(ginit, commutator(set, H_A_not_BW))
end 

main()