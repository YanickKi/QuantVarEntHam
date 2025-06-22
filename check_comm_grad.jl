using QuantVarEntHam
using BenchmarkTools

function main()
    set = TFIM(8, 4, 1)
    ginit = Float64.([1, 1,2,2,2,3,3])
    #g, c = optimize_LBFGS(g, QCFL(set, 1, H_A_not_BW))
    #println(g/g[1])
    #g = Float64.([1, 1,2,2,2,3,3])
    #g, c = optimize_LBFGS(ginit, QCFL(set, 2, H_A_not_BW, integration_method = midpoint()))
    g, c = optimize_LBFGS(ginit, commutator(set, H_A_not_BW))


end 

main()