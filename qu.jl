using QuantVarEntHam
using BenchmarkTools

function main()
    init = initialize(TFIM(8, 4, 1, 2), H_A_not_BW)
    println(init)
    #QuantVarEntHam.cost_count!(g, init)
end 

main()