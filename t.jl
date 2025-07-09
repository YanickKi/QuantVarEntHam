using QuantVarEntHam
using Cthulhu
using BenchmarkTools
function main()
    model = TFIM(8, 4, 1)
    blocks = H_A_BW(model)
    cost = QCFL(model, blocks, 1)
    @descend cost([1,2,3,4])
    #@descend optimize(cost, [1,2,3,4])
end 

main()