using QuantVarEntHam
using Cthulhu
using BenchmarkTools
function main()
    model = TFIM(8, 4, 1)
    cost1 = QCFL(model, H_A_notBW(model),1, integrator = MidPoint(1e-2))
    g = Float64.([1,1,2,2,2,3,3])
    fixed_cost = FixedCost(cost1, [1], [1])
    Float64.([1,1,2,2,2,3,3])
    gfixed = Float64.([1,2,2,2,3,3])
    @descend fixed_cost(gfixed)
end 

main()