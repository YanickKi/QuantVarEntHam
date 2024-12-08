include("src/QuantVarEntHam.jl")

using .QuantVarEntHam
using Profile
#using PProf
function whole()
    init = initialize(TFIM(8, 5, 1., 1.), H_A_BW)
    g = [1., 2., 3., 4., 5.]
    G = rand(length(g))
    QuantVarEntHam.cost_grad!(1.,g, G, init)    
    Profile.Allocs.clear()
    @time Profile.Allocs.@profile sample_rate = 1 QuantVarEntHam.cost_grad!(1.,g, G, init)
    #PProf.Allocs.pprof(from_c = false)
end 
whole() 