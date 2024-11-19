include("./src/QuantVarEntHam.jl")
using .QuantVarEntHam
using BenchmarkTools
using Cthulhu

function whole()
    set = Settings_XXZ(N=10, N_A = 5, T_max = 1.,  Δ = -0.5)
    HAVAR = H_A_BW(set)
    g = [1.,2.,3.,4.,5.]
    @btime QuantVarEntHam.cost_for_grad($g, $set, $HAVAR)
    println(QuantVarEntHam.cost_for_grad(g, set, HAVAR))
end

whole()


function test_owntype()
    set = Settings_XXZ(N=10, N_A = 4, T_max = 1.,  Δ = -0.5)
    g = [1.,2.,3.,4.]
    blks = [hi_XXZ(i, set) for i in 1:4]
    test = cust(blks)
    @btime cost_cust($g, $set, $test)
    println(cost_cust(g, set, test))
end

function seetype()
    set = Settings_XXZ(N=10, N_A = 5, T_max = 1.,  Δ = -0.5)
    H_Var(g,set) = H_A_BW(g, set)
    g = [1.,2.,3.,4.,5.]
    H_A = H_Var(g,set)
    @descend integrand(set, 1.,  H_A)
end

function test_cust_block()
    set = Settings_XXZ(N=10, N_A = 4, T_max = 1.,  Δ = -0.5)
    g = [1.,2.,3.,4.,5.]
    blks = [hi_XXZ(i, set) for i in 1:5]
    H_A = H_A_Var{Float64, 2}(5, rand(5), blks)
    @btime cost_cust($g, $set, $H_A)
    println(cost_cust(g, set, H_A))
end 

#test_cust_block()

#test()




#=
17.894 ms (21648 allocations: 27.80 MiB)
0.0007635403506645518

16.908 ms (21241 allocations: 27.79 MiB)
0.000763540350664588 Improvement: changed Real to Float64 for T_max and added specific density matrix constructor for rhoA

14.044 ms (4537 allocations: 20.61 MiB)
0.0007635403506645726 #Improvement: saved the matrices for the observables in settings (thus type stable integrand as a byprodudct :-) )

13.777 ms (4537 allocations: 20.61 MiB)
0.0007635403506646036 Improvement: added @inbounds for the array element accesses


=#


#=
18.642 ms (191168 allocations: 40.60 MiB)
0.0019478773231168265

16.828 ms (145810 allocations: 37.18 MiB)
0.001947877323116802 WITH COPY RHOA, NOT GOOD FOR MULTITHREADING 

16.608 ms (146110 allocations: 38.75 MiB)
0.0019478773231167998 # Improvement 1: removing creating HA in integrand but in cost function once

16.540 ms (145731 allocations: 38.74 MiB)
0.0019478773231168323 #Improvement 2: wrapping the blocks into a struct, removed it to test to create my own block 

17.201 ms (148031 allocations: 39.11 MiB)
0.0019478773231168104 #"Improvement" 2.1: define own block PROJECT FAILED

16.364 ms (145631 allocations: 38.74 MiB)
0.0019478773231168416 Improvement: using own cust struct and preallocate the time evolution block

16.364 ms (145631 allocations: 38.74 MiB)
0.0019478773231168416 Improvement: using own cust struct and preallocate the time evolution block

16.579 ms (145332 allocations: 37.17 MiB)
    0.001947877323116814 Improvement: using own cust struct and preallocate the time evolution block and give a copy of the density matrix and give rhoA as a copy

=#