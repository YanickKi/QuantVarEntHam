using Distributed
using Yao 
using Parameters
using Cthulhu
include("src/QuantVarEntHam.jl")
using .QuantVarEntHam
using BenchmarkTools
using LinearAlgebra
function test_owntype()
    N_A =  3
    set = Settings_XXZ(N=2*N_A, N_A = N_A, T_max = 5.,  Δ = -0.5)
    g = [1.5, 2.]
    blks = [hi_XXZ(i, set) for i in 1:3]
    test = utilBlocks(blks, Matrix.(blks), time_evolve(put(N_A, 1=>Z) + put(N_A, 1=>Z), 1., check_hermicity = false))
    g,c = optimize_LBFGS(g, set, test, g1=1.0, gtol =  1e-16)
    println(g,c)
    #G = zeros(length(g))
    #println(grad_fixed(G, g, set, test, 1.))
    #println(g_fixed!(G, g, 1., set, eps(Float64)^(1/3), test))
    #
    #@btime grad($G, $g, $set, $test)
    #@btime g_fixed!($G, $g, 1., $set, eps(Float64)^(1/3), $test)


end

function testintegrate()
    set = Settings_XXZ(N=10, N_A = 5, T_max = 5.,  Δ = -0.5)
    HAVAR = H_A_BW(set)
    g = [1.,2.,3.,4.,5.]
    @unpack q = set 
    println(q)
end
#test_owntype()




function testintegrate()
    set = Settings_XXZ(N=10, N_A = 5, T_max = 2.,  Δ = -0.5)
    HAVAR = H_A_BW(set)
    g = [1.,2.,3.,4.,5.]
    te = QuantVarEntHam.tes(zeros(ComplexF64, 2^5, 2^5), zeros(ComplexF64, 2^5, 2^5) ,zeros(ComplexF64, 2^5, 2^5) ,zeros(ComplexF64, 2^5, 2^5) ,
    [zeros(ComplexF64, 2^5, 2^5) for i in 1:4] )
    println(QuantVarEntHam.cost_with_midpointrule(g, set, HAVAR, te))
    println(QuantVarEntHam.cost(g, set, HAVAR))
    H_A = @inbounds sum(g[i]*HAVAR.matrices[i] for i in eachindex(g))
    @btime copyto!($te.Uhalve, exp(-1im*0.001/2*$H_A))

    @btime QuantVarEntHam.cost_with_midpointrule($g, $set, $HAVAR, $te)
    @btime QuantVarEntHam.cost($g, $set, $HAVAR)


end


function dia()
    set = Settings_XXZ(N = 10, N_A = 5, Δ = -0.5, T_max = 1.)
end 

function preintegrate()
    set = Settings_XXZ(N=10, N_A = 5, T_max = 1.,  Δ = -0.5)
    HAVAR = H_A_BW(set)
    g = [1.,2.,3.,4.,5.]
    @btime QuantVarEntHam.cost_for_grad($g, $set, $HAVAR)
    println(QuantVarEntHam.cost_for_grad(g, set, HAVAR))
end

function test_str()
    set = Settings_XXZ(N = 12, N_A = 5, T_max = 5, Δ = -0.5, r_max = 1, periodic = true, signHam = +1)
    HA_VAR = H_A_not_BW(set) 
    g_init = [0.5, 0.5, 1., 1., 1.,1.,0.5,0.5] 
    G = zeros(length(g_init))
    @btime QuantVarEntHam.fg!(1.,$G, $g_init, $set, $HA_VAR)
    #@btime QuantVarEntHam.cost($g_init, $set, $HA_VAR, false)
    #@btime QuantVarEntHam.cost_with_QuadGK($g_init, $set, $HA_VAR,  false)
    #println(QuantVarEntHam.cost(g_init, set, HA_VAR, false))
    #println(QuantVarEntHam.cost_with_QuadGK(g_init, set, HA_VAR,  false))
end 


function typ()
    set = Settings_XXZ(N=4, N_A=2, Δ = 0.5, T_max = 1.)
    g = [1.,2.]
    HAVAR = H_A_BW(set)
    F = 0.
    G = [0.,0.]
    @descend QuantVarEntHam.cost(g, set, HAVAR)
end


function testexmpv()
    set = Settings_XXZ(N=10, N_A = 6, T_max = 10.,  Δ = -0.5)
    g = [1.,2.,3.,4., 5., 6.]
    blks = [hi_XXZ(i, set) for i in 1:6]
    H = Matrix(sum(g.*blks))
    @btime exponential!(-1im*1.0*$H)
    @btime exp(-1im*1.0*$H)
1., sum(g[i]*HAVAR.matrices[i] for i in eachindex(g))
end


function te_state(rhoA::DensityMatrix, H::AbstractBlock, dt::AbstractFloat)
    st = state(rhoA)
    A = Yao.YaoBlocks.BlockMap(ComplexF64, H)
    @inbounds for j = 1:size(st, 2)
        v = view(st, :, j)
        v_out = exponentiate(A, dt, v,
                                   tol=1e-7,
                                   krylovdim=min(1000, size(A,1)),
                                   ishermitian=true,
                                   eager=true,)[1]
        # info.converged is 1 if it converged and 0 otherwise
        v .= v_out
    end
    return st
end 

function te_state_cust(rhoA::DensityMatrix, H::AbstractBlock, dt::AbstractFloat)
    U = exp(-1im*dt*Matrix(H))
    return U*rhoA.state*U'

end

function yaote(rhoA::DensityMatrix, H::AbstractBlock, dt::AbstractFloat)
    U = time_evolve(H, dt, check_hermicity = false)
    apply!(rhoA, U)
end 

function testBlock()
    set = Settings_XXZ(N=10, N_A = 6, T_max = 1.,  Δ = -0.5)
    g = [1.,2.,3.,4.,5.,6.]
    blks = [hi_XXZ(i, set) for i in 1:6]
    H = sum(g.*blks)
    @unpack rhoA = set 
    dt = 1.
    @btime te_state(copy($rhoA), $H, $dt)
    @btime te_state_cust(copy($rhoA), $H, $dt)
    @btime yaote(copy($rhoA), $H, $dt)

end 

using ChainRules, ChainRulesCore

function mul_test()
    init = initialize(TFIM(10, 5, 1., 1., dt = 0.001, r_max = 1), H_A_BW)
    g = [1.,2.,3.,4.,5.]
    G = zeros(length(g))
    g, c = optimize_LBFGS(g, init, gtol = 1e-16)
    κ_var, κ_exact = universal_ratios(g, init)
    κ_var, κ_svd = QuantVarEntHam.universal_ratios_svd(g, init)
    plot_universal_ratios("F.pdf", ["κ_var", "κ_exact", "κ_SVD"], κ_var, κ_exact, κ_svd)
    #println(typeof(init.blks.matrices))
    #@btime QuantVarEntHam.cost_grad!(1., $g, $G, $init)

    #optimize_LBFGS(g, init, gtol = 1e-16, maxiter = 100, print_result = true, show_trace = true)
    #println("QuadGK:" , QuantVarEntHam.cost_grad_quadgk!(1., g, nothing, init))
    #println("tanh-sinh: ", QuantVarEntHam.cost_grad!(1., g, nothing, init))
    #println("hcubature: ",  QuantVarEntHam.cost_grad_hcubature!(1., g, nothing, init))
    #println("midpoint: ", QuantVarEntHam.cost_grad_midpoint!(1., g, nothing, init))
    #println("QuadGK: ",  @btime QuantVarEntHam.cost_grad_quadgk!(1., $g, nothing, $init))
    #println("tanh-sinh: ",  @btime QuantVarEntHam.cost_grad!(1., $g, nothing, $init))
    #println("hcubature: ",  @btime QuantVarEntHam.cost_grad_hcubature!(1., $g, nothing, $init))
    #println("midpoint: ", @btime QuantVarEntHam.cost_grad_midpoint!(1., $g, nothing, $init))

    #QuantVarEntHam.optimize_quadgk(g, init)
    #G = rand(length(g))
    #println(typeof(init.set.mtrxObs))
    #QuantVarEntHam.cost_grad_midpoint!(1., $g, nothing, $init)
    #c = QuantVarEntHam.cost_grad!(1., g, nothing, init)
    #println(c)    
    #@descend QuantVarEntHam.cost_grad!(1., g, G, init)
    #@code_warntype QuantVarEntHam.integrand(1., init)
end 
#34.137 ms (13846 allocations: 64.84 MiB)
#= init = initialize(TFIM(14, 7, 1., 1.), H_A_BW), ginit = [1., 2., 3., 4., 5., 6., 7.]
vor der nutzung von der struktur der matrizen: 746.306 ms (14631 allocations: 1.05 GiB)
668.967 ms (20287 allocations: 1.05 GiB) matrobs mit mat aber ohne TS
661.708 ms (22105 allocations: 1.05 GiB) matrobs mit mat und buffer diagonal ohne TS
680.973 ms (14631 allocations: 1.05 GiB) matrobs mit mat und buffer diagonal mit TS
663.747 ms (16163 allocations: 1.05 GiB) auch für die matrices in blocks die mat funktion benutzt ohne TS
647.123 ms (14695 allocations: 1.05 GiB) done everything and everyting is type stable
679.109 ms (14695 allocations: 1.05 GiB) everyting type stable and parametric buffers, why slower???????
683.174 ms (14695 allocations: 1.05 GiB) before buffer for -1im*t*H_A
680.699 ms (14392 allocations: 1.03 GiB) after buffer for -1im*t*H_A
671.513 ms (14190 allocations: 1.03 GiB) after summing to c without sum and broadcasting 
650.395 ms (14129 allocations: 1.03 GiB) after making blocks dense again 
649.219 ms (14126 allocations: 1.03 GiB) no Hermitian Matrix constructor in cost_grad anymore 
665.716 ms (14108 allocations: 1.03 GiB) loop to make H_A not the sum Function
660.618 ms (14087 allocations: 1.03 GiB) working in place for skalar times matrix
=#
#Everything with T_max = 2. here 
#before optimizing 40.244 ms (17738 allocations: 66.21 MiB) (type instability)
#fixing bug with static arrays but type instable: 39.218 ms (14308 allocations: 66.14 MiB)
#using another buffer 39.284 ms (14294 allocations: 66.14 MiB)
# for T_max = 1 : 18.798 ms (7047 allocations: 32.49 MiB)
mul_test()