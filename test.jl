using Distributed
using Yao 
using Parameters
#addprocs(5)

include("src/QuantVarEntHam.jl")
using .QuantVarEntHam
using BenchmarkTools

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

#test_owntype()


function test_str()
    set = Settings_XXZ(N = 12, N_A = 6, T_max = 5, Δ = -0.5, r_max = 1, periodic = true, signHam = +1)
    HA_VAR = H_A_not_BW(set) 
    g_init = [0.5, 1., 1., 1.5, 1.5, 1.,1.,0.5,0.5] 
    G = zeros(length(g_init))
    #QuantVarEntHam.fg_fixed!(1.,G, g_init, set, HA_VAR, 1.)
    optimize_LBFGS(g_init, set, HA_VAR, g1 = 0.5, gtol = 1e-16)
    #@btime QuantVarEntHam.cost($g_init, $set, $HA_VAR, false)
    #@btime QuantVarEntHam.cost_with_QuadGK($g_init, $set, $HA_VAR,  false)
    #println(QuantVarEntHam.cost(g_init, set, HA_VAR, false))
    #println(QuantVarEntHam.cost_with_QuadGK(g_init, set, HA_VAR,  false))
end 

test_str()

function testexmpv()
    set = Settings_XXZ(N=10, N_A = 6, T_max = 10.,  Δ = -0.5)
    g = [1.,2.,3.,4., 5., 6.]
    blks = [hi_XXZ(i, set) for i in 1:6]
    H = Matrix(sum(g.*blks))
    @btime exponential!(-1im*1.0*$H)
    @btime exp(-1im*1.0*$H)

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


