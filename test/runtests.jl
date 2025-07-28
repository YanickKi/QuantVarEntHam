using QuantVarEntHam
using Test

@testset "Phase transition at Γ=1 for TFIM" begin 
    model = TFIM(8,4,1)
    ansatz = H_A_BWV(model)
    cost = QCFL(model, ansatz, 1, integrator = MidPoint(1e-2))
    g_init = [1,1,2,2,2,3,3]
    g_opt, c_opt = optimize(cost, g_init, print_result=false, show_trace = false)
    @test c_opt ≈ 0. atol = 1e-16
end 

@testset "Classical state at Γ=0 for TFIM" begin
    model = TFIM(8,4, 0)
    ansatz = H_A_BW(model)
    cost = QCFL(model, ansatz, 1)
    for _ in 1:10
        g = rand(0.:1e-16:5., 4)
        @test cost(g) ≈ 0. atol = eps(Float64)
    end
end 

@testset "Classical state at Δ=-2 for XXZ model" begin
    model = XXZ(8,4,-2)
    ansatz = H_A_BW(model)
    cost = QCFL(model, ansatz, 1)
    for _ in 1:10
        g = rand(0.:1e-16:5., 4)
        @test cost(g) ≈ 0. atol = eps(Float64)
    end
end 

@testset "Check cost function results equality at Γ=1 for TFIM" begin
    model = TFIM(8,4,1)
    ansatz = H_A_BWV(model)
    cost1 = QCFL(model, ansatz, 1, integrator = MidPoint(1e-2))
    cost2 = Commutator(model, ansatz)
    cost3 = RelativeEntropy(model, ansatz)
    costs = [cost1, cost2, cost3]
    g_init = [1,1,2,2,2,3,3]
    g_opt_all = Vector{Float64}[]
    for cost in costs 
        g_opt, c_opt = optimize(cost, g_init, print_result=false, show_trace = false)
        push!(g_opt_all, g_opt/g_opt[1])
    end     

    @test g_opt_all[1] ≈ g_opt_all[2] atol = 1e-2
    @test g_opt_all[1] ≈ g_opt_all[3] atol = 1e-2
    @test g_opt_all[2] ≈ g_opt_all[3] atol = 1e-2
end 

@testset "Fixing parameter at Γ = 1 for TFIM" begin 
    model = TFIM(8,4,1)
    ansatz = H_A_BWV(model)
    free_cost = QCFL(model, ansatz, 1, integrator = MidPoint(1e-2))
    fixed_cost = FixedCost(free_cost, 1, 1)
    g_init_free = [1,1,2,2,2,3,3]
    g_init_fixed = [1,2,2,2,3,3]
    g_opt_free, c_opt_free = optimize(free_cost, g_init_free, print_result=false, show_trace = false)
    _g_opt_fixed, c_opt_fixed = optimize(fixed_cost, g_init_fixed, print_result=false, show_trace = false)

    g_opt_fixed = fill_full_g(fixed_cost, _g_opt_fixed)

    @test c_opt_free ≈ c_opt_fixed atol = eps(Float64)
    @test g_opt_free/g_opt_free[1] ≈ g_opt_fixed/g_opt_fixed[1] atol = sqrt(eps(Float64))
end 

@testset "Mid point = Tanh-sinh at Γ=1 for TFIM" begin
    model = TFIM(8,4,1)
    ansatz = H_A_BWV(model)
    cost_mp = QCFL(model, ansatz, 1, integrator = MidPoint(1e-2))
    cost_ts = QCFL(model, ansatz,1)
    g_init = [1,1,2,2,2,3,3]
    g_mp, c_mp = optimize(cost_mp, g_init, print_result=false, show_trace = false)
    g_ts, c_ts = optimize(cost_ts, g_init, print_result=false, show_trace = false)

    @test c_mp ≈ c_ts atol = eps(Float64)
    @test g_mp/g_mp[1] ≈ g_ts/g_ts[1] atol = sqrt(eps(Float64))

end 