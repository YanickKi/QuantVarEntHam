using QuantVarEntHam
using Test

@testset "Classical state at Γ=0 for TFIM" begin 
    model = TFIM(8,4,1)
    ansatz = H_A_BWV(model)
    cost = QCFL(model, ansatz, 1, integrator = MidPoint(1e-2))
    g_init = [1,1,2,2,2,3,3]
    g_opt, c_opt = optimize(cost, g_init, print_result=false, show_trace = false)
    @test c_opt ≈ 0. atol = 1e-16
end 