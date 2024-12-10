include("src/QuantVarEntHam.jl")
using BenchmarkTools
using .QuantVarEntHam
using LaTeXStrings
using CairoMakie
using LinearAlgebra
function whole()
    init = initialize(XXZ(10, 5, 0.5, 1., periodic = false, signHam = +1), H_A_BW)
    g = [3., 2., 3., 4., 5.,]
    #@btime QuantVarEntHam.cost_grad!(1., $g, nothing, $init)
    QuantVarEntHam.optimize_T_max_g1_fixed_custrgad(g, init, 1., gtol = 5e-12)
    #println(init.set.ρ_A.state)
    #println(QuantVarEntHam.cost_T_max(g, init))
    #init = initialize(TFIM(13, 4, 1., 3., periodic = true, signHam = +1), H_A_not_BW)
    #g1 = 0.46945880128558726
    #g = vcat(5., [0.4185400764298917,	0.6120789944202124,	0.6037784972153571,	0.4694588012855859,	0.41854007642989177])
    #println(QuantVarEntHam.cost_T_max(g, init))
    #g, c = QuantVarEntHam.optimize_T_max_g1_fixed(g, init, gtol = 1e-12)
    #g = vcat(3.000101936542932, 0.46945880128558726, [0.4328373623866064, 0.6381525828493494, 0.6008223362929461, 0.4755267947353266, 0.3945138071779576])
    #G = zeros(length(g)-1)
    #QuantVarEntHam.cost_grad_Tmax_g1_fixed!(1.,g, G, init )
    #println(G)
    #g, c = QuantVarEntHam.comm_opt_fixed(g, init, 1., 1e-12, 100)
    #g, C = optimize_LBFGS(g, init, gtol = 1e-16)
    #print_H_A(g, init)
    #κ_var, κ_exact =  universal_ratios(g, init)
    #labels = ["var", "exact"]
    #plot_universal_ratios("F.pdf", labels, κ_var, κ_exact )
    #init = initialize(XXZ(10, 5, 0.5, 2., r_max = 2,periodic = true), H_A_not_BW)
    #ginit = vcat([1., 1., 1.5, 1.5, 1.5, 1.5,1., 1.], zeros(6))
    #g, C = optimize_LBFGS(ginit, init, gtol = 1e-25)
    #print_H_A(g, init)
    #κ_var, κ_exact =  universal_ratios(g, init)
    #init = initialize(XXZ(10, 5, 1., 2., 4, true), H_A_not_BW)
    #ginit = [1., 1., 1.5, 1.5, 1.5, 1.5,1., 1., 0., 0., 0., 0., 0., 0., 0., 0.,0.,0.,0.,0.]
    #g, C = optimize_LBFGS(ginit, init, gtol = 1e-16)
    ##print_H_A(g, init)
    #κ_var_corr, κ_exact =  universal_ratios(g, init)
    #F = Figure()
    #ax = Axis(F[1,1])
    #scatter!(ax, κ_var, label = "var")
    #scatter!(ax, κ_var_corr, label = "corr")
    #scatter!(ax, κ_exact, label = "exact")
    #axislegend(ax)
    #save("F.pdf", F) 
end 

whole()