include("src/QuantVarEntHam.jl")

using .QuantVarEntHam
using LaTeXStrings
using CairoMakie

function whole()
    init = initialize(TFIM(8, 4, 1., 10., r_max = 1,rtol = 1e-7), H_A_BW)
    g = rand(0:1e-16:4., 4)
    QuantVarEntHam.comm_opt(g, init)
    #g, C = optimize_LBFGS(g, init, gtol = 1e-16)
    #print_H_A(g, init)
    #κ_var, κ_exact =  universal_ratios(g, init)
    #labels = ["var", "exact"]
    #plot_universal_ratios("F.pdf", labels, κ_var, κ_exact )
    #=init = initialize(XXZ(10, 5, 1., 2., 1, true), H_A_not_BW)
    ginit = [1., 1., 1.5, 1.5, 1.5, 1.5,1., 1.]
    g, C = optimize_LBFGS(ginit, init, gtol = 1e-16)
    #print_H_A(g, init)
    κ_var, κ_exact =  universal_ratios(g, init)
    init = initialize(XXZ(10, 5, 1., 2., 4, true), H_A_not_BW)
    ginit = [1., 1., 1.5, 1.5, 1.5, 1.5,1., 1., 0., 0., 0., 0., 0., 0., 0., 0.,0.,0.,0.,0.]
    g, C = optimize_LBFGS(ginit, init, gtol = 1e-16)
    #print_H_A(g, init)
    κ_var_corr, κ_exact =  universal_ratios(g, init)
    F = Figure()
    ax = Axis(F[1,1])
    scatter!(ax, κ_var, label = "var")
    scatter!(ax, κ_var_corr, label = "corr")
    scatter!(ax, κ_exact, label = "exact")
    axislegend(ax)
    save("F.pdf", F) =#
end 

whole()