using QuantVarEntHam
using DelimitedFiles

function opt(N, ginit, dt)
    lines = readlines("vectors/N=$(N)_J=1_Bx=0.0_Uzz=0.0_periodic.txt")
    gs_v = ComplexF64[eval(Meta.parse(line)) for line in lines]
    A = [i for i in (N ÷ 2)+1:N ]
    ρ_A = QuantVarEntHam.partial_trace_pure_state(gs_v, 3^(N ÷2), 3^(N ÷ 2))
    init = initialize(pollmann(N, N ÷ 2, 1, 0, 0, 2, periodic = false, ρ_A = ρ_A), H_A_not_BW)
    g,c = optimize_LBFGS(ginit, init, methods = midpoint(dt = dt))
    writedlm("data/g_N=$(N)_J=1_Bx=0_Uzz=0_dt=dt.txt", g)
    writedlm("data/C_N=$(N)_J=1_Bx=0_Uzz=0_dt=dt.txt", c)

end 

function main()
    N = 18
    g = [1,1,1]
    dt = 0.1
    opt(N, g, dt)
end 

main()