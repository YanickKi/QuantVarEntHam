using QuantVarEntHam
using QuantumOpticsBase
using DelimitedFiles

function opt(N, ginit, dt)
    lines = readlines("vectors/N=$(N)_J=1_Bx=0.0_Uzz=0.0_open.txt")
    gs_v = ComplexF64[eval(Meta.parse(line)) for line in lines]
    A = [i for i in (N ÷ 2)+1:N ]
    b = SpinBasis(1)
    OPH = [b for i in 1:N]
    B = tensor(OPH...)
    gs = Ket(B, gs_v)
    ρ_A = ptrace(gs ⊗ dagger(gs),setdiff(1:N,A))
    init = initialize(pollmann(N, N ÷ 2, 1, 0, 0, 2, periodic = false, ρ_A = ρ_A), H_A_not_BW)
    g,c = optimize_LBFGS(ginit, init, methods = midpoint(dt = dt))
    writedlm("data/g_N=$(N)_J=1_Bx=0_Uzz=0_dt=dt.txt", g)
    writedlm("data/C_N=$(N)_J=1_Bx=0_Uzz=0_dt=dt.txt", c)

end 

function main()
    N = 8 
    g = [1,1,1,2,2,2,3,3,3]
    dt = 0.1
    opt(N, g, dt)
end 

main()