using DelimitedFiles
using QuantVarEntHam
using LinearAlgebra

function main(N)
    lines = readlines("vectors/N=$(N)_J=1_Bx=0.0_Uzz=0.0_periodic.txt")
    gs_v = ComplexF64[eval(Meta.parse(line)) for line in lines]
    ρ_A = QuantVarEntHam.partial_trace_pure_state(gs_v, 3^(N ÷ 2), 3^(N ÷ 2), trace_subsystem = :A)
    p, v = eigen(Hermitian(ρ_A))
    if p[1] < 0. 
        p .+= abs(p[1]) + eps(Float64)
    end 
    ξ_exact= -log.(p)
    ξ_exact = reverse(ξ_exact)
    ur =  QuantVarEntHam.calc_universal_ratios(ξ_exact, 1, 5)
    writedlm("N=$(N)_J=1_Bx=0.0_Uzz=0.0_periodic.txt", ur)
end 

main(4)