using LinearAlgebra
using BenchmarkTools


#function dia(F)
#    U = F.vectors * Diagonal(cis.(-F.values)) * F.vectors'
#end 

function dia(A)
    F = eigen(A)
    #U = F.vectors * Diagonal(cis.(-F.values)) * F.vectors'
end 

function main()
    n = 2^10
    A = rand(ComplexF64, n,n)
    B = Hermitian(A + A') 
    F = eigen(B)
    @btime dia($B)
    #@btime cis(-$B)
end 

main()