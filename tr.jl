using LinearAlgebra
using BenchmarkTools

function main()
    n = 100
    N = 40
    As = [rand(n,n) for _ in 1:N]
    g = [rand() for _ in 1:N] 
    As_ntuple = tuple(As...)
    C = zero(As[1])
    @btime lincomb($C, $As, $g)
    @btime lincomb($C, $As_ntuple, $g)
end 


function lincomb(C, As, g::Vector{<:AbstractFloat})
    C .= g[1] .* As[1]

    @fastmath @inbounds @simd for i in 2:length(g)
        C .+= g[i].* As[i]
    end 
    return C
end 
main()

