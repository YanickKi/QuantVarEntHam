using LinearAlgebra
using BenchmarkTools


function nobuffer(A, B)
    return sum(A .* B)
end 

function buffer(C, A, B)
    return sum(had!(C, A, B))
end 

function had!(buff::AbstractMatrix, A::AbstractMatrix,B::AbstractMatrix)
    n = size(A)[1]
    @fastmath @inbounds @simd for j in 1:n
       for i in 1:n
         @inbounds buff[i,j] = A[i,j] *B[i,j]
       end
    end
    return buff
end

function main()
    n = 1000
    A = rand(ComplexF64, n, n)
    B  = rand(ComplexF64, n, n)
    C = zero(A)
    @btime nobuffer($A,$B)
    @btime buffer($C, $A, $B)
end 

main()