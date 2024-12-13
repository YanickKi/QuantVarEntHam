ENV["JULIA_CUDA_CPU_MOCK"] = true
using CUDA
using LinearAlgebra


function whole()
    N = 5
    d = 2^N 
    A = cu(rand(ComplexF64, d,d))
    B = cu(rand(ComplexF64, d,d))
    C = cu(zeros(ComplexF64, d,d))
    println(typeof(A))
    mul!(C, A, B)
end

whole()