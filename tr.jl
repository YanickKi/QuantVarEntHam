using BenchmarkTools
using LinearAlgebra

function trace_two_matrices(A::AbstractMatrix, B::AbstractMatrix)
    n = size(A)[1]

    trace = zero(promote_type(eltype(A),eltype(B)))

    @inbounds for i in 1:n, j in 1:n 
        trace += A[i,j]*B[j,i]
    end 

    return trace
end 

function trace_two_matrices_good(A::AbstractMatrix, B::AbstractMatrix)
    n = size(A)[1]

    B_ = transpose(B)

    trace = zero(promote_type(eltype(A),eltype(B_)))

    @inbounds for j in 1:n, i in 1:n 
        trace += A[i,j]*B_[i,j]
    end 

    return trace
end 

function expect(Op::AbstractMatrix, ρ::AbstractMatrix)
    n = size(Op)[1]

    E = zero(promote_type(eltype(Op),eltype(ρ)))

    @inbounds for j in 1:n, i in j+1:n 
        E += Op[i,j]*conj(ρ[i,j])
    end 

    E *= 2 

    @inbounds for i in 1:n 
        E += Op[i,i]*ρ[i,i]
    end 


    #if abs(imag(E)) > eps(Float64)
    #    error("Imaginary part too large for expectation value!")
    #end 

    return real(E)
end 

function main()
    n = 2^8
    A = rand(ComplexF64, n, n)
    B = rand(ComplexF64, n, n)
    C = rand(ComplexF64, n, n)

    D = zero(A)

    A_ = A + A'
    B_ = B + B'
    C_ = C + C'    
    #@btime trace_two_matrices($A_,$B_)
    #@btime trace_two_matrices_good($A_,$B_)
    A__ = Hermitian(A_)
    B__ = Hermitian(B_)
    @btime dot($A_,$B_,$C_)
    @btime dot($A_, mul!($D, $B_, $C_))
    println(dot(A_,B_,C_))
    println(dot(A_, mul!(D, B_, C_)))
    #@btime expect($A_, $B_)
    #println(trace_two_matrices(A_,B_))
    #println(trace_two_matrices_good(A_,B_))
    #println(tr(A_*B_))
    #println(dot(A_,B_))

end 

main()