using BenchmarkTools
using LinearAlgebra

function ownexpect_memory_corr(Op::AbstractMatrix, ρ::AbstractMatrix)
    n = size(Op)[1]

    E = 0.

    for j in 1:n 
        for i in 1:n 
            E += Op[i,j]*conj(ρ[i,j])
        end 
    end 

    return real(E)
end 

function fastest(Op, ρ)
    n = size(Op)[1]

    E = 0.

    @inbounds for j in 1:n, i in 1:n 
        E += Op[i,j]*conj(ρ[i,j])
    end 

    return real(E)

end 

function fastest_upper_wo(Op, ρ)
    n = size(Op)[1]

    E = 0.

    @inbounds for j in 1:n, i in j+1:n 
        E += Op[i,j]*ρ[j,i]
    end 

    E *= 2 

    @inbounds for i in 1:n 
        E += Op[i,i]*ρ[i,i]
    end 

    return real(E)

end 

function fastest_upper(Op, ρ)
    n = size(Op)[1]

    E = 0.

    @inbounds for j in 1:n, i in j+1:n 
        E += Op[i,j]*conj(ρ[i,j])
    end 

    E *= 2 

    @inbounds for i in 1:n 
        E += Op[i,i]*ρ[i,i]
    end 

    return real(E)

end 

function fastest_upper_coloumnmajor(Op, ρ)
    n = size(Op)[1]

    E_off = 0.
    E_diag = 0.
    @inbounds for j in 1:n
        @inbounds for i in j+1:n 
            E_off += Op[i,j]*conj(ρ[i,j])
        end 
        E_diag += Op[j,j]*conj(ρ[j,j])
    end 

    E = 2 * E_off + E_diag

    return real(E)

end 

function ownexpect(Op::AbstractMatrix, ρ::AbstractMatrix)
    n = size(Op)[1]


    E = 0.

    for j in 1:n 
        for i in 1:n 
            E += Op[i,j]*ρ[j,i]
        end 
    end 
 
    return real(E)
end



function main()
    n = 1000
    A = rand(ComplexF64, n, n)
    B = rand(ComplexF64, n, n)
    A_ = A + A'
    B_ = B + B'
     @btime fastest_upper($A_,$B_)
     @btime fastest_upper_wo($A_,$B_)
        println(fastest_upper(A_,B_))
            println(fastest_upper_wo(A_,B_))
    println(tr(A_*B_))
end 

main()
