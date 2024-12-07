using Yao 
using BenchmarkTools
using LinearAlgebra
using SparseArrays
function whole()
    d = 10
    N = 2^d
    #Ap = sum(rand()*mat(repeat(d, Z,(i,i+1))) for i in 1:d-1)
    Ad = sum(rand()*mat(repeat(d, X,(i,i+1))) for i in 1:d-1)
    #Ap = [mat(repeat(d, Z,(i,i+1))) for i in 1:d-1]
    #Ad = [mat(repeat(d, X,(i,i+1))) for i in 1:d-1]
    println(typeof(Ad))

    #Acombined = Ap
    #Ae = sum(rand()*mat(put(d, i=>X)) for i in 1:d)
    Ac = Ap+Ad
    #Bp = mat(put(10, 1=>Z))
    #A = Matrix(repeat(10, X,(1,2)))
    #B = Matrix(put(10, 1=>Z))
    #C = zeros(ComplexF64, 2^10, 2^10)
    #D = rand(ComplexF64, 2^10, 2^10)
    #@btime mul!($C, $A, $B)
    #@btime $Ap * $Bp
    #@btime $D * $Ap
    #println(typeof(Ad))
    A = rand(ComplexF64, N, N)
    B = sprand(ComplexF64, N, N, 0.0078125)
    C = rand(ComplexF64, N, N)
    D = zeros(ComplexF64, N, N)
    #println(typeof(Ap*A))
    #println(typeof(Ap))
    #fill!(Ap, 0)
    #println(typeof(Ap))
    println(typeof(Acombined))
    Asum = sum(rand()*Acombined[i] for i in eachindex(Acombined))
    println(typeof(Asum))
    #@btime mul!($D, $A,$Ap)
    @btime $A*$Asum 
    @btime mul!($D, $Asum, $A)
    #@btime mul!($D, $Ac, $A)
    #@btime mul!($D, $A, $C)
    ##@btime mul!($D,$A,$C)
end 
whole()