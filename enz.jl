using Enzyme
using BenchmarkTools
using LinearAlgebra

const x0 = 1.

function expv_taylor(t, A, B, degree_max)
    F = Z = B
    #norm_tail_old = _opnormInf(Z)
    for j in 1:degree_max
        Z = (A * Z) * (t / j)  # (t A)ʲ/j! * B
        F += Z
        # check if ratio of norm of tail and norm of series is below tolerance
        #norm_tail = _opnormInf(Z)
        #norm_tail_tot = norm_tail_old + norm_tail
        #norm_tail_tot ≤ tol * _opnormInf(F) && break
        #norm_tail_old = norm_tail
    end
    return F
end

function f(g)
    ψ = [1, 0]
    A = [g 0;0 -1]
    B = [1 0; 0 0]
    U = expv_taylor(-1im, A, B, 20)
    σ = [1 0; 0 -1]
    #C = expv(1, A, B)
    #B = [1 0; 0 1]
    #ϕ  = expv(-1im, A, ψ)
    #U = exp(-1im*A)
    #mul!(C, U, B)
    #println(ψ'*C*ψ)
    ϕ = U*ψ
    return ϕ'*σ*ϕ
end 

function whole()
    res = autodiff(Reverse, f, Active, Active(1.0))
    println(res)
    #@btime  ReverseDiff.gradient(f, [1.])
end

whole()
println(f(1.))
h = eps(Float64)^(1/3)
println((f(x0+h)-f(x0-h))/(2*h))