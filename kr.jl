using KrylovKit: eigsolve
using LuxurySparse

d = 2^4
Atemp = rand(ComplexF64, d,d)

A = (Atemp + Atemp')/2

for i in 1:3
    values, vectors = eigsolve(A ,1 ,:SR, ishermitian=true)
    #println(vectors[1,1])
    phase = angle(vectors[1][1])
    v_normalized = vectors[1] * exp(-im * phase)
    println(v_normalized)
end 
0.03514129052816156 - 0.0020962782953862376im
-0.03434793101362745 + 0.007715200435412491im
0.032346155206583495 - 0.01389355699696334im