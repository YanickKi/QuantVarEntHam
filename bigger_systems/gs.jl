using DelimitedFiles
using QuantVarEntHam 
using KrylovKit: eigsolve

function get_gs(H, N, J, Bx, Uzz) 
    values, vectors = eigsolve(H ,1 ,:SR, ishermitian=true, tol = 1e-16)
    writedlm("vectors/N=$(N)_J=$(J)_Bx=$(Bx)_Uzz=$(Uzz)_periodic.txt", vectors[1])
end

function whole(N)
    J = 1
    Bx = 0.
    Uzz = 0.
    H =  H_pollmann(N, J, Bx, Uzz, periodic = true)
    get_gs(H,N, J, Bx, Uzz)
end 

whole(4)