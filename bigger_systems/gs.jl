using DelimitedFiles
using QuantVarEntHam 
using KrylovKit: eigsolve
using QuantumOpticsBase
function get_gs(H, N, J, Bx, Uzz) 
    values, vectors = eigsolve(H.data ,1 ,:SR, ishermitian=true, tol = 1e-16)
    println(N)
    writedlm("vectors/N=$(N)_J=$(J)_Bx=$(Bx)_Uzz=$(Uzz)_open.txt", vectors[1])
end

function whole(N)
    J = 1
    Bx = 0.
    Uzz = 0.
    H =  H_pollmann(N, J, Bx, Uzz, periodic = false)
    get_gs(H,N, J, Bx, Uzz)
end 

whole(8)