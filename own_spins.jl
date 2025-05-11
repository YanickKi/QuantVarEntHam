using LinearAlgebra
using SparseArrays

function sigmax()
    d = 3

    Sx = spzeros(ComplexF64, d, d)
        
    Sx[1,2] = 1. 
    Sx[2,1] = 1.
    Sx[2,3] = 1.
    Sx[3,2] = 1.

    Sx ./=  sqrt(2.)
    return Sx
end 

function sigmay()
    d = 3

    Sy = spzeros(ComplexF64, d, d)
        
    Sy[1,2] = 1. 
    Sy[2,3] = 1.
    Sy[2,1] = -1.
    Sy[3,2] = -1.

    Sy ./= sqrt(2.) * 1im
    return Sy

end 

function sigmaz()
    d = 3

    Sz = spzeros(ComplexF64, d, d)
    Sz[1,1] = 1.0
    Sz[3,3] = -1.0

    return Sz
end


function identity()
    d = 3

    I = spzeros(ComplexF64, d, d)

    I[1,1] = 1.
    I[2,2] = 1.
    I[3,3] = 1.

    return I

end 


function repeat_cust(N, sig::Function, loc::AbstractVector{<:Int})
    
    M = spzeros(ComplexF64, 3, 3)

    if 1 ∉ loc 
        M = identity()
    else 
        M = sig()
    end 

    for qudit in 2:N
        if qudit ∉ loc 
            M = kron(M, identity())
        end
        if qudit ∈ loc 
            M = kron(M, sig())
        end
    end
    return M
end 

function main()

    M = repeat_cust(2, sigmaz, [1])
    println(M)
end 

main()