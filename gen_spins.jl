using SparseArrays
using LinearAlgebra
function X(S::Union{Rational, Int})
    d = Int64(2*S+1)

    multi = S*(S+1)

    Sx = spzeros(ComplexF64, d, d)
    
    for m in -S:S
        _m = Int64(S-m+1)
        if m  < S
            m´ = m+1
            _m´ = Int64(S-m´+1)
            Sx[_m,_m´] = 1/2 * sqrt(multi - m´*(m´-1))
        end 
        if m  > -S
            m´ = m-1
            _m´ = Int64(S-m´+1)
            Sx[_m,_m´] = 1/2 * sqrt(multi - m´*(m´+1))
        end
    end 

    return Sx
end 

function Y(S::Union{Rational, Int})
    d = Int64(2*S+1)

    multi = S*(S+1)

    Sy = spzeros(ComplexF64, d, d)
    
    for m in -S:S
        _m = Int64(S-m+1)
        if m  < S
            m´ = m+1
            _m´ = Int64(S-m´+1)
            Sy[_m,_m´] = - 1/2 * sqrt(multi - m´*(m´-1))
        end 
        if m  > -S
            m´ = m-1
            _m´ = Int64(S-m´+1)
            Sy[_m,_m´] = 1/2 * sqrt(multi - m´*(m´+1))
        end
    end 

    Sy ./= 1im

    return Sy

end 

function Z(S::Union{Rational, Int})
    d = Int64(2*S+1)

    Sz = spzeros(ComplexF64, d, d)

    for m in -S:S
        _m = Int64(m+S+1) 
        Sz[_m,_m] = -m
    end 

    return Sz
end


function identity(S::Union{Rational, Int})
    d = Int64(2*S+1)

    I = spzeros(ComplexF64, d, d)

    for i in 1:d
        I[i,i] = 1.
    end 

    return I

end 


repeat(N, sig::Function, locs::Int; S::Union{Rational, Int} = 1//2) = repeat(N, sig, (locs,), S = S) 

function repeat(N, sig::Function, locs::NTuple{l,Int}; S::Union{Rational, Int} = 1//2) where {l}
    if maximum(locs) > N || minimum(locs) < 1
        error("Location out of bounds!")
    end 
    
    M = spzeros(ComplexF64, 3, 3)

    if 1 ∉ locs 
        M = identity(S)
    else 
        M = sig(S)
    end 

    for qudit in 2:N
        if qudit ∉ locs 
            M = kron(M, identity(S))
        end
        if qudit ∈ locs 
            M = kron(M, sig(S))
        end
    end
    return M
end  


function main()
    M = repeat(1, Z, (1,), S=1//2)
    values, vectors = eigen(Matrix(M))
    println(values)
end 

main()