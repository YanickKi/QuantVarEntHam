using SparseArrays

function X(S::Union{Rational, Int})
    d = Int64(2*S+1)

    multi = S*(S+1)

    Sx = spzeros(Float64, d, d)
    
    for m in -S:S
        _m = Int64(S-m+1)
        if m  < S
            m´ = m+1
            _m´ = Int64(S-m´+1)
            Sx[_m,_m´] = sqrt(multi - m´*(m´-1))
        end 
        if m  > -S
            m´ = m-1
            _m´ = Int64(S-m´+1)
            Sx[_m,_m´] = sqrt(multi - m´*(m´+1))
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
            Sy[_m,_m´] = - sqrt(multi - m´*(m´-1))
        end 
        if m  > -S
            m´ = m-1
            _m´ = Int64(S-m´+1)
            Sy[_m,_m´] = sqrt(multi - m´*(m´+1))
        end
    end 

    Sy ./= 1im

    return Sy

end 

function Z(S::Union{Rational, Int})
    d = Int64(2*S+1)

    Sz = Diagonal(zeros(d))

    for m in -S:S
        _m = Int64(m+S+1) 
        Sz[_m,_m] = - Float64(2*m)
    end 

    return Sz
end


function identity(S::Union{Rational, Int})
    d = Int64(2*S+1)

    I = Diagonal(zeros(Int64,d))

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