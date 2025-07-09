using SparseArrays
using LinearAlgebra

import Base: *

struct PauliString{S, N}
    sig::String
    locs::Vector{Int}
end

mutable struct Block{S,N} 
    prefactors::Vector{Float64}
    pauli_strings::Vector{PauliString{S,N}}
end 


# B = \sum_i a_i P_i 

mutable struct AddBlock{S, N}
    blocks::Vector{Block{S, N}}
end

function X(S::Union{Rational, Int})
    d = Int64(2*S+1)

    multi = S*(S+1)

    Sx = spzeros(ComplexF64, d, d)
    
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

    Sz = Diagonal{ComplexF64}(zeros(d))

    for m in -S:S
        _m = Int64(m+S+1) 
        Sz[_m,_m] = - Float64(2*m)
    end 

    return Sz
end



repeat(N, sig::Function, locs::Int; S::Union{Rational, Int} = 1//2) = repeat(N, sig, (locs,), S = S) 

function repeat(N, sig::Function, locs::NTuple{l,Int}; S::Union{Rational, Int} = 1//2) where {l}
    if maximum(locs) > N || minimum(locs) < 1
        error("Location out of bounds!")
    end 
    
    d = Int64(2*S+1)
    
    M = spzeros(ComplexF64, d, d)

    Identity = sparse((1.0+0.0im)*I,d,d)

    if 1 ∉ locs 
        M = Identity
    else 
        M = sig(S)
    end 

    for qudit in 2:N
        if qudit ∉ locs 
            M = kron(M, Identity)
        end
        if qudit ∈ locs 
            M = kron(M, sig(S))
        end
    end
    return sparse(M)
end 

PauliString(N, sig::Function, locs::Int; S::Union{Rational, Int} = 1//2) = PauliString(N, sig, (locs, ), S = S) 


function PauliString(N, sig::Function, locs::NTuple{l,Int}; S::Union{Rational, Int} = 1//2) where {l}
    
    sig_string = string(Symbol(sig))
    
    block = PauliString{Rational(S),N}(sig_string, [loc for loc in locs])
    return block
end 



function Base.show(io::IO, ps::PauliString{S, N}) where {S,N}
    println("Spin " , S)
    println("Number of spins: ", N)
    println()
    for loc in ps.locs 
        print(ps.sig, loc)
    end 
end


#=
function Base.show(io::IO, block::Block{S, N}) where {S,N}
    println("Spin " , S)
    println("Number of spins: ", N)
    println()
    block.prefactor == 1 ? nothing : print(string(block.prefactor), " *")
    for loc in block.locs 
        print(" ", block.sig, loc)
    end 
    #string_block = join.(string(block.prefactor), block.sig, block.locs)    
    #print(io, string_block)
end


function *(number::Real, pauli_string{S,N}::PauliString) where {S,N}
    block = Block{S,N}([number], [pauli_string])
    return block
end 

*(block::Block, number::Real) = *(number, block)

=#