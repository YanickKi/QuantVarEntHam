export PauliString
export X, Y, Z
export mat
using SparseArrays
using LinearAlgebra

import Base: *
import Base: +
import Base: ^


struct PauliString{S, N, L}
    sig::String
    locs::NTuple{L, Int}
    power::Int
end

struct Block{S,N} 
    prefactors::Vector{Float64}
    pauli_strings::Vector{PauliString{S,N}}
end 

struct BlocksMatrices{S,N, M<:AbstractMatrix}
    blocks::Vector{Block{S,N}}
    matrices::Vector{M}
end 

function BlocksMatrices(blocks::Vector{<:Block})
    return BlocksMatrices(blocks, mat.(blocks))
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


PauliString(N, sig::String, locs::Int; S::Union{Rational, Int} = 1//2) = PauliString(N, sig, (locs, ), S = S) 


function PauliString(N, sig::String, locs::NTuple{L,Int}; S::Union{Rational, Int} = 1//2) where {L}
    ps = PauliString{Rational(S),N, L}(sig, locs, 1)
    return ps
end 


function print_pauli_string(ps::PauliString)
    for _ in 1:ps.power
        for loc in ps.locs
            print(ps.sig, loc)
        end
    end 
end 

function Base.show(io::IO, ps::PauliString{S, N}) where {S,N}
    println("Spin " , S)
    println("Number of spins: ", N)
    println()
    print_pauli_string(ps)
end

function Base.show(io::IO, block::Block{S, N}) where {S,N}
    println("Spin " , S)
    println("Number of spins: ", N)
    println()
    for index in eachindex(block.pauli_strings) 
        block.prefactors[index] == 1 ? nothing : print(block.prefactors[index], " * ")
        print_pauli_string(block.pauli_strings[index])
        index == length(block.pauli_strings) ? nothing : print(" + ")
    end 
end

function *(number::Real, ps::PauliString{S,N}) where {S,N}
    block = Block{S,N}([Float64(number)], [ps])
    return block
end 

*(ps::PauliString{S,N}, number::Real) where {S,N} = *(number, ps)

function *(number::Real, block::Block{S,N}) where {S,N}
    new_block = Block{S,N}(Float64(number)*block.prefactors, block.pauli_strings)
    return new_block
end 

*(block::Block, number::Real) = *(number, block)

function +(ps1::PauliString{S,N}, ps2::PauliString{S,N}) where {S, N}
    block = Block{S,N}([1., 1.], [ps1, ps2])
    return block
end

function +(block::Block{S,N}, ps::PauliString{S,N}) where {S,N}
    new_block = Block{S,N}(vcat(block.prefactors, 1.), vcat(block.pauli_strings, ps))
    return new_block
end 

function +(ps::PauliString{S,N}, block::Block{S,N}) where {S,N}
    new_block = Block{S,N}(vcat(1., block.prefactors), vcat(ps,block.pauli_strings))
    return new_block
end 
function +(block1::Block{S,N}, block2::Block{S,N}) where {S,N}
    new_block =  Block{S,N}(vcat(block1.prefactors, block2.prefactors), vcat(block1.pauli_strings, block2.pauli_strings))
    return  new_block
end 

function ^(ps::PauliString{S,N,L}, pow::Int) where {S,N, L}
    return PauliString{S,N, L}(ps.sig, ps.locs, ps.power*pow)
end

function mat(block::Block{S,N}) where {S,N}
    d = Int((2*S+1)^N)

    M = spzeros(ComplexF64, d, d)

    for index in eachindex(block.pauli_strings)
        M .+= block.prefactors[index] .* (repeat(N, getfield(Main, Symbol(block.pauli_strings[index].sig)), block.pauli_strings[index].locs, S = S))^block.pauli_strings[index].power
    end 
    return M
end

function mat(ps::PauliString{S,N}) where {S,N}

    return repeat(N, getfield(Main, Symbol(ps.sig)), ps.locs, S = S)^(ps.power)

end