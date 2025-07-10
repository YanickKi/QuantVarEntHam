export PauliString
export X, Y, Z
export mat
using SparseArrays
using LinearAlgebra

import Base: *
import Base: +
import Base: ^
import Base: convert

abstract type AbstractBlock{S,N} end

struct PauliString{S, N, L} <: AbstractBlock{S,N}
    sig::String
    locs::NTuple{L, Int}
    power::Int
end

struct Block{S,N} <: AbstractBlock{S,N}
    prefactors::Vector{Float64}
    pauli_strings::Vector{PauliString{S,N}}
end 


PauliString(N, sig::String, locs::Int; S::Union{Rational, Int} = 1//2) = PauliString(N, sig, (locs, ), S = S) 


function PauliString(N, sig::String, locs::NTuple{L,Int}; S::Union{Rational, Int} = 1//2) where {L}
    ps = PauliString{Rational(S),N, L}(sig, locs, 1)
    return ps
end 



function convert(::Type{Block{S,N}}, ps::PauliString{S,N, L}) where {S,N, L}
    Block{S,N}([1.],[ps])
end 

function convert(::Type{Block{S,N}}, block::Block{S,N}) where {S,N}
    return block
end 


include("algebra.jl")
include("matrices.jl")
include("printing.jl")