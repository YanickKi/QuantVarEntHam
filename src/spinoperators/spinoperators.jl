export AbstractBlock, PauliString, Block, x, y, z
export mat

using SparseArrays

import Base: *
import Base: +
import Base: ^
import Base: -
import Base: convert



# NOT USED CURRENTLY, WAY TOO SLOW
function check_for_duplicate_locs(locs::NTuple{L, <:Any}) where {L}
    for index in 1:L-1
        indices = findall(isequal(locs[index]), locs)
        if !isnothing(indices)
            length(indices) > 1 && throw(ArgumentError("locs must not contain duplicate locations!"))
        end
    end 
end 

""" 
    AbstractBlock{S,N}

Abstract type for Pauli strings and blocks. 
The parameter S is for the spin number and N for the number of spins.
"""
abstract type AbstractBlock{S,N} end

""" 
    PauliString{S, N, L} <: AbstractBlock{S,N}
    x(N::Int, locs; S::Union{Rational,Int}=1//2)
    y(N::Int, locs; S::Union{Rational,Int}=1//2)
    z(N::Int, locs; S::Union{Rational,Int}=1//2)

Pauli String of spin number `S` with `N` spins on locs.

`x`, `y` and `z` create Pauli X, Y and Z operators on the `locs`, repspectively.

`locs` can be either a tuple of Int or just a single Int.

Checks wether S is divisble by `1/2`. 

Currently, only one pauli operator per `PauliString` is supported.

!!! note 
    The definition here slightly deviates from the definition found in most literature. 
    A `PauliString` saves its power.

# Example

Constructing a Pauli string with 4 spins, where the pauli Z operator acts on spin 1 and 2 
```jldoctest PauliStrings
julia> using QuantVarEntHam

julia> ps = z(4, (1,2))
Pauli string
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂
```
Squaring it constructs a new `PauliString` object with power of two
```jldoctest PauliStrings
julia> ps^2
Pauli string
Spin 1//2
Number of spins: 4

(Z₁⊗ Z₂)²
```
"""
struct PauliString{S,N,L} <: AbstractBlock{S,N}
    sig::String
    locs::NTuple{L,Int}
    power::Int
    function PauliString{S,N,L}(sig, locs, power) where {S,N,L}
        divisible_by_half(S)
        #check_for_duplicate_locs(locs)  MAYBE FIND FASTER WAY TO CHECK THIS, TAKES LONG TIME
        isa(L, Int) ||
            throw(ArgumentError("wrong type; L must be of type Int$(Sys.WORD_SIZE)"))
        sig ∈ ["X", "Y", "Z"] || throw(ArgumentError("sig must be \"X\", \"Y\" or \"Z\""))
        if maximum(locs) > N || minimum(locs) < 1
            error("Location out of bounds!")
        end
        new{S,N,L}(sig, locs, power)
    end
end

function PauliString(N::Int, sig::String, locs::Int; S::Union{Rational,Int}=1//2)
    PauliString(N, sig, (locs,); S=S)
end

function PauliString(
    N::Int, sig::String, locs::NTuple{L,Int}; S::Union{Rational,Int}=1//2
) where {L}
    ps = PauliString{Rational(S),N,L}(sig, locs, 1)
    return ps
end

"""
    Block{S,N} <: AbstractBlock{S,N}
    Block(prefactors::Vector{<:Real}, pauli_strings::Vector{<:PauliString{S,N}}) where {S,N}

Linear combinations of pauli strings aka a `Block` with the given `prefactors` and `pauli_strings`. 

# Example 

Adding two `PauliString`s with `prefactors` returns a `Block`.

```jldoctest Block
julia> using QuantVarEntHam

julia> ps1 = z(4, (1,2))
Pauli string
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂

julia> ps2 = x(4, (3,4))
Pauli string
Spin 1//2
Number of spins: 4

X₃⊗ X₄

julia> ps1+2*ps2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + 2*X₃⊗ X₄
```

Ultimately, this could have been achieved by using the constructor for `Block`

```jldoctest Block
julia> Block([1,2], [ps1,ps2])
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + 2*X₃⊗ X₄
```

In this example, `1` and `2` would be the `prefactors` of  the repspective `PauliString`s.
"""
struct Block{S,N} <: AbstractBlock{S,N}
    prefactors::Vector{Float64}
    pauli_strings::Vector{PauliString{S,N}}
    function Block{S,N}(prefactors, pauli_strings) where {S,N}
        divisible_by_half(S)
        new{S,N}(prefactors, pauli_strings)
    end
end

function Block(
    prefactors::Vector{<:Real}, pauli_strings::Vector{<:PauliString{S,N}}
) where {S,N}
    @assert length(prefactors) == length(pauli_strings) "The number of prefactors and Pauli strings need to be equal!"
    Block{S,N}(Float64.(prefactors), pauli_strings)
end

function convert(::Type{Block{S,N}}, ps::PauliString{S,N,L}) where {S,N,L}
    Block{S,N}([1.0], [ps])
end

function convert(::Type{Block{S,N}}, block::Block{S,N}) where {S,N}
    return block
end


function x(N::Int, locs::NTuple{L,Int}; S::Union{Rational,Int}=1//2
) where {L}
    ps = PauliString{Rational(S),N,L}("X", locs, 1)
    return ps
end 

function y(N::Int, locs::NTuple{L,Int}; S::Union{Rational,Int}=1//2
) where {L}
    ps = PauliString{Rational(S),N,L}("Y", locs, 1)
    return ps
end 

function z(N::Int, locs::NTuple{L,Int}; S::Union{Rational,Int}=1//2
) where {L}
    ps = PauliString{Rational(S),N,L}("Z", locs, 1)
    return ps
end 

x(N::Int, loc::Int; S::Union{Rational,Int}=1//2) = x(N, (loc,), S=S)
y(N::Int, loc::Int; S::Union{Rational,Int}=1//2) = y(N, (loc,), S=S)
z(N::Int, loc::Int; S::Union{Rational,Int}=1//2) = z(N, (loc,), S=S)

include("algebra.jl")
include("matrices.jl")
include("printing.jl")
