export AbstractBlock, PauliString, Block
export mat

using SparseArrays

import Base: *
import Base: +
import Base: ^
import Base: -
import Base: convert

""" 
    AbstractBlock{S,N}

Abstract type for Pauli strings and blocks. 
The parameter S is for the spin number and N for the number of spins.
"""
abstract type AbstractBlock{S,N} end

""" 
    PauliString{S, N, L} <: AbstractBlock{S,N}
    PauliString(N, sig::String, locs::NTuple{L,Int}; S::Union{Rational, Int} = 1//2) where {L}
    PauliString(N, sig::String, locs::Int; S::Union{Rational, Int} = 1//2)

Pauli String of spin number `S` with `N` spins on `L` locations.

Checks wether `sig` is either `"X"`, `"Y"` or `"Z"` and 
if `S`, `N` and `L` are of type `Rational{Int}` and divisible by `1/2`, `Int`, and `Int`, respectively (error thrown if not). 

Currently, only one pauli operator per `PauliString` is supported.

!!! note 
    The definition here slightly deviates from the definition found in most literature. 
    A `PauliString` saves its power.

# Example

Constructing a Pauli string with 4 spins, where the pauli Z operator acts on spin 1 and 2 
```jldoctest PauliStrings
julia> using QuantVarEntHam

julia> ps1 = PauliString(4,"Z", (1,2))
Pauli string
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂
```
Squaring it constructs a new `PauliString` object with power of two
```jldoctest PauliStrings
julia> ps1^2
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
        check_S_N_type(S, N)
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

julia> ps1 = PauliString(4,"Z", (1,2))
Pauli string
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂

julia> ps2 = PauliString(4,"X", (3,4))
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
        check_S_N_type(S, N)
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

include("algebra.jl")
include("matrices.jl")
include("printing.jl")
