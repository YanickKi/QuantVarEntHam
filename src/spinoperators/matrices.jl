function X(S::Union{Rational,Int})
    d = Int64(2*S+1)

    multi = S*(S+1)

    Sx = spzeros(ComplexF64, d, d)

    for m in (-S):S
        _m = Int64(S-m+1)
        if m < S
            m´ = m+1
            _m´ = Int64(S-m´+1)
            Sx[_m, _m´] = sqrt(multi - m´*(m´-1))
        end
        if m > -S
            m´ = m-1
            _m´ = Int64(S-m´+1)
            Sx[_m, _m´] = sqrt(multi - m´*(m´+1))
        end
    end

    return Sx
end

function Y(S::Union{Rational,Int})
    d = Int64(2*S+1)

    multi = S*(S+1)

    Sy = spzeros(ComplexF64, d, d)

    for m in (-S):S
        _m = Int64(S-m+1)
        if m < S
            m´ = m+1
            _m´ = Int64(S-m´+1)
            Sy[_m, _m´] = - sqrt(multi - m´*(m´-1))
        end
        if m > -S
            m´ = m-1
            _m´ = Int64(S-m´+1)
            Sy[_m, _m´] = sqrt(multi - m´*(m´+1))
        end
    end

    Sy ./= 1im

    return Sy
end

function Z(S::Union{Rational,Int})
    d = Int64(2*S+1)

    Sz = spzeros(ComplexF64, d, d)

    for m in (-S):S
        _m = Int64(m+S+1)
        Sz[_m, _m] = - Float64(2*m)
    end

    return Sz
end

function repeat(N, sig::Function, locs::Int; S::Union{Rational,Int}=1//2)
    repeat(N, sig, (locs,); S=S)
end

function repeat(
    N, sig::Function, locs::NTuple{l,Int}; S::Union{Rational,Int}=1//2
) where {l}
    d = Int64(2*S+1)

    M = spzeros(ComplexF64, d, d)

    Identity = sparse((1.0+0.0im)*I, d, d)

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
    return M
end

function repeat(N, sig::String, locs::NTuple{L,Int}; S::Union{Rational,Int}=1//2) where {L}
    if sig == "X"
        return repeat(N, X, locs; S=S)
    end
    if sig == "Y"
        return repeat(N, Y, locs; S=S)
    end

    return repeat(N, Z, locs; S=S)
end

"""
    mat(block::AbstractBlock)

Return a complex, sparse matrix of an [`AbstractBlock`](@ref).

!!! note 
    The matrix elements look very similar to the well known matrix elements of the spin operators in the z-basis but they differ 
    by a factor of 2 s.t. the matrices of a [`PauliString`](@ref) for a single spin half will give the well known pauli matrices.  

The matrix elements are given by
```math   
\\begin{aligned} 
\\langle m' | X | m\\rangle &=
(\\delta_{m',m+1}+\\delta_{m'+1,m})
\\sqrt{S(S+1)-m'm}\\\\
\\langle m' | Y | m\\rangle &=
i(\\delta_{m'+1,m}-\\delta_{m',m+1})
\\sqrt{S(S+1)-m'm}\\\\
\\langle m' | Z | m\\rangle &=
2\\delta_{m',m}m
\\end{aligned}
```

# Example 

Matrix representation of the [`PauliString`](@ref)s for a single spin half.

```jldoctest
julia> using QuantVarEntHam

julia> paulistrings = [x(1,1),y(1,1),z(1,1)];

julia> sparse_mat = mat.(paulistrings)
3-element Vector{SparseArrays.SparseMatrixCSC{ComplexF64, Int64}}:
 sparse([2, 1], [1, 2], ComplexF64[1.0 + 0.0im, 1.0 + 0.0im], 2, 2)
 sparse([2, 1], [1, 2], ComplexF64[0.0 + 1.0im, 0.0 - 1.0im], 2, 2)
 sparse([1, 2], [1, 2], ComplexF64[1.0 + 0.0im, -1.0 + 0.0im], 2, 2)

julia> dense_mat = Matrix.(sparse_mat)
3-element Vector{Matrix{ComplexF64}}:
 [0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im]
 [0.0 + 0.0im 0.0 - 1.0im; 0.0 + 1.0im 0.0 + 0.0im]
 [1.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im -1.0 + 0.0im]
```

There we see the pauli matrices.
"""
mat(block::AbstractBlock) = mat(block)

function mat(block::Block{S,N}) where {S,N}
    d = Int((2*S+1)^N)

    M = spzeros(ComplexF64, d, d)

    for index in eachindex(block.pauli_strings)
        M .+=
            block.prefactors[index] .*
            (repeat(
                N, block.pauli_strings[index].sig, block.pauli_strings[index].locs; S=S
            ))^block.pauli_strings[index].power
    end
    return M
end

function mat(ps::PauliString{S,N}) where {S,N}
    return repeat(N, ps.sig, ps.locs; S=S)^(ps.power)
end