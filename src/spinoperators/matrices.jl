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

    Sz = Diagonal{ComplexF64}(zeros(d))

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
    return sparse(M)
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

Return a complex, sparse matrix of a subtype of [`AbstractBlock`](@ref) (can be either a [`PauliString`](@ref) or a [`Block`](@ref)).
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
