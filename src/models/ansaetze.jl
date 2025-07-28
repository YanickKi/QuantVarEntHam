"""
    AbstractAnsatz{S,N_A}

Abstract type for the variational ansätze.
"""
abstract type AbstractAnsatz{S,N_A} end

"""
    H_A_BW{S,N_A} <: AbstractAnsatz{S,N_A}
    H_A_BW(model::AbstractModel{S,N_A}, r_max::Int=1) where {S,N_A}

BW like ansatz for a given `model`.

The maximum range of long range corrections is given by `r_max` (`r_max=1` for no corrections).

# Example 
BW Ansatz for the [`TFIM`](@ref) for `N=8`, `N_A`, `Γ =1` and `J=-1`.
```jlcon
julia> model = TFIM(8,4,1);
Diagonalizing the Hamitlonian via exact diagonalization for constructing the ground state density matrix

julia> ansatz = H_A_BW(model)
H_A_BW
Number of blocks: 4

Spin 1//2
Number of spins: 4

Block 1: 
	-X₁ - 0.5*Z₁⊗ Z₂

Block 2: 
	-X₂ - 0.5*Z₁⊗ Z₂ - 0.5*Z₂⊗ Z₃

Block 3: 
	-X₃ - 0.5*Z₂⊗ Z₃ - 0.5*Z₃⊗ Z₄

Block 4: 
	-X₄ - 0.5*Z₃⊗ Z₄
```

"""
struct H_A_BW{S,N_A} <: AbstractAnsatz{S,N_A}
    r_max::Int
    blocks::Vector{Block{S,N_A}}
end

function H_A_BW(model::AbstractModel{S,N_A}, r_max::Int=1) where {S,N_A}
    N = model.N
    periodic = model.periodic
    J = model.J

    if 2*N_A != N && periodic == false
        @warn "Be aware: The Bisognano-Wichmann theorem for the case of open boundary conditions is only valid for N = 2*N_A i.e. for a half plane!"
    end

    blocks = Block{S,N_A}[]

    H_A_BW_wo_corrections!(blocks, model)

    blocks *= J

    if r_max > 1
        corrections!(blocks, model, r_max)
    end
    return H_A_BW(r_max, blocks)
end

"""
    H_A_BWV{S,N_A} <: AbstractAnsatz{S,N_A}
    H_A_BWV(model::AbstractModel{S,N_A}, r_max::Int=1) where {S,N_A}

BW violating Ansatz for a given `model`.

The maximum range of long range corrections is given by `r_max` (`r_max=1` for no corrections).

# Example  
BW violating Ansatz for the [`TFIM`](@ref) for `N=8`, `N_A`, `Γ =1` and `J=-1`.
```jlcon
julia> model = TFIM(8,4,1);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix

julia> ansatz = H_A_BWV(model)
H_A_BWV
Number of blocks: 7

Spin 1//2
Number of spins: 4

Block 1: 
	-X₁

Block 2: 
	-Z₁⊗ Z₂

Block 3: 
	-X₂

Block 4: 
	-Z₂⊗ Z₃

Block 5: 
	-X₃

Block 6: 
	-Z₃⊗ Z₄

Block 7: 
	-X₄
```
"""
struct H_A_BWV{S,N_A} <: AbstractAnsatz{S,N_A}
    r_max::Int
    blocks::Vector{Block{S,N_A}}
end

function H_A_BWV(model::AbstractModel{S,N_A}, r_max::Int=1) where {S,N_A}
    J = model.J

    blocks = Block{S,N_A}[]

    H_A_BWV_wo_corrections!(blocks, model)

    blocks *= J

    if r_max > 1
        corrections!(blocks, model, r_max)
    end

    return H_A_BWV(r_max, blocks)
end

function H_A_BW_wo_corrections!(
    blocks::Vector{Block{S,N_A}}, model::AbstractModel{S,N_A}
) where {S,N_A}
    for i in 1:N_A
        push!(blocks, hi(model, i))
    end
end

function corrections!(
    blocks::Vector{Block{S,N_A}}, model::AbstractModel{S,N_A}, r_max::Int
) where {S,N_A}
    for r in 2:r_max
        for i in 1:(N_A - r)
            correction!(blocks, model, i, r)
        end
    end
end
