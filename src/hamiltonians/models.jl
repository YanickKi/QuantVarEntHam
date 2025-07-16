using KrylovKit: eigsolve
using LinearAlgebra
using SparseArrays

export AbstractModel
export H_A_BW, H_A_notBW
export TFIM, XXZ, Pollmann
export H_XXZ, H_TFIM, H_pollmann
export getrho_A

"""
    AbstractModel{S,N_A} end

Abstract type to dispatch on the concrete types for the correct variational Ansätze.

Spin number `S` and number of spins in the subsystem `N_A`.

All physical `AbstractModel`s have their own concrete (e.g. [`TFIM`](@ref)).
In general the concrete types will have the same fields besides the model specific Hamiltonian parameters, which are: 
 
- `N::Int`: number of sites in composite system.
- `J::Float64`: global prefactor in Hamiltonian
- `S::Rational`: spin number
- `r_max::Int`: range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) `r_max = N_A-1` corresponds to maximum order.
- `periodic::Bool`: boundary conditions for the system Hamiltonian, false for open and true for periodic boundary conditions, obsolete if an own reduced density matrix ρ_A is provided.
- `ρ_A::Matrix{ComplexF64}`: reduced density matrix of ground state of the composite system on subsystem A.
"""
abstract type AbstractModel{S,N_A} end


function rho_A(H::AbstractBlock{S,N}, N_A::Int) where {S,N} 

    H_mat = mat(H)
    d = (Int64(2*S+1))^(N)
    if d > 1024
        println("Diagonalizing the Hamitlonian via Krylov subspace method for constructing the ground state density matrix")
        _, vectors = eigsolve(H_mat, 1 ,:SR, ishermitian=true, tol = 1e-16)
        ρ_A = partial_trace_pure_state(vectors[1], (Int64(2*S+1))^(N-N_A), (Int64(2*S+1))^(N_A), trace_subsystem = :B) 
    return ρ_A
    else 
        println("Diagonalizing the Hamitlonian via exact diagonalization for constructing the ground state density matrix")
        vectors = eigvecs(Hermitian(Matrix(H_mat)))
        ρ_A = partial_trace_pure_state(vectors[:,1],(Int64(2*S+1))^(N-N_A), (Int64(2*S+1))^(N_A), trace_subsystem = :B) 
        return ρ_A
    end
end 


"""
    H_A_BW(model::AbstractModel) 

Return a vector with the blocks as its entries as `Block`s.

The variational Ansatz follows the Bisognano-Wichmann-theorem.

# Example 
BW Ansatz for the [`TFIM`](@ref) for `N=8`, `N_A`, `Γ =1` and `J=-1`.
```jlcon
julia> model = TFIM(8,4,1);
Diagonalizing the Hamitlonian via exact diagonalization for constructing the ground state density matrix

julia> blocks = H_A_BW(model)
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
function H_A_BW(model::AbstractModel{S,N_A}) where {S,N_A}
    N = model.N
    r_max = model.r_max
    periodic = model.periodic
    J = model.J
    
    if 2*N_A != N && periodic == false 
        @warn "Be aware: The Bisognano-Wichmann theorem for the case of open boundary conditions is only valid for N = 2*N_A i.e. for a half plane!" 
    end 
    
    blocks = Block{S, N_A}[]
    
    H_A_BW_wo_corrections!(blocks, model)
    
    if r_max > 1
        corrections!(blocks, model)
    end 
    return J*blocks
end 

"""
    H_A_notBW(model::AbstractModel) 

Return a vector with the blocks as its entries, which are complex dense matrices.

The variational Ansatz does not follow the Bisognano-Wichmann-theorem.

# Example  
BW violating Ansatz for the [`TFIM`](@ref) for `N=8`, `N_A`, `Γ =1` and `J=-1`.
```jlcon
julia> model = TFIM(8,4,1);
Diagonalizing the Hamitlonian via exact diagonalization for constructing the ground state density matrix

julia> H_A_notBW(model)
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
function H_A_notBW(model::AbstractModel{S,N_A}) where {S,N_A}
    r_max = model.r_max
    J = model.J

    blocks = Block{S, N_A}[]
    
    H_A_notBW_wo_corrections!(blocks, model)
    
    if r_max > 1
        corrections!(blocks, model)
    end 

    return J*blocks
end 


include("xxz.jl")
include("tfim.jl")
include("pollmann.jl")
include("utils.jl")
include("utils_entanglement_hamiltonians.jl")
include("printing.jl")
include("getter.jl")