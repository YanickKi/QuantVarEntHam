using KrylovKit: eigsolve
using LinearAlgebra
using SparseArrays

export AbstractModel
export TFIM, XXZ, Pollmann
export H_XXZ, H_TFIM, H_pollmann
export AbstractAnsatz, H_A_BW, H_A_BWV
export getrho_A, getblocks, getr_max
export H_A

"""
    AbstractModel{S,N_A}

Abstract type to dispatch on the concrete types for the correct variational Ansätze.

Spin number `S` and number of spins in the subsystem `N_A`.

All physical `AbstractModel`s have their own concrete (e.g. [`TFIM`](@ref)).
In general the concrete types will have the same fields besides the model specific Hamiltonian parameters, which are: 
 
- `N::Int`: number of sites in composite system.
- `J::Float64`: global prefactor in Hamiltonian
- `periodic::Bool`: boundary conditions for the system Hamiltonian, false for open and true for periodic boundary conditions, obsolete if an own reduced density matrix ρ_A is provided.
- `ρ_A::Matrix{ComplexF64}`: reduced density matrix of ground state of the composite system on subsystem A.
"""
abstract type AbstractModel{S,N_A} end


function rho_A(H::AbstractBlock{S,N}, N_A::Int) where {S,N} 

    H_mat = mat(H)
    d = (Int64(2*S+1))^(N)
    if d > 1024
        println("Diagonalizing the Hamiltonian via Krylov subspace method for constructing the ground state density matrix")
        _, vectors = eigsolve(H_mat, 1 ,:SR, ishermitian=true, tol = 1e-16)
        ρ_A = partial_trace_pure_state(vectors[1], (Int64(2*S+1))^(N-N_A), (Int64(2*S+1))^(N_A), trace_subsystem = :B) 
    return ρ_A
    else 
        println("Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix")
        vectors = eigvecs(Hermitian(Matrix(H_mat)))
        ρ_A = partial_trace_pure_state(vectors[:,1],(Int64(2*S+1))^(N-N_A), (Int64(2*S+1))^(N_A), trace_subsystem = :B) 
        return ρ_A
    end
end 


include("xxz.jl")
include("tfim.jl")
include("pollmann.jl")
include("ansaetze.jl")
include("utils.jl")
include("printing.jl")
include("getter.jl")