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

!!! note
    All constructors for the concrete types will have the keyword argument `ϵ_reg`.
    Due to numerical imprecisions, the resulting exact reduced density matrices will often times not be positive semidefinite, which will in turn 
    have a huge impact when computing the exact Entanglement Hamiltonian via ``H_\\text{A} = - \\log (ρ_\\text{A})`` and its corresponding 
    Entanglement spectrum.
    If `ϵ_reg` is larger than zero and the smallest eigenvalue of the intermediate exact reduced density matrix is smaller than zero, the reduced density matrix will be regularized by 
    adding `ϵ_reg` and the smallest eigenvalue of the intermediate exact reduced density matrix to the diagonal, i.e. performing the computation 
    ```math
    \\rho_\\text{A} = \\rho'_\\text{A} + (\\epsilon_\\text{reg} + \\text{min} (\\text{eigenvalues} (\\rho'_\\text{A} )))I,
    ```
    where ``I`` is the unit matrix, ``\\rho'_\\text{A}`` the intermediate exact reduced density matrix after the exact diagonalization of the system Hamiltonian and 
    tracing out subsystem B and ``\\rho_\\text{A}`` the returned exact reduced density matrix, which will be actually saved and used for the upcoming computations.
"""
abstract type AbstractModel{S,N_A} end

function regulate!(ρ_A::Matrix, ϵ_reg::Real)
    vals, _ = eigen(Hermitian(ρ_A))
    lowest_val = minimum(vals)
    if lowest_val < 0. 
        sum_correction = ϵ_reg + lowest_val
         println(
            "Regularizing aka making the exact reduced density matrix positive semidefinite with ϵ_reg=$(ϵ_reg)",
        )
        ρ_A[diagind(ρ_A)] .+= sum_correction
    end 
    ρ_A ./= tr(ρ_A)
    return ρ_A
end 

function rho_A(H::AbstractBlock{S,N}, N_A::Int, ϵ_reg::Real=1e-8) where {S,N}
    H_mat = mat(H)
    d = (Int64(2*S+1))^(N)
    if d > 1024
        println(
            "Diagonalizing the Hamiltonian via Krylov subspace method for constructing the ground state density matrix",
        )
        _, vectors = eigsolve(H_mat, 1, :SR; ishermitian=true, tol=1e-16)
        ρ_A = partial_trace_pure_state(
            vectors[1], (Int64(2*S+1))^(N-N_A), (Int64(2*S+1))^(N_A); trace_subsystem=:B
        ) 
        ϵ_reg > 0. && regulate!(ρ_A, ϵ_reg) 
        return ρ_A
    else
        println(
            "Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix",
        )
        vectors = eigvecs(Hermitian(Matrix(H_mat)))
        ρ_A = partial_trace_pure_state(
            vectors[:, 1], (Int64(2*S+1))^(N-N_A), (Int64(2*S+1))^(N_A); trace_subsystem=:B
        )
        ϵ_reg > 0. && regulate!(ρ_A, ϵ_reg) 
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
