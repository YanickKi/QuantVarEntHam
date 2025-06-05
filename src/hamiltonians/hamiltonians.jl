using KrylovKit: eigsolve
using LinearAlgebra
using SparseArrays


"""
    Settings{M<:AbstractMatrix}

Abstract type to dispatch on the concrete types for the correct variational Ansätze.

The parametric type `M<:AbstractBlock` is for determining the most efficient representation of observables.
All physical models have their own concrete type of `Settings` (e.g. `Settings_TFIM`).
In general the concrete types will have the same fields besides the model specific Hamiltonian parameters, which are: 
 
- `N::Int`: number of sites in composite system.
- `N_A::Int`: number of sites in subsystem A.
- `J::Float64`: global prefactor in Hamiltonian
- `T_max::Float64`: maximum time for evolving the observables i.e. maximum integration time.
- `S::Rational`: spin number
- `r_max::Int`: range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) r_max = N_A-1 corresponds to maximum order.
- `periodic::Bool`: boundary conditions for the system Hamiltonian, false for open and true for periodic boundary conditions, obsolete if an own reduced density matrix ρ_A is provided.
- `ρ_A::Matrix{ComplexF64}`: reduced density matrix of ground state of the composite system on subsystem A.
- `meas0::Vector{Float64}`: expectation values of `observables` at time ``t=0``.
- `observables::Vector{M}`: matrix representations of the observables
"""
abstract type Settings{M<:AbstractMatrix} end

"""
    get_rhoA(H::AbstractMatrix, A::AbstractVector{Int}, N::Int) 

Return the reduced density matrix as either a real or complex matrix of the ground state of Hamiltonian `H` for `N` sites on subsystem A.

For a Hilbert space dimension of more than 1024, the [krylov subspace method from KrylovKit](https://jutho.github.io/KrylovKit.jl/stable/man/eig/#KrylovKit.eigsolve) is used for 
extracting the ground state, exact diagonalization otherwise.
"""

function get_rhoA(H::AbstractMatrix, A::AbstractVector{Int}, N::Int; S::Union{Rational, Int} = 1//2) 
    d = (Int64(2*S+1))^(N)
    N_A = length(A)
    if d > 1024
        println("Diagonalizing the Hamitlonian via Krylov subspace method for constructing the ground state density matrix")
        _, vectors = eigsolve(H, 1 ,:SR, ishermitian=true, tol = 1e-16)
        ρ_A = partial_trace_pure_state(vectors[1], (Int64(2*S+1))^(N-N_A), (Int64(2*S+1))^(N_A), trace_subsystem = :B) 
    return ρ_A
    else 
        println("Diagonalizing the Hamitlonian via exact diagonalization for constructing the ground state density matrix")
        vectors = eigvecs(Hermitian(Matrix(H)))
        ρ_A = partial_trace_pure_state(vectors[:,1],(Int64(2*S+1))^(N-N_A), (Int64(2*S+1))^(N_A), trace_subsystem = :B) 
        return ρ_A
    end
end 


function partial_trace_pure_state(ψ::Vector{<:Number}, dimA::Int, dimB::Int; trace_subsystem::Symbol = :B)
    ψ_tensor = reshape(ψ, (dimA, dimB)) 

    if trace_subsystem == :B
        return ψ_tensor * ψ_tensor' 
    elseif trace_subsystem == :A
        return transpose(ψ_tensor) * conj(ψ_tensor)
    else
        error("trace_subsystem must be :A or :B")
    end
end

function expect(Op::AbstractMatrix, ρ::AbstractMatrix)
    E = tr(Op*ρ)
    if abs(imag(E)) > eps(Float64)
        error("Imaginary part too large for expectation value!")
    end 
    return real(E)
end 


"""
    H_A_BW(set::Settings) 

Return a vector of the blocks, which are complex dense matrices.

The variational Ansatz follows the Bisognano-Wichmann-theorem.
This function calls lower level functions which dispatch on the concrete subtypes of the abstract type [`Settings`](@ref) to get the correct variational Ansatz for the corresponding model.


# Example 
`H_A_BW(set::Settings_TFIM)` returns the blocks of the variational Ansatz for the TFIM following the Bisognano-Wichmann-theorem.
"""
function H_A_BW(set::Settings)
    @unpack N, N_A, r_max, periodic, J = set
    
    if 2*N_A != N && periodic == false 
        @warn "Be aware: The Bisognano-Wichmann theorem for the case of open boundary conditions is only valid for N = 2*N_A i.e. for a half plane!" 
    end 
    
    blocks = Matrix{ComplexF64}[]
    
    H_A_BW_wo_corrections!(blocks, set)
    
    if r_max > 1
        corrections!(blocks, set)
    end 
    return J*blocks

end 

"""
    H_A_not_BW(set::Settings) 

Return a vector of the blocks, which are complex dense matrices.

The variational Ansatz does not follow the Bisognano-Wichmann-theorem.
This function calls lower level functions which dispatch on the concrete subtypes of the abstract type [`Settings`](@ref) to get the correct variational Ansatz for the corresponding model.

# Example 
`H_A_not_BW(set::Settings_TFIM)` returns the blocks of  the variational Ansatz for the TFIM not following the Bisognano-Wichmann-theorem.
"""
function H_A_not_BW(set::Settings)
    @unpack N_A, r_max, J = set
    
    blocks = Matrix{ComplexF64}[]
    
    H_A_notBW_wo_corrections!(blocks, set)
    
    if r_max > 1
        corrections!(blocks, set)
    end 

    return J*blocks
end 

function H_A_BW_wo_corrections!(blocks::Vector{<:AbstractMatrix}, set::Settings)
    @unpack N_A = set
    for i in 1:N_A 
        push!(blocks, hi(i, set))
    end     
end


function corrections!(blocks::Vector{<:AbstractMatrix}, set::Settings)
    @unpack N_A, r_max = set
    for r in 2:r_max
        for i in 1:N_A-r
            correction!(blocks, i, r, set)
        end
    end 
end 

include("spinoperators.jl")
include("xxz.jl")
include("tfim.jl")
include("pollmann.jl")
#include("toric_code.jl")
#include("kitaev.jl")
#=


include("toric_code.jl")
include("kitaev.jl")



function map_to_subsystem_chain!(objects, A)


    for obj in objects
        for i in eachindex(obj)
            obj[i] = findfirst(==(obj[i]), A)
        end 
    end

end 


function H_A_BW_bad(set::Settings_toric)
    @unpack Nx,Ny, A ,J = set 
    
    blocks = Vector{AbstractBlock}(undef, 0)

    N_A = length(A)

    for qbit in 1:N_A-1
        push!(blocks, repeat(N_A, Z, (qbit, qbit+1)))
    end 
    for qbit in 1:N_A-1
        push!(blocks, repeat(N_A, X, (qbit, qbit+1)))
    end 

    #for plaquette in plaquets_in_A
    #    push!(blocks, repeat(N_A, Z, plaquette))
    #end

    blocks *= -J 
    
    return H_A_Var(blocks, mat.(blocks))
end

function H_A_BW(set::Settings_kitaev)
    @unpack Jz, Jx, Jy = set
    blocks = Vector{AbstractBlock}(undef, 0)

    zterm = -Jz * [
        repeat(6, Z, (2,4)),
        repeat(6, Z, (3,5))
    ]

    xterm = -Jx * [
        repeat(6, X, (1,2)),
        repeat(6, X, (5,6))
        ]

    yterm = -Jy * [
        repeat(6, Y, (1,3)),
        repeat(6, Y, (4,6))
        ]

    for i in 1:2
        push!(blocks, zterm[i])
    end 
    for i in 1:2
        push!(blocks, xterm[i])
    end 
    for i in 1:2
        push!(blocks, yterm[i])
    end

    return H_A_Var(blocks, mat.(blocks)) 
    
end 
function H_A_BW(set::Settings_toric)
    @unpack Nx,Ny, A ,J = set 

    stars, plaquets = make_all_objects(Nx,Ny)

    blocks = Vector{AbstractBlock}(undef, 0)

    N_obj = length(stars)
    stars_in_A = Vector{Int64}[]
    plaquets_in_A = Vector{Int64}[]
    iscomplete = true 

    for i in 1:N_obj    
        for qbit in stars[i]
            iscomplete =  qbit ∈ A
            if iscomplete == false
                break 
            end 
        end 
        if iscomplete == true 
            push!(stars_in_A, stars[i])
        end 
    end
    iscomplete = true 


    for i in 1:N_obj    
        for qbit in plaquets[i]
            iscomplete =  qbit ∈ A
            if iscomplete == false
                break 
            end 
        end 
        if iscomplete == true 
            push!(plaquets_in_A, plaquets[i])
        end 
    end 

    map_to_subsystem_chain!(stars_in_A,A)
    map_to_subsystem_chain!(plaquets_in_A, A)
    
    N_A = length(A)

    for star in stars_in_A
        push!(blocks, repeat(N_A, X, star))
    end 

    for plaquette in plaquets_in_A
        push!(blocks, repeat(N_A, Z, plaquette))
    end

    
    N_A = length(A)

    for qbit in 1:N_A
        push!(blocks, put(N_A, qbit => X))
    end 

    blocks *= -J 


    return H_A_Var(blocks, mat.(blocks))
end

=#