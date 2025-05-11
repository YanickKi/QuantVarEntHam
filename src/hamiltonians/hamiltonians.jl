using KrylovKit: eigsolve
using LinearAlgebra
using QuantumOpticsBase
using QuantumInterface
using SparseArrays


"""
    Settings{T<:AbstractBlock,S<:AbstractMatrix}

Abstract type to dispatch on the concrete types for the correct variational Ansätze.

The parametric type `T<:AbstractBlock` is introduced for determining the correct concrete type of the Yao Blocks, while 
`S<:AbstractMatrix` is needed to determine the concrete type of the matrix representation of the Yao Blocks to prevent 
working with full complex dense matrices.
"""
abstract type Settings{S<:AbstractMatrix} end

"""
    get_rhoA(H::AbstractBlock, A::AbstractVector{Int}, N::Int) 

Return the reduced density matrix of type `YaoAPI.DensityMatrix` of the ground state of Hamitlonian H for N sites on subsystem A.

For composite systems consisting of more than 10 sites, the [krylov subspace method from KrylovKit](https://jutho.github.io/KrylovKit.jl/stable/man/eig/#KrylovKit.eigsolve) is used for 
extracting the ground state.
"""

function get_rhoA(H, A::AbstractVector{Int}, N::Int; S::Union{Rational, Int} = 1//2) 
    d = (Int64(2*S+1))^(N)
    N_A = length(A)
    if d > 1024
        println("Diagonalizing the Hamitlonian via Krylov subspace method for constructing the ground state density matrix")
        values, vectors = eigsolve(H, 1 ,:SR, ishermitian=true, tol = 1e-16)
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
    H_A_Var
    
Struct to save the Yao AbstractBlocks and its matrix representation throughout the optimization.

The optimizing is only done with the matrices in `matrices`. `blocks` is only saved for convenient utiliy functions such as 
printing the Entanglement Hamiltonian.

# Fields
- `blocks::Vector{AbstractBlock}`: Containing each Hamiltonian Block as an Yao AbstractBlock.
- `matrices::Vector{Matrix{ComplexF64}}`: Containing each Hamiltonian Block as a matrix.
"""
struct H_A_Var
    blocks::Vector{Matrix{ComplexF64}}
end 


"""
    H_A_BW(set::Settings) 

Return an instance of type [`H_A_Var`](@ref) containing the Yao Blocks and its matrix representation.

The variational Ansatz follows the Bisognano-Wichmann-theorem.
This function dispatches on the concrete subtypes of the abstract type [`Settings`](@ref) to get the correct variational Ansatz for the corresponding model. 

# Example 
`H_A_BW(set::Settings_TFIM)` returns the variational Ansatz for the TFIM following the Bisognano-Wichmann-theorem.
"""
function H_A_BW(set::Settings)
    @unpack N, N_A, r_max, periodic, signHam = set
    
    if 2*N_A != N && periodic == false 
        @warn "Be aware: The Bisognano-Wichmann theorem for the case of open boundary conditions is only valid for N = 2*N_A i.e. for a half plane!" 
    end 
    
    blks = Matrix{ComplexF64}[]
    
    H_A_BW_wo_corrections!(blks, set)
    
    blks *= signHam


    if r_max > 1
        corrections!(blks, set)
    end 
    return H_A_Var(blks)

end 

"""
    H_A_not_BW(set::Settings) 

Return an instance of type [`H_A_Var`](@ref) containing the Yao Blocks and its matrix representation.

The variational Ansatz does not follow the Bisognano-Wichmann-theorem.
This function dispatches on the concrete subtypes of the abstract type [`Settings`](@ref) to get the correct variational Ansatz for the corresponding model. 

# Example 
`H_A_not_BW(set::Settings_TFIM)` returns the variational Ansatz for the TFIM following not following the Bisognano-Wichmann-theorem.
"""
function H_A_not_BW(set::Settings)
    @unpack N_A, r_max, signHam = set
    
    blks = Matrix{ComplexF64}[]
    
    H_A_notBW_wo_corrections!(blks, set)
    
    blks *= signHam

    if r_max > 1
        corrections!(blks, set)
    end 

    return H_A_Var(blks)
end 

 

function H_A_BW_wo_corrections!(blks::Vector{<:AbstractMatrix}, set::Settings)
    @unpack N_A = set
    for i in 1:N_A 
        push!(blks, hi(i, set))
    end     
end


function corrections!(blks::Vector{<:AbstractMatrix}, set::Settings)
    @unpack N_A, r_max = set
    for r in 2:r_max
        for i in 1:N_A-r
            correction!(blks, i, r, set)
        end
    end 
end 


function H_A_not_BW_I(set::Settings)
    @unpack N, N_A, r_max, periodic, signHam = set
    
    blks = Matrix{ComplexF64}[]

    H_A_notBW_wo_corrections_I!(blks, set)
    
    blks *= signHam


    if r_max > 1
        corrections!(blks, set)
    end 

    return H_A_Var(blks)
end 


function H_A_BW_wo_corrections!(blks::AbstractVector, set::Settings)
    @unpack N_A = set
    for i in 1:N_A 
        push!(blks, hi(i, set))
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
    
    blks = Vector{AbstractBlock}(undef, 0)

    N_A = length(A)

    for qbit in 1:N_A-1
        push!(blks, repeat(N_A, Z, (qbit, qbit+1)))
    end 
    for qbit in 1:N_A-1
        push!(blks, repeat(N_A, X, (qbit, qbit+1)))
    end 

    #for plaquette in plaquets_in_A
    #    push!(blks, repeat(N_A, Z, plaquette))
    #end

    blks *= -J 
    
    return H_A_Var(blks, mat.(blks))
end

function H_A_BW(set::Settings_kitaev)
    @unpack Jz, Jx, Jy = set
    blks = Vector{AbstractBlock}(undef, 0)

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
        push!(blks, zterm[i])
    end 
    for i in 1:2
        push!(blks, xterm[i])
    end 
    for i in 1:2
        push!(blks, yterm[i])
    end

    return H_A_Var(blks, mat.(blks)) 
    
end 
function H_A_BW(set::Settings_toric)
    @unpack Nx,Ny, A ,J = set 

    stars, plaquets = make_all_objects(Nx,Ny)

    blks = Vector{AbstractBlock}(undef, 0)

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
        push!(blks, repeat(N_A, X, star))
    end 

    for plaquette in plaquets_in_A
        push!(blks, repeat(N_A, Z, plaquette))
    end

    
    N_A = length(A)

    for qbit in 1:N_A
        push!(blks, put(N_A, qbit => X))
    end 

    blks *= -J 


    return H_A_Var(blks, mat.(blks))
end

=#