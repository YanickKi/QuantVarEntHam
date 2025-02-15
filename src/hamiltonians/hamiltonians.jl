using KrylovKit: eigsolve
using LinearAlgebra

"""
    Settings{T<:AbstractBlock,S<:AbstractMatrix}

Abstract type to dispatch on the concrete types for the correct variational Ansätze.

The parametric type `T<:AbstractBlock` is introduced for determining the correct concrete type of the Yao Blocks, while 
`S<:AbstractMatrix` is needed to determine the concrete type of the matrix representation of the Yao Blocks to prevent 
working with full complex dense matrices.
"""
abstract type Settings{T<:AbstractBlock,S<:AbstractMatrix} end



"""
    get_rhoA(H::AbstractBlock, A::AbstractVector{Int}, N::Int) 

Return the reduced density matrix of type `YaoAPI.DensityMatrix` of the ground state of Hamitlonian H for N sites on subsystem A.

For composite systems consisting of more than 10 sites, the [krylov subspace method from KrylovKit](https://jutho.github.io/KrylovKit.jl/stable/man/eig/#KrylovKit.eigsolve) is used for 
extracting the ground state.
"""
function get_rhoA(H::AbstractBlock, A::AbstractVector{Int}, N::Int) 
    if N>10
        println("Diagonalizing the Hamitlonian via Krylov subspace method for constructing the ground state density matrix")
        values, vectors = eigsolve(mat(H) ,1 ,:SR, ishermitian=true, tol = 1e-16)
        rhoA = density_matrix(ArrayReg(vectors[1]), A)
    return rhoA
    else 
        println("Diagonalizing the Hamitlonian via exact diagonalization for constructing the ground state density matrix")
        vectors = eigvecs(Hermitian(Matrix(H)))
        rhoA = density_matrix(ArrayReg(vectors[:,1]), A)
        return rhoA
    end 
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
    blocks::Vector{AbstractBlock}
    matrices::Vector{Matrix{ComplexF64}}
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
    
    blks = Vector{AbstractBlock}(undef, 0)
    
    H_A_BW_wo_corrections!(blks, set)
    
    blks *= signHam

    blks = convert(Vector{AbstractBlock}, blks)

    if r_max > 1
        corrections!(blks, set)
    end 
    return H_A_Var(blks,  mat.(blks))

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
    
    blks = Vector{AbstractBlock}(undef, 0)
    
    H_A_notBW_wo_corrections!(blks, set)
    
    blks *= signHam

    blks = convert(Vector{AbstractBlock}, blks)

    if r_max > 1
        corrections!(blks, set)
    end 

    return H_A_Var(blks, mat.(blks))
end 


function H_A_XYZ(set::Settings)
    @unpack N_A, r_max, signHam = set
    
    blks = Vector{AbstractBlock}(undef, 0)
    
    H_A_XYZ_wo_corrections!(blks, set)
    
    blks *= signHam

    blks = convert(Vector{AbstractBlock}, blks)

    if r_max > 1
        corrections_XYZ!(blks, set)
    end 

    return H_A_Var(blks, mat.(blks))
end 

function H_A_BW_wo_corrections!(blks::Vector{<:AbstractBlock}, set::Settings)
    @unpack N_A = set
    for i in 1:N_A 
        push!(blks, hi(i, set))
    end     
end


function corrections_XYZ!(blks::Vector{<:AbstractBlock}, set::Settings)
    @unpack N_A, r_max = set
    for r in 2:r_max
        for i in 1:N_A-r
            correction_XYZ!(blks, i, r, set)
        end
    end 
end 

function corrections!(blks::Vector{<:AbstractBlock}, set::Settings)
    @unpack N_A, r_max = set
    for r in 2:r_max
        for i in 1:N_A-r
            correction!(blks, i, r, set)
        end
    end 
end 


include("xxz.jl")
include("tfim.jl")
