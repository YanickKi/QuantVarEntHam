abstract type Settings{T<:AbstractBlock,S<:AbstractMatrix} end

using KrylovKit: eigsolve
using LinearAlgebra


"""
    get_rhoA(H::AbstractBlock, A::AbstractRange, N::Int) 

Return the reduced density matrix on subsystem A.
For composite systems consisting of more than 10 sites, the krylov subspace method from KrylovKit is used for 
extracting the ground state.
"""
function get_rhoA(H::AbstractBlock, A::AbstractRange, N::Int) 
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

struct H_A_Var
    blocks::Vector{AbstractBlock}
    matrices::Vector{Matrix{ComplexF64}}
end 

function H_A_BW(set::Settings) 
    @unpack N, N_A, r_max, periodic = set
    
    if 2*N_A != N && periodic == false 
        @warn "Be aware: The Bisognano-Wichmann theorem for the case of open boundary conditions is only valid for N = 2*N_A i.e. for a half plane!" 
    end 
    
    blks = Vector{AbstractBlock}(undef, 0)
    
    H_A_BW_wo_corrections!(blks, set)
    
    if r_max > 1
        corrections!(blks, set)
    end 
    return H_A_Var(blks, mat.(blks))

end 

function H_A_not_BW(set::Settings) 
    @unpack N_A, r_max = set
    
    blks = Vector{AbstractBlock}(undef, 0)
    
    H_A_notBW_wo_corrections!(blks, set)
    
    if r_max > 1
        corrections!(blks, set)
    end 

    return H_A_Var(blks, mat.(blks))

end 


function H_A_BW_wo_corrections!(blks::Vector{AbstractBlock}, set::Settings)
    @unpack N_A = set
    for i in 1:N_A 
        push!(blks, hi(i, set))
    end     
end


function corrections!(blks::Vector{AbstractBlock}, set::Settings)
    @unpack N_A, r_max = set
    for r in 2:r_max
        for i in 1:N_A-r
            correction!(blks, i, r, set)
        end
    end 
end 


include("xxz.jl")
include("tfim.jl")
