export H_TFIM_blocks

"""
    TFIM <: AbstractModel

Object containing the settings for the TFIM.

# Fields 
- see [`AbstractModel`](@ref)
- `Γ::Real`: transverse field strength 
"""
@with_kw struct TFIM{S,N_A} <: AbstractModel{S,N_A}
    N::Int
    Γ::Float64
    J::Float64
    r_max::Int
    periodic::Bool
    ρ_A::Matrix{ComplexF64} 
end

"""
    TFIM(N::Int, N_A::Int, Γ::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool=false,
    J::Real=-1, ρ_A::AbstractMatrix=get_ρ_A(H_TFIM(N, Γ, periodic = periodic, J=J, S=S),  N-N_A+1:N, N))

Convenient constructor for [`TFIM`](@ref) containing settings for the TFIM.
The default values are often used and the density matrix is automatically constructed.

# Required Arguments
- `N`: number of sites in the composite system
- `Γ`: transversal field strength 
- `N_A`: number of sites in subsystem A

# Keyword arguments
- `S`: spin number.
- `J`: global prefactor in the Hamiltonian.
- `r_max`: maximum range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) `r_max = N_A-1` is maximally possible.
- `periodic`: boundary conditions for the system Hamiltonian, false for open and true for periodic boundary conditions.
- `ρ_A`: reduced density matrix of ground state of the composite system on subsystem A, by default the subsystem is on the right border.
"""
function TFIM(N::Int, N_A::Int, Γ::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool=false,
    J::Real=-1, ρ_A::AbstractMatrix=get_rhoA(H_TFIM(N, Γ, periodic = periodic, J=J, S=S),  N-N_A+1:N, N))

    return TFIM{S, N_A}(
        N = N,
        Γ = Γ, J = J, 
        r_max = r_max, periodic = periodic,
        ρ_A = ρ_A
    ) 
end 

"""
    H_TFIM(N::Int, Γ::Real; J::Real = -1, periodic::Bool=false, S::Union{Int64, Rational} = 1//2)

Return the TFIM Hamiltonian ``H= J(\\sum_{i=1}^{N-1} Z_{i}Z_{i+1} + Γ \\sum_{i=1}^{N} X_i)`` with `N` sites, transversal field `Γ` and 
prefactor `J` as a sparse matrix.

Set `periodic` as true for PBC or as false for OBC and 
`S` for the spin number.
"""
function H_TFIM(N::Int, Γ::Real; J::Real = -1, periodic::Bool=false, S::Union{Int64, Rational} = 1//2)

    ising_term = map(1:(periodic ? N : N-1)) do i
        repeat(N,Z,(i,i%N+1), S=S)
    end |> sum

    transversal_term = map(1:N) do i
        repeat(N, X, i, S=S)
    end |> sum

    return Float64(J)*(ising_term + Float64(Γ)*transversal_term) 
end 


function hi(model::TFIM{S,N_A}, i::Int) where {S,N_A}
    @unpack Γ = model
     
    hi = Γ * PauliString(N_A, "X", i, S=S)
    
    if i > 1
        hi += 1/2 * PauliString(N_A, "Z", (i-1, i), S=S)
    end
    if i < N_A
        hi += 1/2 * PauliString(N_A, "Z", (i,i+1), S=S)
    end    
    return hi
end

function correction!(blks::Vector{<:Block}, ::TFIM{S,N_A}, i::Int, r::Int) where {S,N_A}
    push!(blks, PauliString(N_A,"Z",(i,i+r), S=S))
end

function H_A_notBW_wo_corrections!(blks::Vector{<:Block}, model::TFIM{S,N_A}) where {S, N_A}
    @unpack Γ = model
    
    push!(blks, Γ*PauliString(N_A, "X", 1, S=S))

    for i ∈ 1:N_A-1 
        push!(blks, PauliString(N_A,"Z",(i,i+1), S=S))
        push!(blks, Γ*PauliString(N_A, "X", i+1, S=S))
    end 
    
end
