"""
    TFIM

Concrete type of [`AbstractModel`](@ref), containing the modeltings for the TFIM

!!! tip
    Use the constructor [`TFIM`](@ref) to instantiate this struct since the type for `observables` is automatically inferred
    and the default values in [`TFIM`](@ref) are highly recommended.

# Fields 
- see [`AbstractModel`](@ref)
- `Γ::Real`: transverse field strength 
"""
@with_kw struct TFIM <: AbstractModel
    N::Int
    N_A::Int
    Γ::Float64
    J::Float64
    S::Rational
    r_max::Int
    periodic::Bool
    ρ_A::Matrix{ComplexF64} 
end

"""
    TFIM(N::Int, N_A::Int, Γ::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool=false,
    J::Real=-1, ρ_A::AbstractMatrix=get_ρ_A(H_TFIM(N, Γ, periodic = periodic, J=J, S=S),  N-N_A+1:N, N))

Convenient constructor for [`TFIM`](@ref) containing modeltings for the TFIM

# Required Arguments
- `N::Int`: number of sites in the composite system.
- `Γ::Real`: transversal field strength 
- `N_A::Int`: number of sites in subsystem A.

# Keyword arguments
- `S::Union{Int64, Rational} = 1//2`: spin number.
- `J::Real=-1`: global prefactor in the Hamiltonian.
- `r_max::Int=1`: maximum range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) r_max = N_A-1 is maximally possible.
- `periodic::Bool=false`: boundary conditions for the system Hamiltonian, false for open and true for periodic boundary conditions.
- `ρ_A::AbstractMatrix=get_ρ_A(H_TFIM(N, Γ, periodic = periodic, J=J),  N-N_A+1:N, N)`: reduced density matrix of ground state of the composite system on subsystem A, by default the subsystem is on the right border.

# Recommendations
- use only Z-gates (or composition of these) as observables since these are diagonal thus save computation time the most. 
"""
function TFIM(N::Int, N_A::Int, Γ::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool=false,
    J::Real=-1, ρ_A::AbstractMatrix=get_ρ_A(H_TFIM(N, Γ, periodic = periodic, J=J, S=S),  N-N_A+1:N, N))

    return TFIM(
        N = N, N_A = N_A,
        S = S,
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

    return J*(ising_term + Γ*transversal_term) 
end 


function hi(model::TFIM, i::Int)
    @unpack N_A, Γ, S = model
     
    hi = Γ * repeat(N_A, X, i, S=S)
    
    if i > 1
        hi += 1/2 * repeat(N_A, Z, i-1, S=S) * repeat(N_A, Z, i, S=S)
    end
    if i < N_A
        hi += 1/2 * repeat(N_A, Z, i, S=S) * repeat(N_A, Z, i+1, S=S)
    end    
    return hi
end

function correction!(blks::Vector{<:AbstractMatrix}, model::TFIM, i::Int, r::Int)
    @unpack N_A, S = model   
    push!(blks, repeat(N_A,Z,(i,i+r), S=S))
end

function H_A_notBW_wo_corrections!(blks::Vector{<:AbstractMatrix}, model::TFIM)
    @unpack N_A, Γ, S = model
    
    push!(blks, Γ*repeat(N_A, X, 1, S=S))

    for i ∈ 1:N_A-1 
        push!(blks, repeat(N_A,Z,(i,i+1), S=S))
        push!(blks, Γ*repeat(N_A, X, i+1, S=S))
    end 
    
end
