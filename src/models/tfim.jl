export H_TFIM_blocks

"""
    TFIM{S,N_A} <: AbstractModel{S,N_A}
    TFIM(N::Int, N_A::Int, Γ::Real; S::Union{Int, Rational} = 1//2, periodic::Bool=false, J::Real=-1)

Object containing the settings for the TFIM.

# Fields 
- see [`AbstractModel`](@ref)
- `Γ::Real`: transverse field strength 
"""
struct TFIM{S,N_A} <: AbstractModel{S,N_A}
    N::Int
    Γ::Float64
    J::Float64
    periodic::Bool
    ρ_A::Matrix{ComplexF64}
    function TFIM{S,N_A}(N, Γ, J, periodic, ρ_A) where {S,N_A}
        divisible_by_half(S)
        new{S,N_A}(N, Γ, J, periodic, ρ_A)
    end
end

function TFIM(
    N::Int, N_A::Int, Γ::Real; S::Union{Int,Rational}=1//2, periodic::Bool=false, J::Real=-1
)
    ρ_A = rho_A(H_TFIM(N, Γ; periodic=periodic, J=J, S=S), N_A)

    return TFIM{Rational(S),N_A}(N, Γ, J, periodic, ρ_A)
end

"""
    H_TFIM(N::Int, Γ::Real; J::Real = -1, periodic::Bool=false, S::Union{Int64, Rational} = 1//2)

Return the TFIM Hamiltonian ``H= J(\\sum_{i=1}^{N-1} Z_{i}Z_{i+1} + Γ \\sum_{i=1}^{N} X_i)`` with `N` sites, transverse field `Γ` and 
prefactor `J` as a [`Block`](@ref).

Set `periodic` as true for PBC or as false for OBC and 
`S` for the spin number.
"""
function H_TFIM(
    N::Int, Γ::Real; J::Real=-1, periodic::Bool=false, S::Union{Int64,Rational}=1//2
)
    ising_term = sum(
        map(1:(periodic ? N : N - 1)) do i
            z(N, (i, i%N+1); S=S)
        end,
    )

    transverse_term = sum(
        map(1:N) do i
            x(N, i; S=S)
        end,
    )

    return Float64(J)*(ising_term + Float64(Γ)*transverse_term)
end

function hi(model::TFIM{S,N_A}, i::Int) where {S,N_A}
    hi = model.Γ * x(N_A, i; S=S)

    if i > 1
        hi += 1/2 * z(N_A, (i-1, i); S=S)
    end
    if i < N_A
        hi += 1/2 * z(N_A, (i, i+1); S=S)
    end
    return hi
end

function correction!(
    blocks::Vector{Block{S,N_A}}, ::TFIM{S,N_A}, i::Int, r::Int
) where {S,N_A}
    push!(blocks, z(N_A, (i, i+r); S=S))
end

function H_A_BWV_wo_corrections!(
    blocks::Vector{Block{S,N_A}}, model::TFIM{S,N_A}
) where {S,N_A}
    !iszero(model.Γ) && push!(blocks, model.Γ*x(N_A, 1; S=S))

    for i in 1:(N_A - 1)
        push!(blocks, z(N_A, (i, i+1); S=S))
        !iszero(model.Γ) && push!(blocks, model.Γ*x(N_A, i+1; S=S))
    end
end
