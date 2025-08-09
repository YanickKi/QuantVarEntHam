"""
    XXZ{S,N_A} <: AbstractModel{S,N_A}
    XXZ(N::Int, N_A::Int, Δ::Real; S::Union{Int, Rational} = 1//2, periodic::Bool = false, J::Real=+1, ϵ_reg::Real=1e-8)

Object containing the settings for the XXZ model

# Fields 
- see [`AbstractModel`](@ref)
- `Δ::Float64`: anisotropy 
"""
struct XXZ{S,N_A} <: AbstractModel{S,N_A}
    N::Int
    Δ::Float64
    J::Float64
    periodic::Bool
    ρ_A::Matrix{ComplexF64}
    function XXZ{S,N_A}(N, Δ, J, periodic, ρ_A) where {S,N_A}
        divisible_by_half(S)
        new{S,N_A}(N, Δ, J, periodic, ρ_A)
    end
end

function XXZ(
    N::Int, N_A::Int, Δ::Real; S::Union{Int,Rational}=1//2, periodic::Bool=false, J::Real=+1, ϵ_reg::Real=1e-8
)
    ρ_A = rho_A(H_XXZ(N, Δ; periodic=periodic, J=J, S=S), N_A, ϵ_reg)
    return XXZ{Rational(S),N_A}(N, Δ, J, periodic, ρ_A)
end

"""
    H_XXZ(N::Int, Δ::Real; periodic::Bool=false, J::Real = +1, S::Union{Int64, Rational} = 1//2)

Return the XXZ Hamiltonian ``H=J(\\sum_{i=1}^{N-1}( X_{i}X_{i+1} + Y_{i}Y_{i+1} + Δ Z_{i}Z_{i+1}))`` with `N` sites, anisotropy `Δ` and global prefactor `J` 
as a [`Block`](@ref).

Set `periodic` as true for PBC or as false for OBC and 
`S` for the spin number.
"""
function H_XXZ(
    N::Int, Δ::Real; periodic::Bool=false, J::Real=+1, S::Union{Int64,Rational}=1//2
)
    XX_term = sum(
        map(1:(periodic ? N : N - 1)) do i
            x(N, (i, i%N+1); S=S) + y(N, (i, i%N+1); S=S)
        end,
    )

    Z_term = sum(
        map(1:(periodic ? N : N - 1)) do i
            z(N, (i, i%N+1); S=S)
        end,
    )
    return Float64(J)*(XX_term+Float64(Δ)*Z_term)
end

function hi(model::XXZ{S,N_A}, i::Int) where {S,N_A}
    if i > 1 && i < N_A
        return 1/2*(
            x(N_A, (i-1, i); S=S) +
            y(N_A, (i-1, i); S=S) +
            model.Δ*z(N_A, (i-1, i); S=S)
        ) +
               1/2*(
            x(N_A, (i, i+1); S=S) +
            y(N_A, (i, i+1); S=S) +
            model.Δ*z(N_A, (i, i+1); S=S)
        )
    elseif i > 1
        return 1/2*(
            x(N_A, (i-1, i); S=S) +
            y(N_A, (i-1, i); S=S) +
            model.Δ*z(N_A, (i-1, i); S=S)
        )
    elseif i < N_A
        return 1/2*(
            x(N_A, (i, i+1); S=S) +
            y(N_A, (i, i+1); S=S) +
            model.Δ*z(N_A, (i, i+1); S=S)
        )
    end
end

function correction!(
    blocks::Vector{<:Block{S,N_A}}, model::XXZ{S,N_A}, i::Int, r::Int
) where {S,N_A}
    push!(
        blocks, x(N_A, (i, i+r); S=S) + y(N_A, (i, i+r); S=S)
    )
    !iszero(model.Δ) && push!(blocks, z(N_A, (i, i+r); S=S))
end

function H_A_BWV_wo_corrections!(
    blocks::Vector{<:Block{S,N_A}}, model::XXZ{S,N_A}
) where {S,N_A}
    for i in 1:(N_A - 1)
        push!(
            blocks,
            x(N_A, (i, i+1); S=S) + y(N_A, (i, i+1); S=S),
        )
        !iszero(model.Δ) && push!(blocks, model.Δ*z(N_A, (i, i+1); S=S))
    end
end
