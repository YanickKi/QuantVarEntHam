"""
    Pollmann{S,N_A} <: AbstractModel{S,N_A}
    Pollmann(N::Int, N_A::Int, J_Heis::Real, Bx::Real, Uzz::Real; S::Union{Int, Rational}=1, periodic::Bool=false, J::Real=+1, ϵ_reg::Real=1e-8)
    
Object containing the settings for the Pollmann model

# Fields 
- see [`AbstractModel`](@ref)
- `J_Heis::Float64`: Heisenberg coupling strength.
- `Bx::Float64`: transverse field strength
- `Uzz::Float64`: square term prefactor
"""
struct Pollmann{S,N_A} <: AbstractModel{S,N_A}
    N::Int
    J_Heis::Float64
    Bx::Float64
    Uzz::Float64
    J::Float64
    periodic::Bool
    ρ_A::Matrix{ComplexF64}
    function Pollmann{S,N_A}(N, J_Heis, Bx, Uzz, J, periodic, ρ_A) where {S,N_A}
        divisible_by_half(S)
        new{S,N_A}(N, J_Heis, Bx, Uzz, J, periodic, ρ_A)
    end
end

function Pollmann(
    N::Int,
    N_A::Int,
    J_Heis::Real,
    Bx::Real,
    Uzz::Real;
    S::Union{Int,Rational}=1,
    periodic::Bool=false,
    J::Real=+1,
    ϵ_reg::Real=1e-8
)
    ρ_A = rho_A(H_pollmann(N, J_Heis, Bx, Uzz; periodic=periodic, J=J, S=S), N_A, ϵ_reg)

    return Pollmann{Rational(S),N_A}(N, J_Heis, Bx, Uzz, J, periodic, ρ_A)
end

"""
    H_pollmann(N::Int, J_Heis::Real, Bx::Real, Uzz::Real; periodic::Bool=false, J::Real = 1, S::Union{Int64, Rational}=1)

Return the pollmann Hamiltonian ``H= J(J_\\text{Heis} \\sum_{i=1}^{N-1} \\vec{S}_i \\cdot \\vec{S}_{i+1} + B_x \\sum_{i=1}^{N} X_i + U_{zz} \\sum_{i=1}^{N} (Z_i)^2 )`` with `N` sites, Heisenberg coupling `J_Heis`,  
, transverse field strength `B_x`, square term prefactor `Uzz` and global prefactor `J` as a [`Block`](@ref).

Set `periodic` as true for PBC or as false for OBC and 
`S` for the spin number.
"""
function H_pollmann(
    N::Int,
    J_Heis::Real,
    Bx::Real,
    Uzz::Real;
    periodic::Bool=false,
    J::Real=1,
    S::Union{Int64,Rational}=1,
)
    heisenberg_term = sum(
        map(1:(periodic ? N : N - 1)) do i
            x(N, (i, i%N+1); S=S) +
            y(N, (i, i%N+1); S=S) +
            z(N, (i, i%N+1); S=S)
        end,
    )

    transverse_term = sum(
        map(1:N) do i
            x(N, i; S=S)
        end,
    )

    square_term = sum(
        map(1:N) do i
            z(N, i; S=S)^2
        end,
    )
    return Float64(
        J
    )*(J_Heis*heisenberg_term + Float64(Bx)*transverse_term + Float64(Uzz)*square_term)
end

function hi(model::Pollmann{S,N_A}, i::Int) where {S,N_A}
    hi = model.Bx*x(N_A, i; S=S) + model.Uzz*z(N_A, i; S=S)^2

    sigs = [x, y, z]

    if i > 1
        for sig in sigs
            hi += model.J_Heis/2 * sig(N_A, (i, i-1); S=S)
        end
    end
    if i < N_A
        for sig in sigs
            hi += model.J_Heis/2 * sig(N_A, (i, i+1); S=S)
        end
    end
    return hi
end

function H_A_BWV_wo_corrections!(
    blocks::Vector{Block{S,N_A}}, model::Pollmann{S,N_A}
) where {S,N_A}
    !iszero(model.Bx) && push!(blocks, model.Bx*x(N_A, 1; S=S))
    !iszero(model.Uzz) && push!(blocks, model.Uzz*z(N_A, 1; S=S)^2)

    sigs = [x, y, z]

    for i in 1:(N_A - 1)
        for sig in sigs
            !iszero(model.J_Heis) &&
                push!(blocks, model.J_Heis*sig(N_A, (i, i+1); S=S))
        end
        !iszero(model.Bx) && push!(blocks, model.Bx*x(N_A, (i+1); S=S))
        !iszero(model.Uzz) && push!(blocks, model.Uzz*z(N_A, (i+1); S=S)^2)
    end
end

function correction!(
    blocks::Vector{Block{S,N_A}}, model::Pollmann{S,N_A}, i::Int, r::Int
) where {S,N_A}
    iszero(model.J_Heis) && return nothing
    for sig in [x, y, z]
        push!(blocks, sig(N_A, (i, i+r); S=S))
    end
end
