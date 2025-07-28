"""
    Pollmann{S,N_A} <: AbstractModel{S,N_A}
    Pollmann(N::Int, N_A::Int, J_Heis::Real, Bx::Real, Uzz::Real; S::Union{Int, Rational}=1, periodic::Bool=false, J::Real=+1)
    
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
    r_max::Int
    periodic::Bool
    ρ_A::Matrix{ComplexF64}
    function Pollmann{S,N_A}(N, J_Heis, Bx, Uzz, J, periodic, ρ_A) where {S,N_A}
        check_S_N_type(S, N_A)
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
)
    ρ_A = rho_A(H_pollmann(N, J_Heis, Bx, Uzz; periodic=periodic, J=J, S=S), N_A)

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
    heisenberg_term = sum(map(1:(periodic ? N : N - 1)) do i
        PauliString(N, "X", (i, i%N+1); S=S) +
        PauliString(N, "Y", (i, i%N+1); S=S) +
        PauliString(N, "Z", (i, i%N+1); S=S)
    end)

    transverse_term = sum(map(1:N) do i
        PauliString(N, "X", i; S=S)
    end)

    square_term = sum(map(1:N) do i
        PauliString(N, "Z", i; S=S)^2
    end)
    return Float64(
        J
    )*(J_Heis*heisenberg_term + Float64(Bx)*transverse_term + Float64(Uzz)*square_term)
end

function hi(model::Pollmann{S,N_A}, i::Int) where {S,N_A}
    hi = model.Bx*PauliString(N_A, "X", i; S=S) + model.Uzz*PauliString(N_A, "Z", i; S=S)^2

    sigs = ["X", "Y", "Z"]

    if i > 1
        for sig in sigs
            hi += model.J_Heis/2 * PauliString(N_A, sig, (i, i-1); S=S)
        end
    end
    if i < N_A
        for sig in sigs
            hi += model.J_Heis/2 * PauliString(N_A, sig, (i, i+1); S=S)
        end
    end
    return hi
end

function H_A_BWV_wo_corrections!(
    blocks::Vector{Block{S,N_A}}, model::Pollmann{S,N_A}
) where {S,N_A}
    !iszero(model.Bx) && push!(blocks, model.Bx*PauliString(N_A, "X", 1; S=S))
    !iszero(model.Uzz) && push!(blocks, model.Uzz*PauliString(N_A, "Z", 1; S=S)^2)

    sigs = ["X", "Y", "Z"]

    for i in 1:(N_A - 1)
        for sig in sigs
            !iszero(model.J_Heis) &&
                push!(blocks, model.J_Heis*PauliString(N_A, sig, (i, i+1); S=S))
        end
        !iszero(model.Bx) && push!(blocks, model.Bx*PauliString(N_A, "X", (i+1); S=S))
        !iszero(model.Uzz) && push!(blocks, model.Uzz*PauliString(N_A, "Z", (i+1); S=S)^2)
    end
end

function correction!(
    blocks::Vector{Block{S,N_A}}, model::Pollmann{S,N_A}, i::Int, r::Int
) where {S,N_A}
    iszero(model.J_Heis) && return nothing
    for sig in ["X", "Y", "Z"]
        push!(blocks, PauliString(N_A, sig, (i, i+r); S=S))
    end
end

#=
function H_A_notBW_wo_corrections_I!(blocks::Vector{<:AbstractMatrix}, model::Pollmann)
    @unpack N_A, J_Heis, Bx, Uzz, S = model

    if iszero(Bx) == false  
        push!(blocks, Bx*repeat(N_A, X, 1, S=S))
    end 
    if iszero(Uzz) == false 
        push!(blocks, Uzz*repeat(N_A, Z, 1, S=S)^2)
    end 

    for i ∈ 1:N_A-1 
        push!(blocks, J_Heis*(repeat(N_A,Z,(i,i+1)) +repeat(N_A,X,(i,i+1), S=S) +repeat(N_A,Y,(i,i+1), S=S)))

        if iszero(Bx) == false 
            push!(blocks, Bx*repeat(N_A, X, (i+1), S=S))
        end 
        if iszero(Uzz) == false 
            push!(blocks, Uzz*repeat(N_A, Z, (i+1), S=S)^2)
        end
    end 

end

=#
