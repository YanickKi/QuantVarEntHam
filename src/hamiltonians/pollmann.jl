"""
    Settings_pollmann

Concrete type of [`Settings`](@ref), containing the settings for the Pollmann model

!!! tip
    Use the constructor [`pollmann`](@ref) to instantiate this struct since the type for `observables` is automatically inferred
    and the default values in [`pollmann`](@ref) are highly recommended.

# Fields 
- see [`Settings`](@ref)
- `J_Heis::Float64`: Heisenberg coupling strength.
- `Bx::Float64`: transverse field strength
- `Uzz::Float64`: square term prefactor
"""
@with_kw struct Settings_pollmann <: Settings
    N::Int
    N_A::Int
    J_Heis::Float64
    Bx::Float64
    Uzz::Float64
    J::Float64
    S::Rational
    r_max::Int
    periodic::Bool
    ρ_A::Matrix{ComplexF64}
end

"""
    pollmann(N::Int, N_A::Int, J_Heis::Real, Bx::Real, Uzz::Real; S::Union{Int64, Rational}=1//1, r_max::Int=1, periodic::Bool=false,
    J::Real=+1, ρ_A::Matrix{ComplexF64}=get_rhoA(H_pollmann(N, J_Heis , Bx, Uzz, periodic = periodic, J=J, S = S),  N-N_A+1:N, N, S=S))

Convenient constructor for [`Settings_pollmann`](@ref) containing settings for the pollmann model

# Required Arguments
- `N::Int`: number of sites in the composite system.
- `N_A::Int`: number of sites in subsystem A.
- `J_Heis::Real`: Heisenberg coupling strength.
- `Bx::Real`: transverse field strength
- `Uzz::Real`: square term prefactor

# Keyword arguments
- `S::Union{Int64, Rational} = 1`: spin number.
- `J::Real=+1`: global prefactor in the Hamiltonian.
- `r_max::Int=1`: maximum range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) r_max = N_A-1 is maximally possible.
- `periodic::Bool=false`: boundary conditions for the system Hamiltonian, false for open and true for periodic boundary conditions.
- `ρ_A::AbstractMatrix=get_rhoA(H_TFIM(N, Γ, periodic = periodic, J=J),  N-N_A+1:N, N)`: reduced density matrix of ground state of the composite system on subsystem A, by default the subsystem is on the right border.

# Recommendations
- use only Z-gates (or composition of these) as observables since these are diagonal thus save computation time the most. 
"""
function pollmann(N::Int, N_A::Int, J_Heis::Real, Bx::Real, Uzz::Real; S::Union{Int64, Rational}=1, r_max::Int=1, periodic::Bool=false,
    J::Real=+1, ρ_A::Matrix{ComplexF64}=get_rhoA(H_pollmann(N, J_Heis , Bx, Uzz, periodic = periodic, J=J, S = S),  N-N_A+1:N, N, S=S))
    
    return Settings_pollmann(
        N = N, N_A = N_A,
        S=S, 
        J_Heis = J_Heis, Bx = Bx, Uzz = Uzz, J = J,
        r_max = r_max, periodic = periodic,
        ρ_A = ρ_A
    ) 
end 

"""
    H_pollmann(N::Int, J_Heis::Real, Bx::Real, Uzz::Real; periodic::Bool=false, J::Real = 1, S::Union{Int64, Rational}=1)

Return the pollmann Hamiltonian ``H= J(J_\text{Heis} \\sum_{i=1}^{N-1} \\vec{S}_i \\cdot \\vec{S}_{i+1} + B_x \\sum_{i=1}^{N} X_i + U_{zz} \\sum_{i=1}^{N} (Z_i)^2 )`` with `N` sites, Heisenberg coupling `J_Heis`,  
, transverse field strength `B_x`, quare term prefactor `U_{zz}` and global prefactor `J` as a sparse matrix.

Set `periodic` as true for PBC or as false for OBC and 
`S` for the spin number.

"""
function H_pollmann(N::Int, J_Heis::Real, Bx::Real, Uzz::Real; periodic::Bool=false, J::Real = 1, S::Union{Int64, Rational}=1)

    heisenberg_term = map(1:(periodic ? N : N-1)) do i
        repeat(N,X,(i,i%N+1), S=S) + repeat(N,Y,(i,i%N+1),S=S) +repeat(N,Z,(i,i%N+1), S=S)
    end |> sum

    transversal_term = map(1:N) do i
        repeat(N, X, i, S=S)
    end |> sum

    square_term = map(1:N) do i
        repeat(N, Z, i, S=S)^2
    end |> sum
    return J*(J_Heis*heisenberg_term + Bx*transversal_term + Uzz*square_term)         
end 


function hi(i::Int, set::Settings_pollmann)
    @unpack N_A, J_Heis, Bx, Uzz, S = set
     
    hi = Bx*repeat(N_A, X, i, S=S) + Uzz*repeat(N_A, Z, i, S=S)^2
    
    Ops = [X, Y, Z]

    if i > 1
        for Op in Ops
            hi += J_Heis/2  * repeat(N_A, Op, (i,i-1), S=S)
        end
    end
    if i < N_A
        for Op in Ops
            hi += J_Heis/2  * repeat(N_A, Op, (i,i+1), S=S)
        end
    end    
    return hi
end

function H_A_notBW_wo_corrections!(blks::Vector{<:AbstractMatrix}, set::Settings_pollmann)
    @unpack N_A, J_Heis, Bx, Uzz, S = set
    
    if iszero(Bx) == false  
        push!(blks, Bx*repeat(N_A, X, 1, S=S))
    end 
    if iszero(Uzz) == false 
        push!(blks, Uzz*repeat(N_A, Z, 1, S=S)^2)
    end 

    for i ∈ 1:N_A-1 
        push!(blks, J_Heis*(repeat(N_A,Z,(i,i+1), S=S)))
        push!(blks, J_Heis*(repeat(N_A,X,(i,i+1), S=S)))
        push!(blks, J_Heis*(repeat(N_A,Y,(i,i+1), S=S)))
        if iszero(Bx) == false 
            push!(blks, Bx*repeat(N_A, X, (i+1), S=S))
        end 
        if iszero(Uzz) == false 
            push!(blks, Uzz*repeat(N_A, Z, (i+1), S=S)^2)
        end
    end 
    
end

function H_A_notBW_wo_corrections_I!(blks::Vector{<:AbstractMatrix}, set::Settings_pollmann)
    @unpack N_A, J_Heis, Bx, Uzz, S = set
    
    
    if iszero(Bx) == false  
        push!(blks, Bx*repeat(N_A, X, 1, S=S))
    end 
    if iszero(Uzz) == false 
        push!(blks, Uzz*repeat(N_A, Z, 1, S=S)^2)
    end 

    for i ∈ 1:N_A-1 
        push!(blks, J_Heis*(repeat(N_A,Z,(i,i+1)) +repeat(N_A,X,(i,i+1), S=S) +repeat(N_A,Y,(i,i+1), S=S)))
        
        if iszero(Bx) == false 
            push!(blks, Bx*repeat(N_A, X, (i+1), S=S))
        end 
        if iszero(Uzz) == false 
            push!(blks, Uzz*repeat(N_A, Z, (i+1), S=S)^2)
        end
    end 
    
end

#=

function correction!(blks::Vector{<:AbstractBlock}, i::Int, r::Int, set::Settings_TFIM)
    @unpack N_A = set   
    push!(blks, repeat(N_A,Z,(i,i+r)))
end

=#