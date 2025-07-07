"""
    Pollmann <: AbstractModel

Object containing the settings for the Pollmann model

# Fields 
- see [`AbstractModel`](@ref)
- `J_Heis::Float64`: Heisenberg coupling strength.
- `Bx::Float64`: transverse field strength
- `Uzz::Float64`: square term prefactor
"""
@with_kw struct Pollmann <: AbstractModel
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
    Pollmann(N::Int, N_A::Int, J_Heis::Real, Bx::Real, Uzz::Real; S::Union{Int64, Rational}=1//1, r_max::Int=1, periodic::Bool=false,
    J::Real=+1, ρ_A::Matrix{ComplexF64}=get_ρ_A(H_pollmann(N, J_Heis , Bx, Uzz, periodic = periodic, J=J, S = S),  N-N_A+1:N, N, S=S))

Convenient constructor for [`Pollmann`](@ref) containing settings for the Pollmann model.
The default values are often used and the density matrix is automatically constructed.

# Required Arguments
- `N`: number of sites in the composite system.
- `N_A`: number of sites in subsystem A.
- `J_Heis`: Heisenberg coupling strength.
- `Bx`: transverse field strength
- `Uzz`: square term prefactor

# Keyword arguments
- `S`: spin number.
- `J`: global prefactor in the Hamiltonian.
- `r_max`: maximum range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) `r_max = N_A-1` is maximally possible.
- `periodic`: boundary conditions for the system Hamiltonian, false for open and true for periodic boundary conditions.
- `ρ_A`: reduced density matrix of ground state of the composite system on subsystem A, by default the subsystem is on the right border.
"""
function Pollmann(N::Int, N_A::Int, J_Heis::Real, Bx::Real, Uzz::Real; S::Union{Int64, Rational}=1, r_max::Int=1, periodic::Bool=false,
    J::Real=+1, ρ_A::Matrix{ComplexF64}=get_rhoA(H_pollmann(N, J_Heis , Bx, Uzz, periodic = periodic, J=J, S = S),  N-N_A+1:N, N, S=S))
    
    return Pollmann(
        N = N, N_A = N_A,
        S=S, 
        J_Heis = J_Heis, Bx = Bx, Uzz = Uzz, J = J,
        r_max = r_max, periodic = periodic,
        ρ_A = ρ_A
    ) 
end 

"""
    H_pollmann(N::Int, J_Heis::Real, Bx::Real, Uzz::Real; periodic::Bool=false, J::Real = 1, S::Union{Int64, Rational}=1)

Return the pollmann Hamiltonian ``H= J(J_\\text{Heis} \\sum_{i=1}^{N-1} \\vec{S}_i \\cdot \\vec{S}_{i+1} + B_x \\sum_{i=1}^{N} X_i + U_{zz} \\sum_{i=1}^{N} (Z_i)^2 )`` with `N` sites, Heisenberg coupling `J_Heis`,  
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


function hi(model::Pollmann, i::Int)
    @unpack N_A, J_Heis, Bx, Uzz, S = model
     
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

function H_A_notBW_wo_corrections!(blks::Vector{<:AbstractMatrix}, model::Pollmann)
    @unpack N_A, J_Heis, Bx, Uzz, S = model
    
    if iszero(Bx) == false  
        push!(blks, Bx*repeat(N_A, X, 1, S=S))
    end 
    if iszero(Uzz) == false 
        push!(blks, Uzz*repeat(N_A, Z, 1, S=S)^2)
    end 

    for i ∈ 1:N_A-1 
        push!(blks, J_Heis*(repeat(N_A,X,(i,i+1), S=S)))
        push!(blks, J_Heis*(repeat(N_A,Y,(i,i+1), S=S)))
        push!(blks, J_Heis*(repeat(N_A,Z,(i,i+1), S=S)))
        if iszero(Bx) == false 
            push!(blks, Bx*repeat(N_A, X, (i+1), S=S))
        end 
        if iszero(Uzz) == false 
            push!(blks, Uzz*repeat(N_A, Z, (i+1), S=S)^2)
        end
    end 
    
end


function correction!(blks::Vector{<:AbstractMatrix}, model::Pollmann, i::Int, r::Int)
    @unpack N_A, S = model   
    for σ in [X,Y,Z]
        push!(blks, repeat(N_A,σ,(i,i+r), S=S))

    end
end

#=
function H_A_notBW_wo_corrections_I!(blks::Vector{<:AbstractMatrix}, model::Pollmann)
    @unpack N_A, J_Heis, Bx, Uzz, S = model
    
    
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

=#