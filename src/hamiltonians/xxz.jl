"""
    XXZ

Concrete type of [`AbstractModel`](@ref), containing the settings for the XXZ model

!!! tip
    Use the constructor [`XXZ`](@ref) to instantiate this struct since the type for `observables` is automatically inferred
    and the default values in [`XXZ`](@ref) are highly recommended.

# Fields 
- see [`AbstractModel`](@ref)
- `Δ::Real`: anisotropy 
"""
@with_kw struct XXZ <: AbstractModel
    N::Int
    N_A::Int
    Δ::Float64
    J::Float64
    S::Rational
    r_max::Int
    periodic::Bool 
    ρ_A::Matrix{ComplexF64}
end
"""
    XXZ(N::Int, N_A::Int, Δ::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool = false,
    J::Real=+1, ρ_A::Matrix{ComplexF64}=get_ρ_A(H_XXZ(N, Δ, periodic=periodic, J=J, S = S),  N-N_A+1:N, N))

Convenient constructor for [`XXZ`](@ref) containing settings for the XXZ Model 

# Required Arguments
- `N::Int`: number of sites in the composite system.
- `Δ::Real`: anisotropy 
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
function XXZ(N::Int, N_A::Int, Δ::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool = false,
    J::Real=+1, ρ_A::Matrix{ComplexF64}=get_ρ_A(H_XXZ(N, Δ, periodic=periodic, J=J, S = S),  N-N_A+1:N, N))
    
    return XXZ(
        N = N, N_A = N_A,
        S = S,
        Δ = Δ, J = J, 
        r_max = r_max, periodic = periodic,
        ρ_A = ρ_A
    ) 
end 

"""
    H_XXZ(N::Int, Δ::Real; periodic::Bool=false, J::Real = +1, S::Union{Int64, Rational} = 1//2)

Return the XXZ Hamiltonian ``H=J(\\sum_{i=1}^{N-1}( X_{i}X_{i+1} + Y_{i}Y_{i+1} + Δ Z_{i}Z_{i+1}))`` with `N` sites, anisotropy `Δ` and global prefactor `J` as a sparse matrix.

Set `periodic` as true for PBC or as false for OBC and 
`S` for the spin number.
"""
function H_XXZ(N::Int, Δ::Real; periodic::Bool=false, J::Real = +1, S::Union{Int64, Rational} = 1//2)
    
    XX_term = map(1:(periodic ? N : N-1)) do i
        repeat(N,X,(i,i%N+1), S=S) + repeat(N,Y,(i,i%N+1), S=S)
    end |> sum

    Z_term = map(1:(periodic ? N : N-1)) do i
        repeat(N,Z,(i,i%N+1), S=S)
    end |> sum
    return J*(XX_term+Δ*Z_term)
end 

function hi(model::XXZ, i::Int)
    @unpack N_A, Δ, S= model

    if i > 1 && i < N_A 
        return 1/2*(repeat(N_A,X,(i-1,i), S=S) + repeat(N_A,Y,(i-1,i), S=S)+ Δ*repeat(N_A,Z,(i-1,i), S=S)) + 1/2*(repeat(N_A,X,(i,i+1), S=S) + repeat(N_A,Y,(i,i+1), S=S) + Δ*repeat(N_A,Z,(i,i+1), S=S))
    elseif i > 1
        return 1/2*(repeat(N_A,X,(i-1,i), S=S) + repeat(N_A,Y,(i-1,i), S=S)+ Δ*repeat(N_A,Z,(i-1,i), S=S)) 
    elseif i < N_A
        return 1/2*(repeat(N_A,X,(i,i+1), S=S) + repeat(N_A,Y,(i,i+1), S=S) + Δ*repeat(N_A,Z,(i,i+1), S=S))
    end    
end


function correction!(blks::Vector{<:AbstractMatrix}, model::XXZ, i::Int, r::Int)
    @unpack N_A, Δ, S = model   
    push!(blks, repeat(N_A,X,(i,i+r), S=S) + repeat(N_A,Y,(i,i+r), S=S))
    push!(blks, Δ*repeat(N_A,Z,(i,i+r), S=S)) 
end

function H_A_notBW_wo_corrections!(blks::Vector{<:AbstractMatrix}, model::XXZ)
    @unpack N_A, Δ, S = model
     
    for i ∈ 1:N_A-1 
        push!(blks, repeat(N_A,X,(i,i+1), S=S) + repeat(N_A,Y,(i,i+1), S=S))
        push!(blks, Δ*repeat(N_A,Z,(i,i+1), S=S))
    end 
end

function H_A_XYZ_wo_corrections!(blks::Vector{<:AbstractMatrix}, model::XXZ)
    @unpack N_A, Δ, S= model
     
    for i ∈ 1:N_A-1 
        push!(blks, repeat(N_A,X,(i,i+1), S=S))
        push!(blks, repeat(N_A,Y,(i,i+1), S=S))
        push!(blks, Δ*repeat(N_A,Z,(i,i+1), S=S))
    end 
end
