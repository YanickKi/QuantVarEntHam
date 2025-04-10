"""
    Settings_TFIM{T<:AbstractBlock, S<:AbstractMatrix} <:Settings{T,S}

Contains the settings for the TFIM

The parametric type `T<:AbstractBlock` is introduced for determining the correct concrete type of the Yao Blocks, while 
`S<:AbstractMatrix` is needed to determine the concrete type of the matrix representation of the Yao Blocks to prevent 
working with full complex dense matrices.

!!! tip
    Use the constructor [`TFIM`](@ref) to instantiate this struct since the types for `observables` and its matrices `mtrxObs` are automatically inferred then
    and the default values in [`TFIM`](@ref) are highly recommended.

# Fields 
- `N::Int`: number of sites in composite system.
- `N_A::Int`: number of sites in subsystem A.
- `Γ::Float64`: Anisotropy.
- `T_max::Float64`: maximum time for evolving the observables i.e. maximum integration time.
- `r_max::Int`: range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) r_max = N_A-1 corresponds to maximum order.
- `periodic::Bool`: boundary conditions for the system Hamitlonian, false for open and true for periodic boundary conditions, obsolete if an own reduced density matrix ρ_A is provided.
- `signHam::Int`: global minus sign of Hamiltonian
- `ρ_A::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}}`: reduced density matrix of ground state of the composite system on subsystem A.
- `observables::Vector{T}`: monitored observables in the cost function.
- `meas0::Vector{Float64} = [expect(observables[i], ρ_A) for i in eachindex(observables)]`: expectation values of `observables` at time ``t=0``.
- `mtrxObs::Vector{S}`: matrix representations for the observables
- `dt::Float64=0.01`: time step for evaluating the cost function via midpoint rule, obsolete if other integration techniques are used. 
"""
@with_kw mutable struct Settings_TFIM{T<:AbstractBlock, S<:AbstractMatrix} <:Settings{T, S}
    N::Int
    N_A::Int
    Γ::Float64
    T_max::Float64
    r_max::Int
    periodic::Bool
    signHam::Int
    ρ_A::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}} 
    observables::Vector{T}
    meas0::Vector{Float64}  = [expect(observables[i], ρ_A) for i in eachindex(observables)]
    mtrxObs::Vector{S}
    dt::Float64
end

"""
    TFIM(N::Int, N_A::Int, Γ::Real, T_max::Real; r_max::Int=1, periodic::Bool=false, atol::Real=0.0, rtol::Real=atol>0 ? 0. : sqrt(eps(Float64)),
    signHam::Integer=-1, ρ_A::DensityMatrix{2}=get_rhoA(H_TFIM(N, Γ, periodic = periodic, signHam=signHam), N-N_A+1:N, N),
    observables::Vector{<:AbstractBlock}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1], dt::Real=0.01) 

Convenient constructor for [`Settings_TFIM`](@ref) containing settings for the Transversal Field Ising Model

# Required Arguments
- `N::Int`: number of sites in the composite system.
- `Γ::Real`: transversal field strength 
- `N_A::Int`: number of sites in subsystem A.
- `T_max::Real`: maximum time for evolving the observables i.e. maximum integration time.

# Keyword arguments
- `r_max::Int=1`: range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) r_max = N_A-1 corresponds to maximum order.
- `periodic::Bool=false`: boundary conditions for the system Hamitlonian, false for open and true for periodic boundary conditions.
- `atol::Real=0.0`: absolute tolerance for the integrator.
- `rtol::Real=atol>0 ? 0.: sqrt(eps(Float64))`: relative tolerance for the integrator.
- `signHam::Integer=-1`: global sign of Hamiltonian, obsolete if an own reduced density matrix ρ_A is provided.
- `ρ_A::DensityMatrix{2}=get_rhoA(H_TFIM(N, Γ, periodic = periodic, signHam=signHam), N-N_A+1:N, N)`: reduced density matrix of ground state of the composite system on subsystem A, by default the subsystem is on the right border.
- `observables::Vector{<:AbstractBlock}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1]`: monitored observables in the cost function.
- `dt::Real=0.01`: time step for evaluating the cost function via midpoint rule, obsolete if other integration techniques are used.

# Recommendations
- use only Z-gates (or composition of these) as observables since these are diagonal thus save computation time the most. 
- use only one type of observables (e.g. X_i X_i+1) since these are then stored as sparse matrices, otherwise dense which leads to higher computation time.
- be carefull when changing the tolerances for integration (not recommended), a relative tolerance higher than ≈ 1e-7 is not recommended since this can lead to wrong results.
"""
function TFIM(N::Int, N_A::Int, Γ::Real, T_max::Real; r_max::Int=1, periodic::Bool=false,
    signHam::Integer=-1, ρ_A::DensityMatrix{2}=get_rhoA(H_TFIM(N, Γ, periodic = periodic, signHam=signHam),  N-N_A+1:N, N),
    observables::Vector{<:AbstractBlock}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1], dt::Real=0.01)

    mtrxObs = mat.(observables)

    return Settings_TFIM{eltype(observables), eltype(mtrxObs)}(
        N = N, N_A = N_A, Γ = Γ, T_max = T_max, r_max = r_max,  periodic = periodic,
        ρ_A = ρ_A, observables = observables,
        mtrxObs = mtrxObs,
        dt = dt,
        signHam = signHam
    ) 
end 

"""
    H_TFIM(N::Int, Γ::Real; periodic::Bool=false, signHam::Integer = -1)

Return the TFIM Hamiltonian ``H=\\sum_{i=1}^{N-1} Z_{i}Z_{i+1} + Γ \\sum_{i=1}^{N} X_i`` with `N` sites and transversal field `Γ` as an AbstractBlock.

Set `periodic` as true for PBC or as false for OBC and `signHam` for a global sign.

"""
function H_TFIM(N::Int, Γ::Real; periodic::Bool=false, signHam::Integer = -1)

    ising_term = map(1:(periodic ? N : N-1)) do i
        repeat(N,Z,(i,i%N+1))
    end |> sum

    transversal_term = map(1:N) do i
        put(N, i =>X)
    end |> sum

    return signHam*(ising_term + Γ*transversal_term) 
end 


function hi(i::Int, set::Settings_TFIM)
    @unpack N_A, Γ = set
     
    hi = Γ * put(N_A, i=>X)
    
    if i > 1
        hi += 1/2 * put(N_A, i-1=>Z) * put(N_A, i=>Z)
    end
    if i < N_A
        hi += 1/2 * put(N_A, i=>Z) * put(N_A, i+1=>Z)
    end    
    return hi
end

function correction!(blks::Vector{<:AbstractBlock}, i::Int, r::Int, set::Settings_TFIM)
    @unpack N_A = set   
    push!(blks, repeat(N_A,Z,(i,i+r)))
end

function H_A_notBW_wo_corrections!(blks::Vector{<:AbstractBlock}, set::Settings_TFIM)
    @unpack N_A, Γ = set
    
    push!(blks, Γ*put(N_A, 1 => X))

    for i ∈ 1:N_A-1 
        push!(blks, repeat(N_A,Z,(i,i+1)))
        push!(blks, Γ*put(N_A, i+1 => X))
    end 
    
end
