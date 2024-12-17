"""
    Settings_XXZ{T<:AbstractBlock, S<:AbstractMatrix} <:Settings{T,S}

Contains the settings for the XXZ model

The parametric type `T<:AbstractBlock` is introduced for determining the correct concrete type of the Yao Blocks, while 
`S<:AbstractMatrix` is needed to determine the concrete type of the matrix representation of the Yao Blocks to prevent 
working with full complex dense matrices.

!!! tip
    Use the constructor [`XXZ`](@ref) to instantiate this struct since the types for `observables` and its matrices `mtrxObs` are automatically inferred then
    and the default values in [`XXZ`](@ref) are highly recommended.

# Fields 
- `N::Int`: number of sites in composite system.
- `N_A::Int`: number of sites in subsystem A.
- `Δ::Float64`: Anisotropy.
- `T_max::Float64`: maximum time for evolving the observables i.e. maximum integration time.
- `r_max::Int`: range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) r_max = N_A-1 corresponds to maximum order.
- `periodic::Bool`: boundary conditions for the system Hamitlonian, false for open and true for periodic boundary conditions, obsolete if an own reduced density matrix ρ_A is provided.
- `ρ_A::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}}`: reduced density matrix of ground state of the composite system on subsystem A.
- `observables::Vector{T}`: monitored observables in the cost function.
- `meas0::Vector{Float64} = [expect(observables[i], ρ_A) for i in eachindex(observables)]`: expectation values of `observables` at time ``t=0``.
- `mtrxObs::Vector{S}`: matrix representations for the observables
- `atol::Float64=0.0`: absolute tolerance for the integrator.
- `rtol::Float64=atol>0 ? 0.: sqrt(eps(Float64))`: relative tolerance for the integrator.
- `dt::Float64=0.01`: time step for evaluating the cost function via midpoint rule, obsolete if other integration techniques are used.

"""
@with_kw mutable struct Settings_XXZ{T<:AbstractBlock, S<:AbstractMatrix} <:Settings{T,S}
    N::Int
    N_A::Int
    Δ::Float64
    T_max::Float64
    r_max::Int
    periodic::Bool 
    ρ_A::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}} 
    observables::Vector{T}
    meas0::Vector{Float64} = [expect(observables[i], ρ_A) for i in eachindex(observables)]
    mtrxObs::Vector{S}
    atol::Float64
    rtol::Float64 
    dt::Float64
end

"""
    XXZ(N::Int, N_A::Int, Δ::Real, T_max::Real; r_max::Int=1, periodic::Bool = false, atol::Real=0.0, rtol::Real=atol>0 ? 0. : sqrt(eps(Float64)),
    signHam::Integer=+1, ρ_A::DensityMatrix{2}=get_rhoA(H_XXZ(N, Δ, periodic=periodic, signHam=signHam), N-N_A+1:N, N),
    observables::Vector{<:AbstractBlock}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1], dt::Float64 = 0.01)

Convenient constructor for [`Settings_XXZ`](@ref) containing settings for the XXZ Model 

# Required Arguments
- `N::Int`: number of sites in the composite system.
- `Δ::Real`: Anisotropy 
- `N_A::Int`: number of sites in subsystem A.
- `T_max::Real`: maximum time for evolving the observables i.e. maximum integration time.

# Keyword arguments
- `r_max::Int=1`: range of interaction (1 for nearest neighbour, 2 for next nearest neighbour, etc..) r_max = N_A-1 corresponds to maximum order.
- `periodic::Bool=false`: boundary conditions for the system Hamitlonian, false for open and true for periodic boundary conditions, obsolete if an own reduced density matrix ρ_A is provided.
- `atol::Real=0.0`: absolute tolerance for the integrator.
- `rtol::Real=atol>0 ? 0.: sqrt(eps(Float64))`: relative tolerance for the integrator.
- `signHam::Integer=1`: global sign of Hamiltonian, obsolete if an own reduced density matrix ρ_A is provided.
- `ρ_A::DensityMatrix{2}=get_rhoA(H_XXZ(N, Δ, periodic=periodic, signHam=signHam), N-N_A+1:N, N)`: reduced density matrix of ground state of the composite system on subsystem A, by default the subsystem is on the right border.
- `observables::Vector{<:AbstractBlock}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1]`: monitored observables in the cost function.
- `dt::Real=0.01`: time step for evaluating the cost function via midpoint rule, obsolete if other integration techniques are used.

# Recommendations
- use only Z-gates (or composition of these) as observables since these are diagonal thus save computation time the most.
- use only one type of observables (e.g. X_i X_i+1) since these are then stored as sparse or diagonal matrices, otherwise dense which leads to higher computation time.
- be carefull when changing the tolerances for integration (not recommended), a relative tolerance higher than ≈ 1e-7 is not recommended since this can lead to wrong results.
"""
function XXZ(N::Int, N_A::Int, Δ::Real, T_max::Real; r_max::Int=1, periodic::Bool = false, atol::Real=0.0, rtol::Real=atol>0 ? 0. : sqrt(eps(Float64)),
    signHam::Integer=+1, ρ_A::DensityMatrix{2}=get_rhoA(H_XXZ(N, Δ, periodic=periodic, signHam=signHam),  N-N_A+1:N, N),
    observables::Vector{<:AbstractBlock}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1], dt::Float64 = 0.01)
    
    mtrxObs = mat.(observables)

    return Settings_XXZ{eltype(observables), eltype(mtrxObs)}(
        N = N, N_A = N_A, Δ = Δ, T_max = T_max, r_max = r_max, periodic = periodic,
        atol = atol, rtol = rtol,
        ρ_A = ρ_A, observables = observables,
        mtrxObs = mtrxObs,
        dt = dt 
    ) 
end 

"""
    H_XXZ(N::Int,  Δ::Real; periodic::Bool=false, signHam::Integer = +1)

Return the XXZ Hamiltonian ``H=\\sum_{i=1}^{N-1}( X_{i}X_{i+1} + Y_{i}Y_{i+1} + Δ Z_{i}Z_{i+1})`` with `N` sites and Anisotropy `Δ` as an AbstractBlock.

Set `periodic` as true for PBC or as false for OBC and `signHam` for a global sign.
"""
function H_XXZ(N::Int, Δ::Real; periodic::Bool=false, signHam::Integer = +1)
    
    XX_term::Add{2} = map(1:(periodic ? N : N-1)) do i
        repeat(N,X,(i,i%N+1)) + repeat(N,Y,(i,i%N+1))
    end |> sum

    Z_term::Add{2} = map(1:(periodic ? N : N-1)) do i
        repeat(N,Z,(i,i%N+1))
    end |> sum
    return signHam*(XX_term+Δ*Z_term)
end 

function hi(i::Int, set::Settings_XXZ)
    @unpack N_A, Δ = set

    if i > 1 && i < N_A 
        return 1/2*(repeat(N_A,X,(i-1,i)) + repeat(N_A,Y,(i-1,i))+ Δ*repeat(N_A,Z,(i-1,i))) + 1/2*(repeat(N_A,X,(i,i+1)) + repeat(N_A,Y,(i,i+1)) + Δ*repeat(N_A,Z,(i,i+1)))
    elseif i > 1
        return 1/2*(repeat(N_A,X,(i-1,i)) + repeat(N_A,Y,(i-1,i))+ Δ*repeat(N_A,Z,(i-1,i))) 
    elseif i < N_A
        return 1/2*(repeat(N_A,X,(i,i+1)) + repeat(N_A,Y,(i,i+1)) + Δ*repeat(N_A,Z,(i,i+1)))
    end    
end


function correction!(blks::Vector{AbstractBlock}, i::Int, r::Int, set::Settings_XXZ)
    @unpack N_A = set   
    push!(blks, repeat(N_A,X,(i,i+r)) + repeat(N_A,Y,(i,i+r)))
    push!(blks, repeat(N_A,Z,(i,i+r))) 
end

function H_A_notBW_wo_corrections!(blks::Vector{AbstractBlock}, set::Settings_XXZ)
    @unpack N_A, Δ = set
     
    for i ∈ 1:N_A-1 
        push!(blks, repeat(N_A,X,(i,i+1)) + repeat(N_A,Y,(i,i+1)))
        push!(blks, Δ*repeat(N_A,Z,(i,i+1)))
    end 
end
