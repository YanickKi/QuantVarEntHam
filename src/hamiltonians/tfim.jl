@with_kw mutable struct Settings_TFIM{T<:AbstractBlock, S<:AbstractMatrix} <:Settings{T, S}
    N::Int
    N_A::Int
    Γ::Float64
    T_max::Float64
    r_max::Int
    periodic::Bool
    ρ_A::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}} 
    observables::Vector{T}
    meas0::Vector{Float64}  = [expect(observables[i], ρ_A) for i in 1:lastindex(observables)]
    mtrxObs::Vector{S}
    atol::Float64
    rtol::Float64 
end

function TFIM(N::Integer, N_A::Integer, Γ::Real, T_max::Real; r_max::Integer = 1, periodic::Bool = false, atol::Real=0.0, rtol::Real=atol>0 ? 0. : sqrt(eps(Float64)),
    signHam::Int = -1, ρ_A::DensityMatrix{2} = get_rhoA(H_TFIM(N, Γ, periodic = periodic, signHam = signHam),  N-N_A+1:N, N),
    observables::Vector{<:AbstractBlock} = [repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1])

    mtrxObs = mat.(observables)

    return Settings_TFIM{eltype(observables), eltype(mtrxObs)}(
        N = N, N_A = N_A, Γ = Γ, T_max = T_max, r_max = r_max,  periodic = periodic,
        atol = atol, rtol = rtol,
        ρ_A = ρ_A, observables = observables,
        mtrxObs = mtrxObs
    ) 
end 

function H_TFIM(N::Int64, Γ::Real; periodic::Bool=false, signHam::Integer = -1)

    ising_term = map(1:(periodic ? N : N-1)) do i
        repeat(N,Z,(i,i%N+1))
    end |> sum

    transversal_term = map(1:N) do i
        put(N, i =>X)
    end |> sum

    return signHam*(ising_term + Γ*transversal_term) 
end 


function hi(i::Integer, set::Settings_TFIM)
    @unpack N_A, Γ = set
     
    hi = -Γ * put(N_A, i=>X)
    
    if i > 1
        hi -= 1/2 * put(N_A, i-1=>Z) * put(N_A, i=>Z)
    end
    if i < N_A
        hi -= 1/2 * put(N_A, i=>Z) * put(N_A, i+1=>Z)
    end    
    return hi
end

function correction!(blks::Vector{AbstractBlock}, i::Integer, r::Integer, set::Settings_TFIM)
    @unpack N_A = set   
    push!(blks, repeat(N_A,Z,(i,i+r)))
end

function H_A_notBW_wo_corrections!(blks::Vector{AbstractBlock}, set::Settings_TFIM)
    @unpack N_A, Γ = set
    
    push!(blks, -Γ*put(N_A, 1 => X))

    for i ∈ 1:N_A-1 
        push!(blks, -repeat(N_A,Z,(i,i+1)))
        push!(blks, -Γ*put(N_A, i+1 => X))
    end 
    
end
