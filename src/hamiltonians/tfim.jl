@with_kw struct Settings_TFIM <:Settings
    N::Int
    N_A::Int
    Γ::Real
    T_max::Float64
    r_max::Int = 1
    periodic::Bool = false
    A::UnitRange{Int} = N-N_A+1:N      
    rhoA::DensityMatrix = get_rhoA(H_TFIM(N, Γ, periodic), A)
    observables::Vector{RepeatedBlock{2, 2, ZGate}}  = [repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1]
    meas0::Vector{Float64}  = [expect(observables[i], rhoA) for i in 1:lastindex(observables)]
    mtrxObs::Vector{Matrix{ComplexF64}} = Matrix.(observables)
    q::QuadTS = integration_tables()
end


function H_TFIM(N::Int64, Γ::Real, periodic::Bool)

    ising_term = map(1:(periodic ? N : N-1)) do i
        repeat(N,Z,(i,i%N+1))
    end |> sum

    transversal_term = map(1:N) do i
        put(N, i =>X)
    end |> sum
    return -(ising_term + Γ*transversal_term) 
end 


function hi(i::Integer, set::Settings_TFIM)
    @unpack N_A, A, Γ = set
     
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
    
    push!(blks, -Γ*put(N_A, 1 => Z))

    for i ∈ 1:N_A-1 
        push!(blks, -repeat(N_A,Z,(i,i+1)))
        push!(blks, -Γ*put(N_A, i+1 => X))
    end 
    
end
