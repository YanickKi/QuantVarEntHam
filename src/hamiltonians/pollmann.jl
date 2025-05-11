@with_kw mutable struct Settings_Pollmann{M<:AbstractMatrix} <: Settings{M}
    N::Int
    N_A::Int
    J::Float64
    Bx::Float64
    Uzz::Float64
    T_max::Float64
    S::Rational
    r_max::Int
    periodic::Bool
    signHam::Int
    ρ_A::Matrix{ComplexF64}
    meas0::Vector{Float64}
    mtrxObs::Vector{M}
end

function pollmann(N::Int, N_A::Int, J::Real, Bx::Real, Uzz::Real,T_max::Real; S::Rational=1//1, r_max::Int=1, periodic::Bool=false,
    signHam::Integer=+1, ρ_A::Matrix{ComplexF64}=get_rhoA(H_pollmann(N, J , Bx, Uzz, periodic = periodic, signHam=signHam),  N-N_A+1:N, N, S=1),
    observables::Vector{<:AbstractMatrix}=[repeat(N_A, Z, (i,i+1), S=1) for i in 1:N_A-1])
    
    mtrxObs = [Diagonal(diag(obs)) for obs in observables]

    meas0 = [expect(obs, ρ_A) for obs in observables]

    return Settings_Pollmann{eltype(mtrxObs)}(
        N = N, N_A = N_A, J = J, Bx = Bx, Uzz = Uzz, T_max = T_max, S=S, r_max = r_max, periodic = periodic,
        ρ_A = ρ_A,
        mtrxObs = mtrxObs, meas0 = meas0,
        signHam = signHam
    ) 
end 


function H_pollmann(N::Int, J::Real, Bx::Real, Uzz::Real; periodic::Bool=false, signHam::Integer = 1)

    heisenberg_term = map(1:(periodic ? N : N-1)) do i
        repeat(N,X,(i,i%N+1), S=1) + repeat(N,Y,(i,i%N+1),S=1) +repeat(N,Z,(i,i%N+1), S=1)
    end |> sum

    transversal_term = map(1:N) do i
        repeat(N, X, i, S=1)
    end |> sum

    Interaction_term = map(1:N) do i
        repeat(N, Z, i, S=1)^2
    end |> sum
    return signHam*(J*heisenberg_term + Bx*transversal_term + Uzz*Interaction_term)         
end 


function hi(i::Int, set::Settings_Pollmann)
    @unpack N_A, J, Bx, Uzz = set
     
    hi = Bx*repeat(N_A, X, i, S=1) + Uzz*repeat(N_A, Z, i, S=1)^2
    
    Ops = [X, Y, Z]

    if i > 1
        for Op in Ops
            hi += J/2  * repeat(N_A, Op, (i,i-1), S=1)
        end
    end
    if i < N_A
        for Op in Ops
            hi += J/2  * repeat(N_A, Op, (i,i+1), S=1)
        end
    end    
    return hi
end

function H_A_notBW_wo_corrections!(blks::Vector{<:AbstractMatrix}, set::Settings_Pollmann)
    @unpack N_A, J, Bx, Uzz = set
    
    if iszero(Bx) == false  
        push!(blks, Bx*repeat(N_A, X, 1, S=1))
    end 
    if iszero(Uzz) == false 
        push!(blks, Uzz*repeat(N_A, Z, 1, S=1)^2)
    end 

    for i ∈ 1:N_A-1 
        push!(blks, J*(repeat(N_A,Z,(i,i+1), S=1)))
        push!(blks, J*(repeat(N_A,X,(i,i+1), S=1)))
        push!(blks, J*(repeat(N_A,Y,(i,i+1), S=1)))
        if iszero(Bx) == false 
            push!(blks, Bx*repeat(N_A, X, (i+1), S=1))
        end 
        if iszero(Uzz) == false 
            push!(blks, Uzz*repeat(N_A, Z, (i+1), S=1)^2)
        end
    end 
    
end

function H_A_notBW_wo_corrections_I!(blks::Vector{<:AbstractMatrix}, set::Settings_Pollmann)
    @unpack N_A, J, Bx, Uzz = set
    
    
    if iszero(Bx) == false  
        push!(blks, Bx*repeat(N_A, X, 1, S=1))
    end 
    if iszero(Uzz) == false 
        push!(blks, Uzz*repeat(N_A, Z, 1, S=1)^2)
    end 

    for i ∈ 1:N_A-1 
        push!(blks, J*(repeat(N_A,Z,(i,i+1)) +repeat(N_A,X,(i,i+1), S=1) +repeat(N_A,Y,(i,i+1), S=1)))
        
        if iszero(Bx) == false 
            push!(blks, Bx*repeat(N_A, X, (i+1), S=1))
        end 
        if iszero(Uzz) == false 
            push!(blks, Uzz*repeat(N_A, Z, (i+1), S=1)^2)
        end
    end 
    
end

#=

function correction!(blks::Vector{<:AbstractBlock}, i::Int, r::Int, set::Settings_TFIM)
    @unpack N_A = set   
    push!(blks, repeat(N_A,Z,(i,i+r)))
end

=#