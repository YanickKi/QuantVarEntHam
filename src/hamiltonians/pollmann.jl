include("qudits.jl")


@with_kw mutable struct Settings_Pollmann{T<:AbstractOperator, M<:AbstractMatrix} <: Settings_qudits{T, M}
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
    observables::Vector{T}
    meas0::Vector{Float64}
    mtrxObs::Vector{M}
end

function pollmann(N::Int, N_A::Int, J::Real, Bx::Real, Uzz::Real,T_max::Real; S::Rational=1//1, r_max::Int=1, periodic::Bool=false,
    signHam::Integer=+1, ρ_A::Operator=get_rhoA(H_pollmann(N, J , Bx, Uzz, periodic = periodic, signHam=signHam),  N-N_A+1:N, N),
    observables::Vector{<:AbstractOperator}=[repeat_cust(N_A, sigmaz, [i,i+1]) for i in 1:N_A-1])
    
    mtrxObs = [Diagonal(diag(obs.data)) for obs in observables]

    meas0 = [QuantumInterface.expect(obs, ρ_A) for obs in observables]

    return Settings_Pollmann{eltype(observables), eltype(mtrxObs)}(
        N = N, N_A = N_A, J = J, Bx = Bx, Uzz = Uzz, T_max = T_max, S=S, r_max = r_max, periodic = periodic,
        ρ_A = ρ_A.data, observables = observables,
        mtrxObs = mtrxObs, meas0 = meas0,
        signHam = signHam
    ) 
end 


function H_pollmann(N::Int, J::Real, Bx::Real, Uzz::Real; periodic::Bool=false, signHam::Integer = 1)

    heisenberg_term = map(1:(periodic ? N : N-1)) do i
        repeat_cust(N,sigmax,[i,i%N+1]) + repeat_cust(N,sigmay,[i,i%N+1]) +repeat_cust(N,sigmaz,[i,i%N+1])
    end |> sum

    transversal_term = map(1:N) do i
        repeat_cust(N, sigmax, [i])
    end |> sum

    Interaction_term = map(1:N) do i
        repeat_cust(N, sigmaz, [i])^2
    end |> sum
    return signHam*(J*heisenberg_term + Bx*transversal_term + Uzz*Interaction_term) 
end 


function hi(i::Int, set::Settings_Pollmann)
    @unpack N_A, J, Bx, Uzz = set
     
    hi = Bx*repeat_cust(N_A, sigmax, [i]) + Uzz*repeat_cust(N_A, sigmaz, [i])^2
    
    Ops = [sigmax, sigmay, sigmaz]

    if i > 1
        for Op in Ops
            hi += J/2  * repeat_cust(N_A, Op,[i,i-1])          #put(N_A, i-1=>Z) * put(N_A, i=>Z)
        end
    end
    if i < N_A
        for Op in Ops
            hi += J/2  * repeat_cust(N_A, Op,[i,i+1])          #put(N_A, i-1=>Z) * put(N_A, i=>Z)
        end
    end    
    return hi
end

function H_A_notBW_wo_corrections!(blks::Vector{<:AbstractOperator}, set::Settings_Pollmann)
    @unpack N_A, J, Bx, Uzz = set
    
    if iszero(Bx) == false  
        push!(blks, Bx*repeat_cust(N_A, sigmax, [1]))
    end 
    if iszero(Uzz) == false 
        push!(blks, Uzz*repeat_cust(N_A, sigmaz, [1])^2)
    end 

    for i ∈ 1:N_A-1 
        push!(blks, J*(repeat_cust(N_A,sigmaz,[i,i+1])))
        push!(blks, J*(repeat_cust(N_A,sigmax,[i,i+1])))
        push!(blks, J*(repeat_cust(N_A,sigmay,[i,i+1])))
        if iszero(Bx) == false 
            push!(blks, Bx*repeat_cust(N_A, sigmax, [i+1]))
        end 
        if iszero(Uzz) == false 
            push!(blks, Uzz*repeat_cust(N_A, sigmaz, [i+1])^2)
        end
    end 
    
end

function H_A_notBW_wo_corrections_I!(blks::Vector{<:AbstractOperator}, set::Settings_Pollmann)
    @unpack N_A, J, Bx, Uzz = set
    
    
    if iszero(Bx) == false  
        push!(blks, Bx*repeat_cust(N_A, sigmax, [1]))
    end 
    if iszero(Uzz) == false 
        push!(blks, Uzz*repeat_cust(N_A, sigmaz, [1])^2)
    end 

    for i ∈ 1:N_A-1 
        push!(blks, J*(repeat_cust(N_A,sigmaz,[i,i+1]) +repeat_cust(N_A,sigmax,[i,i+1]) +repeat_cust(N_A,sigmay,[i,i+1])))
        
        if iszero(Bx) == false 
            push!(blks, Bx*repeat_cust(N_A, sigmax, [i+1]))
        end 
        if iszero(Uzz) == false 
            push!(blks, Uzz*repeat_cust(N_A, sigmaz, [i+1])^2)
        end
    end 
    
end



#=

function correction!(blks::Vector{<:AbstractBlock}, i::Int, r::Int, set::Settings_TFIM)
    @unpack N_A = set   
    push!(blks, repeat_cust(N_A,Z,(i,i+r)))
end

=#