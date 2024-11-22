@with_kw struct Settings_XXZ <:Settings
    N::Int
    N_A::Int
    Δ::Real
    T_max::Float64
    r_max::Int = 1
    periodic::Bool = false
    A::UnitRange{Int} = N-N_A+1:N      
    signHam::Int = +1
    rhoA::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}} = get_rhoA(H_XXZ(N, Δ, periodic=periodic), A)
    observables::Vector{RepeatedBlock{2}}  = [repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1]
    meas0::Vector{Float64}  = [expect(observables[i], rhoA) for i in 1:lastindex(observables)]
    mtrxObs::Vector{Matrix{ComplexF64}} = Matrix.(observables)
    q::QuadTS = integration_tables()
end

function H_XXZ(N::Integer, Δ::Real; periodic::Bool=false, signHam::Integer = +1)
    
    XX_term::Add{2} = map(1:(periodic ? N : N-1)) do i
        repeat(N,X,(i,i%N+1)) + repeat(N,Y,(i,i%N+1))
    end |> sum

    Z_term::Add{2} = map(1:(periodic ? N : N-1)) do i
        repeat(N,Z,(i,i%N+1))
    end |> sum
    return signHam*(XX_term+Δ*Z_term)
end 

function hi(i::Integer, set::Settings_XXZ)::AbstractBlock
    @unpack N_A, A, Δ = set

    if i > 1 && i < N_A 
        return 1/2*(repeat(N_A,X,(i-1,i)) + repeat(N_A,Y,(i-1,i))+ Δ*repeat(N_A,Z,(i-1,i))) + 1/2*(repeat(N_A,X,(i,i+1)) + repeat(N_A,Y,(i,i+1)) + Δ*repeat(N_A,Z,(i,i+1)))
    elseif i > 1
        return 1/2*(repeat(N_A,X,(i-1,i)) + repeat(N_A,Y,(i-1,i))+ Δ*repeat(N_A,Z,(i-1,i))) 
    elseif i < N_A
        return 1/2*(repeat(N_A,X,(i,i+1)) + repeat(N_A,Y,(i,i+1)) + Δ*repeat(N_A,Z,(i,i+1)))
    end    
end


function correction!(blks::Vector{AbstractBlock}, i::Integer, r::Integer, set::Settings_XXZ)
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
