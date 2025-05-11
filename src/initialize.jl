using Parameters

mutable struct Buffers{T<:AbstractMatrix}
    C_G::Vector{Float64}
    C_G_result::Vector{Float64}
    Σ::Vector{Float64}
    G_buffer::Vector{Float64}
    ρ_A_evolved::Matrix{ComplexF64}
    ρ_A_right::Matrix{ComplexF64}
    dAforpb::Matrix{ComplexF64}
    H_A::Matrix{ComplexF64}
    H_A_forexp::Matrix{ComplexF64}
    sumobs::T
    evolob::Matrix{ComplexF64}
    adjtimesBlockMatrices::Vector{Matrix{ComplexF64}}
    count::Int
end 

function create_buffers(d::Integer, numBlocks::Integer, numObservables::Integer, test_sumobs_type::T) where T<:AbstractMatrix
    return Buffers{T}(
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks),
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d),
    zeros(ComplexF64, d, d),
    zeros(ComplexF64, d, d),
    zero(test_sumobs_type),
    zeros(ComplexF64, d, d), 
    [zeros(ComplexF64, d, d) for i in 1:numBlocks],
    0
    )
end

#=
"""
    initialize(Model::Settings, H_A::Function)
    
WIP
"""
function initialize(Model::Settings, H_A::Function)
    set = Model
    @unpack N_A, observables, r_max, mtrxObs = set
    if r_max > N_A-1
    @warn "You entered r_max=$r_max but you have only $N_A sites (r_max needs to be less than N_A). This will still work and give you the maximum order of corrections."
    end
    blks = H_A(Model)

    numBlocks = length(blks.blocks)
    numObservables = length(observables)

    test_sumobs_type = sum(rand()*mtrxObs[i] for i in eachindex(mtrxObs))


    return Init{typeof(set), typeof(test_sumobs_type)}(
        set,
        blks,
        create_buffers(2^N_A, numBlocks, numObservables, test_sumobs_type),
        ExpBuffer(ComplexF64, 2^N_A)
    )
end

=#
function initialize(set::Settings, H_A::Function)
    @unpack N_A, r_max, S, mtrxObs = set
    if r_max > N_A-1
    @warn "You entered r_max=$r_max but you have only $N_A sites (r_max needs to be less than N_A). This will still work and give you the maximum order of corrections."
    end
    blks = H_A(set)

    d = Int(2*S+1)^N_A # composite Hilbert space dimension

    numBlocks = length(blks.blocks)
    numObservables = length(mtrxObs)

    test_sumobs_type = sum(rand()*mtrxObs[i] for i in eachindex(mtrxObs))


    return Init{typeof(set), typeof(test_sumobs_type)}(
        set,
        blks,
        create_buffers(d, numBlocks, numObservables, test_sumobs_type),
        ExpBuffer(ComplexF64, d)
    )
end

"""
    Init{T<:Settings, S<:AbstractMatrix}

WIP
"""
struct Init{T<:Settings, S<:AbstractMatrix}
    set::T
    blks::H_A_Var
    buff::Buffers{S}
    exp_buf::ExpBuffer{ComplexF64}
end 

#=

function initialize(Model::Settings_kitaev, H_A::Function)
    set = Model
    @unpack N_A, observables, mtrxObs = set

    blks = H_A(Model)
    println("hier")
    numBlocks = length(blks.blocks)
    numObservables = length(observables)

    test_sumobs_type = sum(rand()*mtrxObs[i] for i in eachindex(mtrxObs))


    return Init{typeof(set), typeof(test_sumobs_type)}(
        set,
        blks,
        create_buffers(N_A, numBlocks, numObservables, test_sumobs_type),
        ExpBuffer(ComplexF64, 2^N_A)
    )
end



function initialize(Model::Settings_toric, H_A::Function)
    set = Model
    @unpack N_A, observables, mtrxObs = set

    blks = H_A(Model)

    numBlocks = length(blks.blocks)
    numObservables = length(observables)

    test_sumobs_type = sum(rand()*mtrxObs[i] for i in eachindex(mtrxObs))


    return Init{typeof(set), typeof(test_sumobs_type)}(
        set,
        blks,
        create_buffers(N_A, numBlocks, numObservables, test_sumobs_type),
        ExpBuffer(ComplexF64, 2^N_A)
    )
end


=#