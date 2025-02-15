using Parameters

mutable struct Buffers{T<:AbstractMatrix}
    C_G::Vector{Float64}
    C_G_result::Vector{Float64}
    G_buffer::Vector{Float64}
    dev::Vector{Float64}
    ρ_A_evolved::Matrix{ComplexF64}
    ρ_A_right::Matrix{ComplexF64}
    dAforpb::Matrix{ComplexF64}
    H_A::Matrix{ComplexF64}
    H_A_forexp::Matrix{ComplexF64}
    sumobs::T
    evolobs::Vector{Matrix{ComplexF64}}
    adjtimesBlockMatrices::Vector{Matrix{ComplexF64}}
    count::Int
end 

function create_buffers(N_A::Integer, numBlocks::Integer, numObservables::Integer, test_sumobs_type::T) where T<:AbstractMatrix
    d = 2^N_A
    return Buffers{T}(
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks),
    zeros(Int, numObservables),
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d),
    zeros(ComplexF64, d, d),
    zeros(ComplexF64, d, d),
    zero(test_sumobs_type),
    [zeros(ComplexF64, d, d) for i in 1:numObservables], 
    [zeros(ComplexF64, d, d) for i in 1:numBlocks],
    0
    )
end

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

    q = integration_tables(maxlevel = 12)

    test_sumobs_type = sum(rand()*mtrxObs[i] for i in eachindex(mtrxObs))

    return Init{typeof(set), typeof(test_sumobs_type)}(
        set,
        blks,
        q,
        create_buffers(N_A, numBlocks, numObservables, test_sumobs_type)
    )
end

"""
    Init{T<:Settings, S<:AbstractMatrix}

WIP
"""
struct Init{T<:Settings, S<:AbstractMatrix}
    set::T
    blks::H_A_Var
    q::QuadTS{12}
    buff::Buffers{S}
end 
