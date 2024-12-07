using Parameters

#=
#### backup original 

mutable struct Buffers 
    C_G::Vector{Float64}
    C_G_result::Vector{Float64}
    G_buffer::Vector{Float64}
    dev::Vector{Float64}
    ρ_A_evolved::Matrix{ComplexF64}
    ρ_A_right::Matrix{ComplexF64}
    dAforpb::Matrix{ComplexF64}
    sumobs::Matrix{ComplexF64}
    evolobs::Vector{Matrix{ComplexF64}}
    adjtimesBlockMatrices::Vector{Matrix{ComplexF64}}
end 

function create_buffers(N_A::Integer, numBlocks::Integer, numObservables::Integer)
    d = 2^N_A
    return Buffers(
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks),
    zeros(Int, numObservables),
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d),
    zeros(ComplexF64, d, d), 
    [zeros(ComplexF64, d, d) for i in 1:numObservables], 
    [zeros(ComplexF64, d, d) for i in 1:numBlocks], 
    )
end

function initialize(Model::Settings, H_A::Function)
    set = Model
    @unpack N_A, observables = set 
    blks = H_A(Model)

    numBlocks = length(blks.blocks)
    numObservables = length(observables)

    q = integration_tables()
    buff = create_buffers(N_A, numBlocks, numObservables)

    return Init{typeof(set)}(
        set,
        blks,
        q,
        buff
    )
end

struct Init{T<:Settings}
    set::T
    blks::H_A_Var
    q::QuadTS{12}
    buff::Buffers
end 



=#

#### ohne parametrisch 
#=
mutable struct Buffers 
    C_G::Vector{Float64}
    C_G_result::Vector{Float64}
    G_buffer::Vector{Float64}
    dev::Vector{Float64}
    ρ_A_evolved::Matrix{ComplexF64}
    ρ_A_right::Matrix{ComplexF64}
    dAforpb::Matrix{ComplexF64}
    sumobs::Diagonal{ComplexF64, Vector{ComplexF64}}
    evolobs::Vector{Matrix{ComplexF64}}
    adjtimesBlockMatrices::Vector{Matrix{ComplexF64}}
end 

function create_buffers(N_A::Integer, numBlocks::Integer, numObservables::Integer)
    d = 2^N_A
    return Buffers(
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks),
    zeros(Int, numObservables),
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d),
    Diagonal(zeros(d)),
    [zeros(ComplexF64, d, d) for i in 1:numObservables], 
    [zeros(ComplexF64, d, d) for i in 1:numBlocks], 
    )
end

function initialize(Model::Settings, H_A::Function)
    set = Model
    @unpack N_A, observables = set 
    blks = H_A(Model)

    numBlocks = length(blks.blocks)
    numObservables = length(observables)

    q = integration_tables()
    buff = create_buffers(N_A, numBlocks, numObservables)

    return Init{typeof(set)}(
        set,
        blks,
        q,
        buff
    )
end

struct Init{T<:Settings}
    set::T
    blks::H_A_Var
    q::QuadTS{12}
    buff::Buffers
end 

=#

#### mit parametrisch 


mutable struct Buffers{T<:AbstractMatrix}
    C_G::Vector{Float64}
    C_G_result::Vector{Float64}
    G_buffer::Vector{Float64}
    dev::Vector{Float64}
    ρ_A_evolved::Matrix{ComplexF64}
    ρ_A_right::Matrix{ComplexF64}
    dAforpb::Matrix{ComplexF64}
    H_A::Matrix{ComplexF64}
    sumobs::T
    evolobs::Vector{Matrix{ComplexF64}}
    adjtimesBlockMatrices::Vector{Matrix{ComplexF64}}
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
    zero(test_sumobs_type),
    [zeros(ComplexF64, d, d) for i in 1:numObservables], 
    [zeros(ComplexF64, d, d) for i in 1:numBlocks], 
    )
end

function initialize(Model::Settings, H_A::Function)
    set = Model
    @unpack N_A, observables, r_max, mtrxObs = set
    if r_max > N_A-1
    @warn "You entered r_max=$r_max but you have only $N_A sites (r_max needs to be less than N_A). This will still work and give you the maximum order of corrections."
    end
    blks = H_A(Model)

    numBlocks = length(blks.blocks)
    numObservables = length(observables)

    q = integration_tables()

    test_sumobs_type = sum(rand()*mtrxObs[i] for i in eachindex(mtrxObs))

    return Init{typeof(set), typeof(test_sumobs_type)}(
        set,
        blks,
        q,
        create_buffers(N_A, numBlocks, numObservables, test_sumobs_type)
    )
end

struct Init{T<:Settings, S<:AbstractMatrix}
    set::T
    blks::H_A_Var
    q::QuadTS{12}
    buff::Buffers{S}
end 
