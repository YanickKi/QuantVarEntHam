mutable struct BufferOnlyQCFL{T<:AbstractMatrix}
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
    frechTimesBlock::Matrix{ComplexF64}
end 

function make_qcfl_only_buffer(d::Integer, numBlocks::Integer, observables::Vector{<:AbstractMatrix}) 
    return BufferOnlyQCFL{eltype(observables)}(
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks+1),
    zeros(Float64, numBlocks),
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d),
    zeros(ComplexF64, d, d),
    zeros(ComplexF64, d, d),
    zero(observables[1]),
    zeros(ComplexF64, d, d), 
    zeros(ComplexF64, d, d),
    )
end


struct QCFL_buffer
    qcfl_buff::BufferOnlyQCFL
    exp_buff::Exp_frech_buffer{ComplexF64}
end 

function make_QCFL_buffer(set::Settings, numBlocks::Integer, observables::Vector{<:AbstractMatrix})
    d = size(set.ρ_A)[1] # Hilbert space dimension
    return QCFL_buffer(make_qcfl_only_buffer(d, numBlocks, observables), make_exp_frech_buffer(ComplexF64, d))
end 


struct Commutator_buffer
    comm::Matrix{ComplexF64}
    temp::Matrix{ComplexF64}
    H_A::Matrix{ComplexF64}
end 


function make_commutator_buffer(set::Settings)
    d = size(set.ρ_A)[1]
    return Commutator_buffer(zeros(ComplexF64, d, d), zeros(ComplexF64, d, d), zeros(ComplexF64, d, d))
end 


struct Relative_entropy_buffer
    H_A::Matrix{ComplexF64}
    H_A_forexp::Matrix{ComplexF64}
    ρ_A_times_H_A::Matrix{ComplexF64}
    exp_buff::Exp_buffer{ComplexF64}
end 

function make_relative_entropy_buffer(set::Settings)

    d = size(set.ρ_A)[1]
    return Relative_entropy_buffer(zeros(ComplexF64, d, d), zeros(ComplexF64, d, d), zeros(ComplexF64, d, d), 
            make_exp_buffer(ComplexF64, d))
end 