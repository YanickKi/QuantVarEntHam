struct Commutator_buffer
    comm::Matrix{ComplexF64}
    temp::Matrix{ComplexF64}
    H_A::Matrix{ComplexF64}
end 

struct Commutator{M} <: AbstractCostFunction
    model::M
    blocks::Vector{Matrix{ComplexF64}}
    buff::Commutator_buffer
end 

function Commutator(model::AbstractModel, blocks::Vector{<:AbstractMatrix}; buffer::Union{Nothing, Commutator_buffer} = nothing)
    buffer = something(buffer, Commutator_buffer(model))
    return Commutator(
        model, 
        blocks, 
        buffer
    )
end 

function Commutator_buffer(model::AbstractModel)
    d = size(model.ρ_A)[1]
    return Commutator_buffer(zeros(ComplexF64, d, d), zeros(ComplexF64, d, d), zeros(ComplexF64, d, d))
end 

function (c::Commutator)(g::Vector{<:Real})
    
    comm = c.buff.comm
    H_A = c.buff.H_A
    ρ_A = c.model.ρ_A 

    get_H_A!(c, g)

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm
    
    return norm(comm)/(2*norm(H_A)*norm(ρ_A))
end 


function _gradient!(c::Commutator, G::Vector{<:Real}, g::Vector{<:Real}, free_indices)
    
    get_H_A!(c, g)

    Λ_comm, Λ_H_A, Λ_ρ_A = pre_computations_gradient(c)

    @fastmath @inbounds @simd for i in eachindex(free_indices)
        G[i] = gradient_component(c, Λ_comm, Λ_H_A, Λ_ρ_A, free_indices[i])
    end 
    return Λ_comm/(2*Λ_H_A * Λ_ρ_A)
end 



function pre_computations_gradient(c::Commutator)
    
    comm = c.buff.comm
    H_A = c.buff.H_A
    ρ_A = c.model.ρ_A

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm


    Λ_comm = norm(comm)
    Λ_H_A = norm(H_A)
    Λ_ρ_A = norm(ρ_A)

    return Λ_comm, Λ_H_A, Λ_ρ_A 
end 


function gradient_component(c::Commutator, Λ_comm::Real, Λ_H_A::Real, Λ_ρ_A::Real, index::Integer)
    mul!(mul!(c.buff.temp, c.model.ρ_A, c.blocks[index]), c.blocks[index], c.model.ρ_A, 1, -1) # compute [h_index, ρ_A] and save it in temp
    ∂Λ_comm = 1/Λ_comm * real(sum(had!(c.buff.temp, c.buff.temp, c.buff.comm)))
    ∂Λ_H_A = 1/Λ_H_A * real(sum(had!(c.buff.temp, c.buff.H_A, c.blocks[index])))
    return 1/(2*Λ_ρ_A*Λ_H_A^2)*(∂Λ_comm * Λ_H_A - Λ_comm * ∂Λ_H_A) 
end 

# compute hadamard product of two square matrices A and B and save it in buff
function had!(buff::AbstractMatrix, A::AbstractMatrix,B::AbstractMatrix)
    n = size(A)[1]
    @fastmath @inbounds @simd for j in 1:n
       for i in 1:n
         @inbounds buff[i,j] = A[i,j] *B[i,j]
       end
    end
    return buff
end
