""" 
    Commutator_buffer

Struct containing the buffers for the cost function [`Commutator`](@ref).
"""
struct Commutator_buffer
    comm::Matrix{ComplexF64}
    temp::Matrix{ComplexF64}
    H_A::Matrix{ComplexF64}
end 


""" 
    Commutator_buffer

Outer constructor for [`Commutator_buffer`](@ref) given a `model`.

# Arguments

-`model::AbstractModel`: model
"""
function Commutator_buffer(model::AbstractModel)
    d = size(model.ρ_A)[1]
    return Commutator_buffer(zeros(ComplexF64, d, d), zeros(ComplexF64, d, d), zeros(ComplexF64, d, d))
end 

"""
    Commutator

Commutator as a cost function defined as 
```math
C(\vec{g})  = \frac{||[H_\text{A}^\text{Var}(\vec{g}), ρ_\text{A}]||}{||H_\text{A}^\text{Var}(\vec{g}||||ρ_\text{A}||},
```
where ``H_\\text{A}^\\text{Var}(\\vec{g})`` is the variational Ansatz and ``ρ_\\text{A}`` the exact reduced density matrix. 
The Frobenius norm is used.

# Fields 

-`model::M`: model 
-`blocks::Vector{Matrix{ComplexF64}}`: blocks for the variational Ansatz 
-`buff::Commutator_buffer`: buffers 
"""
struct Commutator{M<:AbstractModel} <: AbstractCostFunction
    model::M
    blocks::Vector{Matrix{ComplexF64}}
    buff::Commutator_buffer
end 

"""
    Commutator

Outer Constructor for [`Commutator`](@ref) s.t. the correct buffers will be automatically constructed 
for the given `model` and `blocks`.  

# Required Arguments 

-`model::AbstractModel`: model 
-`blocks::Vector{<:AbstractMatrix}`: blocks for the variational Ansatz 

# Keyword Arguments

-`buffer::Union{Nothing, Commutator_buffer} = nothing`
"""
function Commutator(model::AbstractModel, blocks::Vector{<:AbstractMatrix}; buffer::Union{Nothing, Commutator_buffer} = nothing)
    buffer = something(buffer, Commutator_buffer(model))
    return Commutator(
        model, 
        blocks, 
        buffer
    )
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
