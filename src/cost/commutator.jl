""" 
    CommutatorBuffer{S,N_A}
    CommutatorBuffer(::AbstractModel{S,N_A}) where {S,N_A}

Buffer for the cost function [`Commutator`](@ref) for a given [`AbstractModel`](@ref).
"""
struct CommutatorBuffer{S,N_A}
    comm::Matrix{ComplexF64}
    temp::Matrix{ComplexF64}
    H_A::Matrix{ComplexF64}
end 


function CommutatorBuffer(::AbstractModel{S,N_A}) where {S,N_A}
    d = Int((2*S+1)^N_A)
    return CommutatorBuffer{S,N_A}(zeros(ComplexF64, d, d), zeros(ComplexF64, d, d), zeros(ComplexF64, d, d))
end 

"""
    Commutator{M<:AbstractModel, B<:CommutatorBuffer, A<:AbstractAnsatz} <: AbstractFreeCostFunction
    Commutator(model::AbstractModel{S,N_A}, ansatz::AbstractAnsatz{S,N_A}; buffer::Union{Nothing, CommutatorBuffer{S,N_A}} = nothing) where {S,N_A}

Contains the model, blocks, observables and buffers for the Commutator as a cost function.

Existing buffers can be provided, and are constructed automatically otherwise. 

The cost function is given by 
```math
\\mathcal{C}(\\vec{g})  = \\frac{||[H_\\text{A}^\\text{Var}(\\vec{g}), ρ_\\text{A}]||}{||H_\\text{A}^\\text{Var}(\\vec{g})||||ρ_\\text{A}||},
```
where ``H_\\text{A}^\\text{Var}(\\vec{g})`` is the variational Ansatz and ``ρ_\\text{A}`` the exact reduced density matrix. 
The Frobenius norm is used.

# Gradient 

The gradient is given by 

```math
\\partial_{g_k} \\mathcal{C}(\\vec{g}) =
\\frac{1}{2||\\rho_\\text{A}|| || H_\\text{A}^\\text{Var}(\\vec{g}) ||^2} 
\\left ( \\frac{|| H_\\text{A}^\\text{Var}(\\vec{g}) ||}{|| [H_\\text{A}^\\text{Var}(\\vec{g}), ρ_\\text{A}] ||} 
\\text{Sum} \\left ( [h_k, ρ_\\text{A}] \\, .\\!* [H_\\text{A}^\\text{Var}(\\vec{g}), ρ_\\text{A}] \\right ) -
\\frac{||[H_\\text{A}^\\text{Var}(\\vec{g}), ρ_\\text{A}] ||}{|| H_\\text{A}^\\text{Var}(\\vec{g}) ||}
\\text{Sum} \\left ( H_\\text{A}^\\text{Var}(\\vec{g}) \\, . \\! * h_k \\right )
\\right ).
```
Here, `` \\text{Sum}(X) = \\sum_{ij} X_{ij}`` denotes the sum of all matrix elements and 
the elementwise product aka Hadamard product is denoted by ``C = A \\, . \\! *B`` s.t. 
``C_{ij} = A_{ij} B_{ij}``.
"""
struct Commutator{M<:AbstractModel, B<:CommutatorBuffer, A<:AbstractAnsatz} <: AbstractFreeCostFunction
    model::M
    blocks_mat::Vector{Matrix{ComplexF64}}
    ansatz::A
    buff::B
end 

function Commutator(model::AbstractModel{S,N_A}, ansatz::AbstractAnsatz{S,N_A}; buffer::Union{Nothing, CommutatorBuffer{S,N_A}} = nothing) where {S,N_A}
    buffer = @something buffer CommutatorBuffer(model)
    
    blocks_mat = Matrix.(mat.(ansatz.blocks))

    return Commutator(
        model, 
        blocks_mat, 
        ansatz,
        buffer
    )
end 

# calc_cost(cost::Commuator, g::Vector{float64})

function (c::Commutator)(g::Vector{Float64})
    
    comm = c.buff.comm
    H_A = c.buff.H_A
    ρ_A = c.model.ρ_A 

    get_H_A!(c, g)

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm
    
    return norm(comm)/(2*norm(H_A)*norm(ρ_A))
end 

(c::Commutator)(g::Vector{<:Real}) = c(Float64.(g))

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
    mul!(mul!(c.buff.temp, c.model.ρ_A, c.blocks_mat[index]), c.blocks_mat[index], c.model.ρ_A, 1, -1) # compute [h_index, ρ_A] and save it in temp
    ∂Λ_comm = 1/Λ_comm * real(sum(had!(c.buff.temp, c.buff.temp, c.buff.comm)))
    ∂Λ_H_A = 1/Λ_H_A * real(sum(had!(c.buff.temp, c.buff.H_A, c.blocks_mat[index])))
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
