struct Commutator{M} <: AbstractCostFunction
    model::M
    blocks::Vector{Matrix{ComplexF64}}
    buff::Commutator_buffer
end 

function Commutator(model::AbstractModel, blocks::Vector{<:AbstractMatrix})
    return Commutator(
        model, 
        blocks, 
        make_commutator_buffer(model)
    )
end 


function (c::Commutator)(g::Vector{<:Real})
    
    comm = c.buff.comm
    H_A = c.buff.H_A
    ρ_A = c.model.ρ_A 

    get_H_A!(H_A, g, c.blocks)

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm
    
    return norm(comm)/(2*norm(H_A)*norm(ρ_A))
end 

function gradient!(G::Vector{<:Real}, c::Commutator, g::Vector{<:Real})
    
    comm = c.buff.comm
    H_A = c.buff.H_A
    temp = c.buff.temp
    ρ_A = c.model.ρ_A
    blocks = c.blocks

    get_H_A!(H_A, g, c.blocks)

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm


    Λ_comm = norm(comm)
    Λ_H_A = norm(H_A)
    Λ_ρ_A = norm(ρ_A)

    @fastmath @inbounds @simd for index in eachindex(G)
        mul!(mul!(temp, ρ_A, blocks[index]), blocks[index], ρ_A, 1, -1) # compute [h_index, ρ_A] and save it in temp
        ∂Λ_comm = 1/Λ_comm * real(sum(had!(temp, temp, comm)))
        ∂Λ_H_A = 1/Λ_H_A * real(sum(had!(temp, H_A, blocks[index])))
        G[index] = 1/(2*Λ_ρ_A*Λ_H_A^2)*(∂Λ_comm * Λ_H_A - Λ_comm * ∂Λ_H_A) 
    end 
    return Λ_comm/(2*Λ_H_A * Λ_ρ_A)
end 

# compute hadamard product of two square matrices A and B and save it in buff
function had!(buff::AbstractMatrix, A::AbstractMatrix,B::AbstractMatrix)
    n = size(A)[1]
    for j in 1:n
       for i in 1:n
         @inbounds buff[i,j] = A[i,j] *B[i,j]
       end
    end
    return buff
end
