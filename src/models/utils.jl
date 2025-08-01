function partial_trace_pure_state(
    ψ::Vector{<:Number}, dimA::Int, dimB::Int; trace_subsystem::Symbol=:B
)
    ψ_tensor = reshape(ψ, (dimA, dimB))

    if trace_subsystem == :B
        return ψ_tensor * ψ_tensor'
    elseif trace_subsystem == :A
        return transpose(ψ_tensor) * conj(ψ_tensor)
    else
        error("trace_subsystem must be :A or :B")
    end
end

function expect(Op::AbstractMatrix, ρ::AbstractMatrix)
    E = dot(Op, ρ)

    if abs(imag(E)) > 10*eps(Float64)
        error("Imaginary part $(imag(E)) too large for expectation value!")
    end

    return real(E)
end

""" 
    H_A(ansatz::AbstractAnsatz, g::Vector{<:Real}; digits::Int = 16)

Return the variational ansatz for a given parameter set `g` as a [`Block`](@ref).

The number of digits in the prefactors can be set via `digits.`

# Example 

```jldoctest
julia> using QuantVarEntHam

julia> model = TFIM(8,4,1);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix

julia> ansatz = H_A_BW(model)
H_A_BW
Number of blocks: 4

Spin 1//2
Number of spins: 4

Block 1: 
	-X₁ - 0.5*Z₁⊗ Z₂

Block 2: 
	-X₂ - 0.5*Z₁⊗ Z₂ - 0.5*Z₂⊗ Z₃

Block 3: 
	-X₃ - 0.5*Z₂⊗ Z₃ - 0.5*Z₃⊗ Z₄

Block 4: 
	-X₄ - 0.5*Z₃⊗ Z₄

julia> g = [1,2,3,4];

julia> H_A_var = H_A(ansatz, g)
Block
Spin 1//2
Number of spins: 4

-X₁ - 1.5*Z₁⊗ Z₂ - 2*X₂ - 2.5*Z₂⊗ Z₃ - 3*X₃ - 3.5*Z₃⊗ Z₄ - 4*X₄
```
"""
function H_A(ansatz::AbstractAnsatz, g::Vector{<:Real}; digits::Int=16)
    @assert length(ansatz.blocks) == length(g) "The number of blocks and parameters need to be equal!"
    return fuse_duplicate_pauli_strings(sum(ansatz.blocks .* round.(g; digits=digits)))
end

function fuse_duplicate_pauli_strings(block::Block{S,N_A}) where {S,N_A}
    pauli_strings = PauliString{S,N_A}[]
    prefactors = Float64[]

    for index in eachindex(block.pauli_strings)
        block.pauli_strings[index] ∈ pauli_strings && continue
        indices = findall(ps -> ps == block.pauli_strings[index], block.pauli_strings)
        push!(pauli_strings, block.pauli_strings[index])
        p = sum(block.prefactors[indices])
        push!(prefactors, p)
    end

    return Block{S,N_A}(prefactors, pauli_strings)
end
