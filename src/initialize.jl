using Parameters


"""
    initialize(Model::Settings, H_A::Function)
    
Convenient constructor for [Init](@ref). 
Automatically pre allocates the correct buffer and constructs the correct blocks for the given model `Model` and given variational Ansatz `H_A`  

# Arguments

- `Model::Settings`: struct containing the settings for the model 
- `H_A::Function`: the variational Ansatz

# Example
Constructing the settings, buffer and blocks for the TFIM with ``N=8``, ``N_\\text{A}=4``, OBC, ``Γ=1`` and ``T_\\text{max} = 2`` with
the variational Ansatz ``H_\\text{A}^\\text{BW}``.

```julia
init = initialize(TFIM(8, 4, 1, 2), H_A_BW)
```

"""
function initialize(set::Settings, H_A::Function)
    @unpack N_A, r_max, S, observables = set
    if r_max > N_A-1
    @warn "You entered r_max=$r_max but you have only $N_A sites (r_max needs to be less than N_A). This will still work and give you the maximum order of corrections."
    end
    blocks = H_A(set)

    d = Int(2*S+1)^N_A # composite Hilbert space dimension

    numBlocks = length(blocks)

    test_sumobs_type = sum(rand()*observables[i] for i in eachindex(observables))


    return Init{typeof(set), typeof(test_sumobs_type)}(
        set,
        blocks,
        create_buffers(d, numBlocks, test_sumobs_type),
        ExpBuffer(ComplexF64, d)
    )
end

"""
    Init{T<:Settings, S<:AbstractMatrix}

Struct containing all settings, buffers and blocks for the variational Ansatz.

!!! tip
    Use the constructor [initialize](@ref) for constructing this struct, since it automatically constructs the correct buffers and blocks for the 
    variational Ansatz.

# Fields 

- `set::T`: settings of the model (see [Settings](@ref))
- `blocks::Vector{Matrix{ComplexF64}}`: the blocks of the variational Ansatz
- `buff::Buffers{S}`: buffer for minimizing runtime during cost function computation (not public)
- `exp_buf::ExpBuffer{ComplexF64}`: buffer for minimizing runtime during computation of matrix exponential and its Frechet derivative (not public)
"""
struct Init{T<:Settings, S<:AbstractMatrix}
    set::T
    blocks::Vector{Matrix{ComplexF64}}
    buff::BufferOnlyQCFL{S}
    exp_buf::BufferExp{ComplexF64}
end