"""
    print_blocks(cost::AbstractCostFunction)

Print the blocks of a given cost function. 
"""
print_blocks(cost::AbstractCostFunction) = print_blocks(cost)
print_blocks(cost::AbstractFreeCostFunction) = show(stdout, MIME"text/plain"(),cost.blocks)
print_blocks(fc::FixedCost) = unwrap(print_blocks, fc)


"""
    print_observables(cost::AbstractCostFunction)

Print the observables of a given cost function.
"""
print_observables(cost::AbstractCostFunction) = print_observables(cost)
print_observables(cost::AbstractFreeCostFunction) = show(stdout, MIME"text/plain"(),cost.observables)
print_observables(fc::FixedCost) = unwrap(print_observables, fc)


"""
    print_model(cost::AbstractCostFunction)

Print the model of a given cost function
"""
print_model(cost::AbstractCostFunction) = print_model(cost)
print_model(cost::AbstractFreeCostFunction) = show(stdout,  MIME"text/plain"(), cost.model)
print_model(fc::FixedCost) = unwrap(print_model, fc)


function Base.show(io::IO, ::MIME"text/plain", cost::QCFL{S,N_A}) where {S,N_A}
    println(io, "QCFL")
    println(io)
    print_model_short(io, cost.model)
    println(io)
    print_BW_BWV(io, cost.blocks)
    println(io)
    print(io, "Integration method: ")
    print_integration_name(io, cost.integrator)
    print_integration_settings_short(io, cost.integrator.scalar_integrate)
    println(io)
    print(io, "T_max=", float_to_int(cost.T_max))
end 

function Base.show(io::IO, ::MIME"text/plain", cost::Commutator)
    println(io, "Commutator")
    println(io)
    print_model_short(io, cost.model)
    println(io)
    print_BW_BWV(io, cost.blocks)
end 

function Base.show(io::IO, ::MIME"text/plain", cost::RelativeEntropy)
    println(io, "Relative entropy")
    println(io)
    print_model_short(io, cost.model)
    println(io)
    print_BW_BWV(io, cost.blocks)
end 

function print_model_short(io::IO, model::AbstractModel{S,N_A}) where {S,N_A}
    print(io, "Model: ")
    print_model_name(io, model)
    bc = (model.periodic == true ? "PBC" : "OBC")
    print(io, " (S=$S, N=$(model.N), $bc, N_A = $N_A, J=$(float_to_int(model.J)), ")
    print_ham_params(io, model)
    print(io, ")")
end 

function print_BW_BWV(io::IO, blocks::Vector{<:AbstractBlock{S,N_A}}) where {S,N_A}
    length(blocks) == N_A ? print(io, "BW Ansatz") : print(io, "BW violating Ansatz")
end


function Base.show(io::IO, ::MIME"text/plain", fc::FixedCost)
    print(io, "Fixed cost: ")
    show(io, fc.c)
    println(io,)
    println(io,)
    println(io, "Fixed indices: ", fc.fixed_indices)
    print(io, "Fixed values: ", float_to_int.(fc.fixed_values))
end 

function Base.show(io::IO, ::MIME"text/plain", ::QCFLBuffer{S,N_A,L}) where {S,N_A,L}
    println(io, "Buffer for QCFL")
    println(io)
    println(io, "Spin number ", S)
    println(io, "Number of spins in subsytem: N_A = ", N_A)
    print(io, "Length of buffer L = ", L)
end 

Base.show(io::IO, buffer::QCFLBuffer) = show(io, MIME"text/plain"(), buffer)

function Base.show(io::IO, ::MIME"text/plain", ::CommutatorBuffer{S,N_A}) where {S,N_A}
    println(io, "Buffer for Commutator")
    println(io)
    println(io, "Spin number ", S)
    print(io, "Number of spins in subsytem: N_A = ", N_A)
end 

Base.show(io::IO, buffer::CommutatorBuffer) = show(io, MIME"text/plain"(), buffer)

function Base.show(io::IO, ::MIME"text/plain", ::RelativeEntropyBuffer{S,N_A}) where {S,N_A}
    println(io, "Buffer for Relative entropy")
    println(io)
    println(io, "Spin number ", S)
    print(io, "Number of spins in subsytem: N_A = ", N_A)
end 

Base.show(io::IO, buffer::RelativeEntropyBuffer) = show(io, MIME"text/plain"(), buffer)

Base.show(io::IO, cost::AbstractCostFunction) = show(io, MIME"text/plain"(), cost)
