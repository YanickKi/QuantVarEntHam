function Base.show(io::IO, cost::QCFL{S,N_A}) where {S,N_A}
    println(io, "QCFL")
    println(io)
    print_model_short(io, cost.model)
    println(io)
    print_BW_BWV(io, cost.str_blocks)
    println(io)
    print(io, "Integration method: ")
    print_integration_name(io, cost.integrator)
    print_integration_settings_short(io, cost.integrator.scalar_integrate)
    println(io)
    print(io, "T_max=", float_to_int(cost.T_max))
end 

function Base.show(io::IO, cost::Commutator)
    println(io, "Commutator")
    println(io)
    print_model_short(io, cost.model)
    println(io)
    print_BW_BWV(io, cost.str_blocks)
end 

function Base.show(io::IO, cost::RelativeEntropy)
    println(io, "Relative entropy")
    println(io)
    print_model_short(io, cost.model)
    println(io)
    print_BW_BWV(io, cost.str_blocks)
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


function print_blocks(cost::AbstractFreeCostFunction)
    show(stdout, MIME"text/plain"(),cost.str_blocks)
end 

function print_observables(cost::AbstractFreeCostFunction)
    show(stdout, MIME"text/plain"(),cost.str_observables)
end 

print_blocks(fc::FixedCost) = unwrap(print_blocks, fc)
print_observables(fc::FixedCost) = unwrap(print_observables, fc)

function Base.show(io::IO, fc::FixedCost)
    print(io, "Fixed cost: ")
    show(io, fc.c)
    println(io,)
    println(io,)
    println(io, "Fixed indices: ", fc.fixed_indices)
    println(io, "Fixed values: ", float_to_int.(fc.fixed_values))
end 