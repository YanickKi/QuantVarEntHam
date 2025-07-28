function print_model(io::IO, model::AbstractModel{S,N_A}) where {S,N_A}
    println()
    println(io, "Spin $S")
    println(io, "Number of spins in the composite system N =  ", model.N)
    println(io, "Number of spins in the subsystem N_A = ", N_A)
    println(io, "Boundary conditions: ", model.periodic ? "periodic" : "open")
    print(io, "Global prefactor J = ", float_to_int(model.J))
end

function Base.show(io::IO, ::MIME"text/plain", model::TFIM{S,N_A}) where {S,N_A}
    print_model_name(io, model)
    println()
    print_model(io::IO, model)
    println()
    print(io, "Transverse field strength Γ = ", float_to_int(model.Γ))
end

function Base.show(io::IO, ::MIME"text/plain", model::XXZ{S,N_A}) where {S,N_A}
    print_model_name(io, model)
    println()
    print_model(io::IO, model)
    println()
    print(io, "Anisotropy = ", float_to_int(model.Δ))
end

function Base.show(io::IO, ::MIME"text/plain", model::Pollmann{S,N_A}) where {S,N_A}
    print_model_name(io, model)
    println()
    print_model(io::IO, model)
    println()
    println(io, "Heisenberg coupling J_Heis = ", float_to_int(model.J_Heis))
    println(io, "Transverse field strength B_x = ", float_to_int(model.Bx))
    print(io, "Square term prefactor Uzz = ", float_to_int(model.Uzz))
end

#Base.show(io::IO, model::TFIM) = show(io, ::MIME"text/plain", model)

function print_model_name(io::IO, ::TFIM)
    print(io, "TFIM")
end

function print_model_name(io::IO, ::XXZ)
    print(io, "XXZ")
end

function print_model_name(io::IO, ::Pollmann)
    print(io, "Pollmann")
end

function print_ham_params(io::IO, model::TFIM)
    print(io, "Γ=", float_to_int(model.Γ))
end

function print_ham_params(io::IO, model::XXZ)
    print(io, "Δ=", float_to_int(model.Δ))
end

function print_ham_params(io::IO, model::Pollmann)
    print(
        io,
        "J_Heis=",
        float_to_int(model.J_Heis),
        ", Bx=",
        float_to_int(model.Bx),
        ", Uzz=",
        float_to_int(model.Uzz),
    )
end

Base.show(io::IO, model::AbstractModel) = show(io, MIME"text/plain"(), model)

function Base.show(io::IO, ::MIME"text/plain", ansatz::AbstractAnsatz)
    print_ansatz_name(io, ansatz)
    println(io)
    ansatz.r_max > 1 && println(io, "r_max = $(ansatz.r_max)")
    show(io, MIME"text/plain"(), ansatz.blocks)
end

Base.show(io::IO, ansatz::AbstractAnsatz) = show(io, MIME"text/plain"(), ansatz)
