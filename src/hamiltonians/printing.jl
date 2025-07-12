function print_model(io::IO, model::AbstractModel{S,N_A}) where {S,N_A}
    println()
    println(io, "Spin $S")
    println(io, "Number of spins in the composite system N =  ", model.N)
    println(io, "Number of spins in the subsystem N_A = ", N_A)
    println(io, "Global prefactor J = ", model.J)
end 

function Base.show(io::IO, model::TFIM{S,N_A}) where {S,N_A}
    println(io, "TFIM")
    print_model(io::IO, model)
    println(io, "Transversal field strength Γ = ", model.Γ)
end 

function Base.show(io::IO, model::XXZ{S,N_A}) where {S,N_A}
    println(io, "XXZ model")
    print_model(io::IO, model)
    println(io, "Anisotropy = ", model.Δ)
end 

function Base.show(io::IO, model::Pollmann{S,N_A}) where {S,N_A}
    println(io, "Pollmann model")
    print_model(io::IO, model)
    println(io, "Heisenberg coupling J_Heis = ", model.J_Heis)
    println(io, "Transversal field strength B_x = ", model.Bx)
    println(io, "Square term prefactor = ", model.Uzz)
end 
