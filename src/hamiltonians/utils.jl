function partial_trace_pure_state(ψ::Vector{<:Number}, dimA::Int, dimB::Int; trace_subsystem::Symbol = :B)
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
