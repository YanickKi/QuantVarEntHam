struct QuadTS{N}
    h0::Float64
    origin::Tuple{Float64,Float64}
    table0::Vector{Tuple{Float64,Float64}}
    tables::NTuple{N,Vector{Tuple{Float64,Float64}}}
end

function integration_tables(;maxlevel::Integer=12, h0::Real=1)
    _h0 = Float64(h0)
    origin = samplepoint(0.)
    table0 = generate_table(_h0, 1)
    tables = Vector{Tuple{Float64,Float64}}[]
    for level in 1:maxlevel
        h = _h0/2^level
        table = generate_table(h, 2)
        push!(tables, table)
    end
    return QuadTS{maxlevel}(_h0, origin, table0, Tuple(tables))
end

struct Tanh_sinh <: AbstractIntegrator
    integration_table::QuadTS
    atol::Float64
    rtol::Float64
    buffer::Vector{Float64}
    scalar_integrate::Function 
    vector_integrate::Function
end

function Tanh_sinh(length_buffer::Int; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Integer=12, h0::Float64=1.0)
    q = integration_tables(maxlevel = maxlevel, h0 = h0)
    buffer = zeros(length_buffer)
    scalar_integrate = (f, T_max) -> tanh_sinh(f, 0, T_max, q, atol=atol, rtol=rtol)
    vector_integrate = (I, f, T_max) -> tanh_sinh!(I, f, 0, T_max, q, buffer, atol=atol, rtol=rtol)

    return Tanh_sinh(q, atol, rtol, buffer, scalar_integrate, vector_integrate)

end 

include("mapsum.jl")
include("mapsum_vector.jl")
include("tanh-sinh_scalar.jl")
include("tanh-sinh_vector.jl")