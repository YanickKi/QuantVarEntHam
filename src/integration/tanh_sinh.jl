struct QuadTS{N}
    h0::Float64
    origin::Tuple{Float64,Float64}
    table0::Vector{Tuple{Float64,Float64}}
    tables::NTuple{N,Vector{Tuple{Float64,Float64}}}
end

struct TanhSinh_scalar <: AbstractScalarIntegrator
    integration_table::QuadTS
    atol::Float64
    rtol::Float64
end 
struct TanhSinh_vector <: AbstractVectorIntegrator
    integration_table::QuadTS
    atol::Float64
    rtol::Float64
    buffer::Vector{Float64}
end 

buffertrait(::TanhSinh_vector) = NeedBuffer()

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


function tanh_sinh(length_buffer::Int; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Integer=12, h0::Float64=1.0)
    q = integration_tables(maxlevel = maxlevel, h0 = h0)
    buffer = zeros(length_buffer)
    scalar_integrate = TanhSinh_scalar(q, atol, rtol)
    vector_integrate = TanhSinh_vector(q, atol, rtol, buffer)

    return Integrator(scalar_integrate, vector_integrate)

end 

include("mapsum.jl")
include("mapsum_vector.jl")
include("tanh-sinh_scalar.jl")
include("tanh-sinh_vector.jl")