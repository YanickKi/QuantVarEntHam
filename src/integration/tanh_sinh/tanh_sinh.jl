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

"""
    TanhSinh

Struct containing the settings for the tanh-sinh quadrature.
    
# Fields 
-`atol::Float64`: absolute tolerance
-`rtol::Float64`: relative tolerance 
-`maxlevel::Int64`: maximum amount of repititions 
-`h0::Float64`: initial step size for trapezoidal integration
"""
struct TanhSinh <: AbstractIntegrator
    atol::Float64
    rtol::Float64
    maxlevel::Int64
    h0::Float64
end 

"""
    TanhSinh

Outer constructor for [`TanhSinh`](@ref) with recommended default values.
    
# Keyword arguments 
-`atol::Real=0.0`: absolute tolerance
-`rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64))`: relative tolerance 
-`maxlevel::Integer=12`: maximum amount of repititions 
-`h0::Real=1`: initial step size for trapezoidal integration
"""
function TanhSinh(; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Integer=12, h0::Real=1)
    
    return TanhSinh(Float64(atol), Float64(rtol), maxlevel, Float64(h0))
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

include("mapsum.jl")
include("mapsum_vector.jl")
include("tanh-sinh_scalar.jl")
include("tanh-sinh_vector.jl")