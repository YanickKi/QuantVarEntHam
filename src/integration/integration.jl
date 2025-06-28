abstract type AbstractVectorIntegrator end
abstract type AbstractScalarIntegrator end 


""" 
    Integrator{S,V}
    
Abstract type for integration.
"""
struct Integrator{S<:AbstractScalarIntegrator,V<:AbstractVectorIntegrator}
    scalar_integrate::S
    vector_integrate::V
end 

abstract type BufferTrait end 
struct NeedBuffer <: BufferTrait end 
struct NoNeedBuffer <: BufferTrait end

shorten_buffer!(::NoNeedBuffer,::AbstractVectorIntegrator, ::AbstractVector) = nothing

function shorten_buffer!(::NeedBuffer, vector_integrate::AbstractVectorIntegrator, fixed_indices::AbstractVector)
    for _ in eachindex(fixed_indices)
            pop!(vector_integrate.buffer)
    end 
end 

export Integrator, midpoint, Tanh_sinh

include("tanh_sinh.jl")
include("midpoint.jl")