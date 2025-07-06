export Integrator, MidPoint, TanhSinh



abstract type AbstractScalarIntegrator end 


abstract type AbstractVectorIntegrator end


""" 
    Integrator{S<:AbstractScalarIntegrator,V<:AbstractVectorIntegrator}
    
Contains the scalar and vector integration method.
Allows unified handling of scalar and vector integration and can be conveniently passed to cost functions.
"""
struct Integrator{S<:AbstractScalarIntegrator,V<:AbstractVectorIntegrator}
    scalar_integrate::S
    vector_integrate::V
end 

abstract type AbstractIntegrator end 

include("traits.jl")
include("tanh_sinh/tanh_sinh.jl")
include("midpoint.jl")