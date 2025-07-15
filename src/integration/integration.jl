export Integrator, MidPoint, TanhSinh
export AbstractIntegrator, AbstractScalarIntegrator, AbstractVectorIntegrator


#= 
    AbstractScalarIntegrator
    
Asbtract type for scalar integration.
Each concrete type of it needs to be callable.
=#
abstract type AbstractScalarIntegrator end 

#= 
    AbstractScalarIntegrator
    
Asbtract type for vector integration.
Each concrete type of it needs to be callable.
=#
abstract type AbstractVectorIntegrator end


#= 
    Integrator{S<:AbstractScalarIntegrator,V<:AbstractVectorIntegrator}

Allows unified handling of scalar and vector integration.
=#
struct Integrator{S<:AbstractScalarIntegrator,V<:AbstractVectorIntegrator}
    scalar_integrate::S
    vector_integrate::V
end 

""" 
    AbstractIntegrator
    
Abstract type for integration.

Concrete types contain the settings (e.g. step size, tolerances, etc...).
"""
abstract type AbstractIntegrator end 

include("traits.jl")
include("tanh_sinh/tanh_sinh.jl")
include("midpoint.jl")
include("printing.jl")