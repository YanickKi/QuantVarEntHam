""" 
    AbstractIntegrator
    
Abstract type for integration.
"""
abstract type AbstractIntegrator end


export AbstractIntegrator, MidPoint, Tanh_sinh

include("tanh_sinh.jl")
include("midpoint.jl")