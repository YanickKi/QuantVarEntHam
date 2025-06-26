"""
    MidPoint <: AbstractIntegrator

Struct for the integration via mid-point rule.

# Fields 
-`dt::Float64` : step size
-`scalar_integrate::Function`: function for integrating scalar functions 
-`vector_integrate::Function`: function for integrating vector functions

Allows unified handling of scalar and vector functions as in e.g. ['QCFL'](@ref)

# Example 


"""
struct MidPoint <: AbstractIntegrator
    dt::Float64
    scalar_integrate::Function 
    vector_integrate::Function 
end 

function MidPoint(dt::Real)
    scalar_integrate = (f, T_max) -> midpoint_scalar(f,0,T_max,dt)
    vector_integrate = (I, f, T_max) -> midpoint_vector!(I, f, 0, T_max, dt)
    return MidPoint(dt, scalar_integrate, vector_integrate) 
end 



function midpoint_scalar(f::Function, a::Real, b::Real, dt::Real)
    _a = Float64(a)
    _b = Float64(b)
    _dt = Float64(dt)
    I = f(_a+_dt/2)
    points = range(start = _a+3*_dt/2, stop = _b-_dt/2, step = _dt) 
   
    for p in points 
        I += f(p)
    end 
    return I*dt
end 

function midpoint_vector!(I::AbstractVector, f::Function, a::Real, b::Real, dt::Real)
    _a = Float64(a)
    _b = Float64(b)
    _dt = Float64(dt)
    I .= f(_a+_dt/2)
    points = range(start = _a+3*_dt/2, stop = _b-_dt/2, step = _dt) 
    for p in points 
        I .+= f(p)
    end 
    I .*= dt
end 