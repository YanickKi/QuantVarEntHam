struct MidPoint_scalar <: AbstractScalarIntegrator
    dt::Float64
end 
struct MidPoint_vector <: AbstractVectorIntegrator
    dt::Float64
end 

buffertrait(::MidPoint_vector) = NoNeedBuffer()

struct MidPoint <: AbstractIntegrator
    dt::Float64
end 
    
#=
"""
    MidPoint <: AbstractIntegrator

Struct for the integration via mid-point rule.

# Fields 
-`dt::Float64` : step size
-`scalar_integrate::Function`: function for integrating scalar functions 
-`vector_integrate::Function`: function for integrating vector functions

Allows unified handling of integrating scalar and vector functions as in e.g. ['QCFL'](@ref).

!!! tip
    Use 

"""
=#

"""
    MidPoint(dt::Real) 

Outer Constructor for [`Integrator`](@ref) to construct the scalar and vector integration via midpoint rule with a given step size `dt`

# Arguments 

-`dt::Real`: step size 
"""
function midpoint(dt::Real)
    #scalar_integrate = (f, T_max) -> midpoint_scalar(f,0,T_max,dt)
    #vector_integrate = (I, f, T_max) -> midpoint_vector!(I, f, 0, T_max, dt)
    return _ -> _midpoint(dt) 
end 

function _midpoint(dt::Real)
    #scalar_integrate = (f, T_max) -> midpoint_scalar(f,0,T_max,dt)
    #vector_integrate = (I, f, T_max) -> midpoint_vector!(I, f, 0, T_max, dt)
    return Integrator(MidPoint_scalar(dt), MidPoint_vector(dt)) 
end 

function (mp_s::MidPoint_scalar)(f::Function, b::Real)
    dt = mp_s.dt
    _a = 0.
    _b = Float64(b)
    _dt = Float64(dt)
    I = f(_a+_dt/2)
    points = range(start = _a+3*_dt/2, stop = _b-_dt/2, step = _dt) 
   
    for p in points 
        I += f(p)
    end 
    return I*dt
end 

function (mp_v::MidPoint_vector)(I::AbstractVector, f::Function, b::Real)
    dt = mp_v.dt
    _a = 0.
    _b = Float64(b)
    _dt = Float64(dt)
    I .= f(_a+_dt/2)
    points = range(start = _a+3*_dt/2, stop = _b-_dt/2, step = _dt) 
    for p in points 
        I .+= f(p)
    end 
    I .*= dt
end 