struct MidPoint_scalar <: AbstractScalarIntegrator
    dt::Float64
end 
struct MidPoint_vector <: AbstractVectorIntegrator
    dt::Float64
end 

buffertrait(::MidPoint_vector) = NoNeedBuffer()


"""
    MidPoint(dt::Real) 

Struct containing the settings for midpoint integration.

# Arguments 

-`dt::Real`: step size 
"""
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