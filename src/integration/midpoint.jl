struct MidPoint_scalar <: AbstractScalarIntegrator
    dt::Float64
end 
struct MidPoint_vector <: AbstractVectorIntegrator
    dt::Float64
end 

buffertrait(::MidPoint_vector) = NoNeedBuffer()

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

Outer Constructor for [`MidPoint`](@ref).

# Arguments 

-`dt::Real`: step size 

# Example 

Integrating ``\\sin (t)`` and `` (\\sin (t), \\cos (t)) `` via midpoint rule with `dt=1e-2` from `0` to `T_max = pi`

Create [`MidPoint`](@ref) struct for scalar and vector integration with `dt=1e-2`

```julia
julia> mp = MidPoint(1e-2)
MidPoint(0.01, QuantVarEntHam.var"#19#21"{Float64}(0.01), QuantVarEntHam.var"#20#22"{Float64}(0.01))
```

Now integrate the scalar function 
```julia
julia> mp.scalar_integrate(t -> sin(t), pi)
2.0000070650798945
``` 

which yields approximately 2 as one would expect.
The vector integration needs a preallocated vector to save the result in

```julia
julia> I = zeros(2)
2-element Vector{Float64}:
 0.0
 0.0

julia> mp.vector_integrate(I, t -> [sin(t), cos(t)], pi)
2-element Vector{Float64}:
 2.0000070650798945
 0.0015926595525598975
``` 
The result is returned and saved in the vector `I` aswell 

```julia
julia> I
2-element Vector{Float64}:
 2.0000070650798945
 0.0015926595525598975
```
"""
function midpoint(dt::Real)
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