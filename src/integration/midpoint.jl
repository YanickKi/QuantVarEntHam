function midpoint(f::Function, a::Real, b::Real, dt::Real)
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

function midpoint!(I::AbstractVector, f::Function, a::Real, b::Real, dt::Real)
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