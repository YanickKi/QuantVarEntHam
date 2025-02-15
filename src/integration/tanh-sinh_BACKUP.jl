using LinearAlgebra

struct QuadTS{N}
    h0::Float64
    origin::Tuple{Float64,Float64}
    table0::Vector{Tuple{Float64,Float64}}
    tables::NTuple{N,Vector{Tuple{Float64,Float64}}}
end

function integration_tables(;maxlevel::Integer=12, h0::Float64=1.)
    origin = samplepoint(0.)
    table0 = generate_table(h0, 1)
    tables = Vector{Tuple{Float64,Float64}}[]
    for level in 1:maxlevel
        h = h0/2^level
        table = generate_table(h, 2)
        push!(tables, table)
    end
    return QuadTS{maxlevel}(h0, origin, table0, Tuple(tables))
end

function integrate(f::Function, q::QuadTS{N}; atol::Real=0.,
    rtol::Real=atol>0 ? 0. : sqrt(eps(Float64))) where {N}

    eval(t) = f(t[1])*t[2] + f(-t[1])*t[2]
    x0, w0 = q.origin
    Σ = f(x0)*w0 + mapsum(eval, q.table0)
    h0 = q.h0
    I = h0*Σ
    E = 0.
    prevI2 = 0.
    prevI = 0.
    for level in 1:N
        table = q.tables[level]
        Σ += mapsum(eval, table)
        h = h0/2^level
        prevI2 = prevI
        prevI = I[1]
        I = h*Σ
        #if level <= 2
        #    E = 1.
        #else 
        #    E = estimate_error(prevI, I[1], prevI2)
        #end 
        #tol = max(norm(I)*rtol, atol)
        #!(E > tol) && break
        if level <= 2
            continue
        end  
        !(abs(I[1]-prevI) > 1e-12) && break
    end
    return I, E
end

function generate_table(h::Float64, step::Int)
    table = Tuple{Float64,Float64}[]
    k = 1
    while true
        t = k*h
        xk, wk = samplepoint(t)
        1 - xk ≤ eps(Float64) && break
        wk ≤ floatmin(Float64) && break
        push!(table, (xk, wk))
        k += step
    end
    reverse!(table)
    return table
end


function samplepoint(t::Float64)
    sinht = sinh(t)
    ϕ = tanh(sinht*π/2)
    ϕ′ = (cosh(t)*π/2)/cosh(sinht*π/2)^2
    return ϕ, ϕ′
end

function estimate_error(prevI::Float64, I::Float64, prevI2::Float64)
    if I-prevI == 0.
        return 0.
    end 
    d1 = log10(norm(I-prevI))
    d2 = log10(norm(I-prevI2))
    d1byd2 = d1^2/d2
    twod1 = 2*d1
    if twod1 >= d1byd2 
        return 10^twod1 
    else 
        return 10^d1byd2 
    end
end


function tanh_sinh(f::Function, a::Real, b::Real, q::QuadTS{N};
    atol::Real=0.0,
    rtol::Real=atol>0 ? 0. : sqrt(eps(Float64))) where {N}
    if a == b
        return 0.0
    else 
        _a = Float64(a)
        _b = Float64(b)
        s = (_b + _a)/2
        t = (_b - _a)/2
        I, E = integrate(u -> f(s + t*u), q; atol=atol/t, rtol=rtol)
        return I*t
    end
end 
