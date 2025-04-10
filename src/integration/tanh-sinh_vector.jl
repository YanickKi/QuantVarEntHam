function integrate(f::Function, q::QuadTS{N}, I::AbstractVector, Σ::AbstractVector; atol::Real=0.,
    rtol::Real=atol>0 ? 0. : sqrt(eps(Float64))) where {N}
    
    function eval(t, Σ)
        Σ .+= f(t[1]).*t[2]
        Σ .+= f(-t[1]).*t[2]
    end

    function eval_withmax(t, Σ)
        f1 = f(t[1])
        f11 = f1[1]
        Σ .+= f1 .* t[2]

        f2 = f(-t[1])
        f21 = f2[2]
        Σ .+= f2 .* t[2]

        maxf = max(f11, f21)
        
        return maxf*t[2]
    end 
    x0, w0 = q.origin
    h0 = q.h0
    Σ .+= f(x0).*w0 
    mapsum(eval, q.table0, Σ)
    I .*= Σ.*h0
    E = 0.
    prevI2 = 0.
    prevI = I[1]
    for level in 1:N
        if level == 1
            table = q.tables[level]
            h = h0/2^level
            mapsum(eval, table, Σ)
            I .= Σ.*h
            E = 1.
        else 
            table = q.tables[level]
            maxlr = eval_withmax(table[1], Σ)
            maxj = mapsum_withmax(eval_withmax, Σ, @view table[2:end])
            h = h0/2^level
            prevI2 = prevI
            prevI = I[1]
            I .= Σ.*h
            E = estimate_error(prevI, I[1], prevI2, maxj, maxlr)
        end 
        tol = max(abs(I[1])*rtol, atol)
        !(E > tol) && break
    end
end

function tanh_sinh!(f::Function, a::Real, b::Real, q::QuadTS{N}, I::AbstractVector, Σ::AbstractVector;
    atol::Real=0.0,
    rtol::Real=atol>0 ? 0. : sqrt(eps(Float64))) where {N}
    if a == b
         fill!(I, 0.)
    else 
        _a = Float64(a)
        _b = Float64(b)
        s = (_b + _a)/2
        t = (_b - _a)/2
        integrate(u -> f(s + t*u), q, I, Σ; atol=atol/t, rtol=rtol)
        I .*= t 
    end
    fill!(Σ, 0.)
end 
