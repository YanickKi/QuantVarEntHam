function integrate!(I::AbstractVector, f::Function, q::QuadTS{N}, Σ::AbstractVector; atol::Real=0.,
    rtol::Real=atol>0 ? 0. : sqrt(eps(Float64))) where {N}
    
    function eval!(Σ, t)
        Σ .+= f(t[1]).*t[2]
        Σ .+= f(-t[1]).*t[2]
    end

    function eval_withmax!(Σ, t)
        f1 = f(t[1])
        f11 = f1[1]
        Σ .+= f1 .* t[2]

        f2 = f(-t[1])
        f21 = f2[2]
        Σ .+= f2 .* t[2]
(1e-2)
        maxf = max(f11, f21)
        
        return maxf*t[2]
    end 
    x0, w0 = q.origin
    h0 = q.h0
    Σ .+= f(x0).*w0 
    mapsum!(Σ, eval!, q.table0)
    I .*= Σ.*h0
    E = 0.
    prevI2 = 0.
    prevI = I[1]
    for level in 1:N
        if level == 1
            table = q.tables[level]
            h = h0/2^level
            mapsum!(Σ, eval!, table)
            I .= Σ.*h
            E = 1.
        else 
            table = q.tables[level]
            maxlr = eval_withmax!(Σ, table[1])
            maxj = mapsum_withmax!(Σ, eval_withmax!, @view table[2:end])
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

function (th_v::TanhSinh_vector)(I::AbstractVector, f::Function, b::Real)
    q       = th_v.integration_table
    atol    = th_v.atol
    rtol    = th_v.rtol
    Σ       = th_v.buffer
    if iszero(b)
         fill!(I, 0.)
    else 
        _a = 0.
        _b = Float64(b)
        s = (_b + _a)/2
        t = (_b - _a)/2
        integrate!(I, u -> f(s + t*u), q, Σ; atol=atol/t, rtol=rtol)
        I .*= t 
    end
    fill!(Σ, 0.)
end 
