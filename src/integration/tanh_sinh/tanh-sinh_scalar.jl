function integrate(f::Function, q::QuadTS{N}; atol::Real=0.,
    rtol::Real=atol>0 ? 0. : sqrt(eps(Float64))) where {N}
    
    eval(t) = f(t[1])*t[2] + f(-t[1])*t[2]

    function eval_withmax(t)
        f1 = copy(f(t[1]))
        f2 = copy(f(-t[1]))
        maxf = max(f1[1], f2[1])
        val = f1*t[2]+f2*t[2]
        return val, maxf*t[2]
    end 
    x0, w0 = q.origin
    Σ = f(x0)*w0 + mapsum(eval, q.table0)
    h0 = q.h0
    I = h0*Σ
    E = 0.
    prevI2 = 0.
    prevI = I[1]
    for level in 1:N
        if level == 1
            table = q.tables[level]
            Σ += mapsum(eval, table)
            h = h0/2^level
            #prevI2 = prevI
            #prevI = I[1]
            I = h*Σ            
            E = 1.
        else 
            table = q.tables[level]
            frfl, maxlr = eval_withmax(table[1])
            Σ_level, maxj = mapsum_withmax(eval_withmax, @view table[2:end])
            Σ += Σ_level+frfl
            h = h0/2^level
            prevI2 = prevI
            prevI = I[1]
            I = h*Σ
            E = estimate_error(prevI, I[1], prevI2, maxj, maxlr)
        end 
        tol = max(abs(I[1])*rtol, atol)
        !(E > tol) && break
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

function estimate_error(prevI::Float64, I::Float64, prevI2::Float64, maxj::Float64, maxlr::Float64)
    if I-prevI == 0.
        return 0.
    end 
    d1 = abs(I-prevI)
    d2 = abs(I-prevI2)
    d3 = 1e-15*maxj
    d4 = maxlr
    #println("d1^2/d2: ", (d1^2/d2))
    #println("2d1: ", (2*d1))
    #println("d3: ", d3)
    #println("d4: ", d4)
    return max(d1^2/d2, 2*d1, d3, d4)
    #d = max(d1^2/d2, 2*d1, d3, d4)
    #return 10^d
end


function (th_s::TanhSinh_scalar)(f::Function, b::Real)
    q = th_s.integration_table
    atol = th_s.atol
    rtol = th_s.rtol
    if iszero(b)
        return zero(b)
    else 
        _a = 0.
        _b = Float64(b)
        s = (_b + _a)/2
        t = (_b - _a)/2
        I, E = integrate(u -> f(s + t*u), q; atol=atol/t, rtol=rtol)
        return I*t
    end
end 
