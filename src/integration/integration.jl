include("tanh_sinh.jl")
include("midpoint.jl")

function tanh_sinh(; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Integer=12, h0::Float64=1.0)

    q = integration_tables(maxlevel = maxlevel, h0 = h0)

    integrators = (
        scalar = (f, T_max) -> tanh_sinh(f, 0, T_max, q, atol=atol, rtol=rtol),
        vector = (I, f, T_max, Σ) -> tanh_sinh!(I, f, 0, T_max, q, Σ, atol=atol, rtol=rtol),
        need_buffer = true
    )
    return integrators
end


function midpoint(; dt::Real = 1e-2)

    integrators = (
        scalar = (f, T_max) -> midpoint(f, 0, T_max, dt)  ,
        vector = (I, f, T_max) -> midpoint!(I, f, 0, T_max, dt),
        need_buffer = false
    )

    return integrators
end 
 