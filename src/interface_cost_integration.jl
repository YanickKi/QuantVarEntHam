function make_integrator(c::AbstractModel, ts::TanhSinh)
    q = integration_tables(maxlevel = ts.maxlevel, h0 = ts.h0)
    buffer = zeros(length(c.blocks+1))
    scalar_integrate = TanhSinh_scalar(q, ts.atol, ts.rtol)
    vector_integrate = TanhSinh_vector(q, ts.atol, ts.rtol, buffer)

    return Integrator(scalar_integrate, vector_integrate)
end         


function make_integrator(::AbstractModel, mp::MidPoint)
    return Integrator(MidPoint_scalar(mp.dt), MidPoint_vector(mp.dt)) 
end 