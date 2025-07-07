function make_integrator(blocks::Vector{<:AbstractMatrix}, ts::TanhSinh)
    q = integration_tables(maxlevel = ts.maxlevel, h0 = ts.h0)
    buffer = zeros(length(blocks)+1)
    scalar_integrate = TanhSinhScalar(q, ts.atol, ts.rtol)
    vector_integrate = TanhSinhVector(q, ts.atol, ts.rtol, buffer)

    return Integrator(scalar_integrate, vector_integrate)
end         


function make_integrator(::Vector{<:AbstractMatrix}, mp::MidPoint)
    return Integrator(MidPointScalar(mp.dt), MidPointVector(mp.dt)) 
end 