function print_integration_name(io::IO, integrator::Integrator{S, V}) where {S<:AbstractScalarIntegrator, V<:AbstractVectorIntegrator}
    print_integration_name(io, integrator.scalar_integrate)
end     

function print_integration_name(io::IO, ::TanhSinhScalar)
    print(io::IO, "Tanh-sinh quadrature")   
end 

function print_integration_name(io::IO, ::MidPointScalar)
    print(io::IO, "Midpoint method")   
end 

function print_integration_settings_short(io::IO, scalar_integrate::TanhSinhScalar{N}) where {N}
    h0   = float_to_int(scalar_integrate.integration_table.h0)
    atol = float_to_int(scalar_integrate.atol)
    rtol = float_to_int(scalar_integrate.rtol)
    print(io, " (atol=", atol, ", rtol=", rtol, ", h0=", h0 , ", maxlevel=", N, ")")
end 

function print_integration_settings_short(io::IO, scalar_integrate::MidPointScalar)
    print(io, " (dt=", scalar_integrate.dt, ")")
end 
