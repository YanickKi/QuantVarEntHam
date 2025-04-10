using LinearAlgebra

abstract type AbstractIntegrator end

struct TanhSinhScalar <: AbstractIntegrator
    q::QuadTS{12}
    atol::Float64
    rtol::Float64
end

struct TanhSinhVector <: AbstractIntegrator
    q::QuadTS{12}
    atol::Float64
    rtol::Float64
end

function integrate(integrator::TanhSinhScalar, f, T_max)
    tanh_sinh(f, 0, T_max, integrator.q, atol=integrator.atol, rtol=integrator.rtol)
end

function integrate(integrator::TanhSinhVector, f, T_max, I, Σ)
    tanh_sinh!(f, 0, T_max, integrator.q, I, Σ, atol=integrator.atol, rtol=integrator.rtol)
end