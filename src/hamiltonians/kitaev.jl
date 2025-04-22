using Yao

@with_kw mutable struct Settings_kitaev{T<:AbstractBlock, S<:AbstractMatrix} <:Settings{T, S}
    N::Int = 16
    N_A::Int = 6
    Jz::Float64
    Jx::Float64
    Jy::Float64
    T_max::Float64
    ρ_A::DensityMatrix{2, ComplexF64, Matrix{ComplexF64}} 
    observables::Vector{T}
    meas0::Vector{Float64}  = [expect(observables[i], ρ_A) for i in eachindex(observables)]
    mtrxObs::Vector{S}
end

function kitaev(Jz::Real, Jx::Real, Jy::Real,T_max::Real; ρ_A::DensityMatrix{2}=get_rhoA(H_kitaev(Jz, Jx, Jy),  [2,5,6,9,10,14], 16),
    observables::Vector{<:AbstractBlock}=vcat([repeat(6, Z, eachindex([2,5,6,9,10, 14]))] , [chain(6, put(qbit => Z), put(qbit => Z)) for qbit in 1:6]))

    mtrxObs = mat.(observables)

    return Settings_kitaev{eltype(observables), eltype(mtrxObs)}(
        Jz = Float64(Jz), Jx = Float64(Jx), Jy = Float64(Jy), T_max = T_max,
        ρ_A = ρ_A, observables = observables,
        mtrxObs = mtrxObs,
    ) 
end 


function H_kitaev(Jz, Jx, Jy)
 
    zterm = sum([
        repeat(16, Z, (1,13)),
        repeat(16, Z, (2,14)),
        repeat(16, Z, (3,15)),
        repeat(16, Z, (4,16)),
        repeat(16, Z, (5,9)),
        repeat(16, Z, (6,10)),
        repeat(16, Z, (7,11)),
        repeat(16, Z, (8,12)),
    ])

    xterm = sum([
        repeat(16, X, (1,8)),
        repeat(16, X, (2,5)),
        repeat(16, X, (3,6)),
        repeat(16, X, (4,7)),
        repeat(16, X, (9,13)),
        repeat(16, X, (10,14)),
        repeat(16, X, (11,15)),
        repeat(16, X, (12,16)),
        ])

    yterm = sum([
        repeat(16, Y, (1,5)),
        repeat(16, Y, (2,6)),
        repeat(16, Y, (3,7)),
        repeat(16, Y, (4,8)),
        repeat(16, Y, (9,14)),
        repeat(16, Y, (10,15)),
        repeat(16, Y, (11,16)),
        repeat(16, Y, (12,13)),
        ])
    return -(Jz*zterm + Jx*xterm + Jy*yterm)

end 