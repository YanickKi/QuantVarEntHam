# Settings

```@docs 
Settings{T<:AbstractBlock,S<:AbstractMatrix} 
```
```@docs
Settings_XXZ
```
```@docs
XXZ(N::Int, N_A::Int, Δ::Real, T_max::Real; r_max::Int=1, periodic::Bool = false, atol::Real=0.0, rtol::Real=atol>0 ? 0. : sqrt(eps(Float64)),
    signHam::Integer=+1, ρ_A::DensityMatrix{2}=get_rhoA(H_XXZ(N, Δ, periodic=periodic, signHam=signHam),  N-N_A+1:N, N),
    observables::Vector{<:AbstractBlock}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1], dt::Float64 = 0.01)
```
```@docs
Settings_TFIM
```
```@docs
TFIM(N::Int, N_A::Int, Γ::Real, T_max::Real; r_max::Int=1, periodic::Bool=false, atol::Real=0.0, rtol::Real=atol>0 ? 0. : sqrt(eps(Float64)),
    signHam::Integer=-1, ρ_A::DensityMatrix{2}=get_rhoA(H_TFIM(N, Γ, periodic = periodic, signHam=signHam), N-N_A+1:N, N),
    observables::Vector{<:AbstractBlock}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1], dt::Real=0.01) 
```