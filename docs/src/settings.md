# Settings

```@docs 
Settings{M<:AbstractMatrix}
```
```@docs
Settings_TFIM
```
```@docs
TFIM(N::Int, N_A::Int, Γ::Real, T_max::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool=false,
J::Real=-1, ρ_A::AbstractMatrix=get_rhoA(H_TFIM(N, Γ, periodic = periodic, J=J, S=S),  N-N_A+1:N, N),
observables::Vector{<:AbstractMatrix}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1])
```
```@docs
Settings_XXZ
```
```@docs
XXZ(N::Int, N_A::Int, Δ::Real, T_max::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool = false,
J::Real=+1, ρ_A::Matrix{ComplexF64}=get_rhoA(H_XXZ(N, Δ, periodic=periodic, J=J, S = S),  N-N_A+1:N, N),
observables::Vector{<:AbstractMatrix}=[repeat(N_A, Z, (i,i+1)) for i in 1:N_A-1])
```
```@docs 
Settings_pollmann
```
```@docs
pollmann(N::Int, N_A::Int, J_Heis::Real, Bx::Real, Uzz::Real,T_max::Real; S::Union{Int64, Rational}=1//1, r_max::Int=1, periodic::Bool=false,
J::Real=+1, ρ_A::Matrix{ComplexF64}=get_rhoA(H_pollmann(N, J_Heis , Bx, Uzz, periodic = periodic, J=J, S = S),  N-N_A+1:N, N, S=S),
observables::Vector{<:AbstractMatrix}=[repeat(N_A, Z, (i,i+1), S=S) for i in 1:N_A-1])
```