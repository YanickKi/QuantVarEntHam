# Models

```@docs 
AbstractModel
```
## TFIM
```@docs 
H_TFIM(N::Int, Γ::Real; J::Real = -1, periodic::Bool=false, S::Union{Int64, Rational} = 1//2)
```
```@docs
TFIM
```
```@docs
TFIM(N::Int, N_A::Int, Γ::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool=false,
J::Real=-1, ρ_A::AbstractMatrix=get_ρ_A(H_TFIM(N, Γ, periodic = periodic, J=J, S=S),  N-N_A+1:N, N))
```
### Ansätze
The BW Ansatz has the blocks 
```math
\hat{h}_i = J \left ( \frac{1}{2} \sum_{j \in \langle j,i \rangle \cap \text{A}} Z_j Z_i - \Gamma X_i \right ) 
```
The BW violating Ansatz has the blocks 
```math
\hat{h}_i \in  \{ J Z_j Z_{j+1} |  1 \leq j < N_\text{A}-1 \} \cup \{ J\Gamma X_j | 1 \leq j < N_\text{A}\}
```
s.t. the BW violating Ansatz reads
 ```math
\hat{H}_\text{A}^\text{BWV} = J \left (\sum_{i=1}^{N_\text{A} -1} J_{i,i+1}Z_i Z_{i+1} - \Gamma \sum_{i=1}^{N_\text{A}} \Gamma_i X_i \right )
```
## XXZ model 
```@docs 
H_XXZ(N::Int, Δ::Real; periodic::Bool=false, J::Real = +1, S::Union{Int64, Rational} = 1//2)
```
```@docs
XXZ
```
```@docs
XXZ(N::Int, N_A::Int, Δ::Real; S::Union{Int64, Rational} = 1//2, r_max::Int=1, periodic::Bool = false,
J::Real=+1, ρ_A::Matrix{ComplexF64}=get_ρ_A(H_XXZ(N, Δ, periodic=periodic, J=J, S = S),  N-N_A+1:N, N))
```
### Ansätze
The BW Ansatz has the blocks 
```math
\hat{h}_i = \frac{J}{2} \sum_{j \in \langle j,i \rangle \cap \text{A}} ( X_i X_j + Y_i Y_j + J \Delta Z_i Z_j)
```
The BW violating Ansatz has the blocks 
```math
\hat{h}_i \in \{J X_j X_{j+1}+J Y_j Y_{j+1} | 1 \leq j < N_\text{A}-1\} \cup \{\Delta J Z_j Z_{j+1} | 1 \leq j < N_\text{A}-1\}
```
s.t. the BW violating Ansatz reads
 ```math
\hat{H}_\text{A}^\text{BWV} = J  \sum_{i=1}^{N_\text{A}-1} \left ( J_{i,i+1}^\text{XX} \left (  X_i X_{i+1} + Y_i Y_{i+1}  \right ) +  J_{i,i+1}^\text{Z} \Delta Z_{i}Z_{i+1}\right )
```
## Pollmann model
```@docs 
H_pollmann(N::Int, J_Heis::Real, Bx::Real, Uzz::Real; periodic::Bool=false, J::Real = 1, S::Union{Int64, Rational}=1)
```
```@docs 
Pollmann
```
```@docs
Pollmann(N::Int, N_A::Int, J_Heis::Real, Bx::Real, Uzz::Real; S::Union{Int64, Rational}=1//1, r_max::Int=1, periodic::Bool=false,
J::Real=+1, ρ_A::Matrix{ComplexF64}=get_ρ_A(H_pollmann(N, J_Heis , Bx, Uzz, periodic = periodic, J=J, S = S),  N-N_A+1:N, N, S=S))
```
### Ansätze
The BW Ansatz has the blocks 
```math
\hat{h}_i = J \left ( \frac{J_\text{Heis}}{2} \sum_{j \in \langle j,i \rangle \cap \text{A}} \vec{S}_i \cdot \vec{S}_{j} + B_x X_i + U_{zz} Z_i^2 \right )
```
The BW violating Ansatz has the blocks 
```math
\hat{h}_i \in \{J J_\text{Heis} X_j X_{j+1} | 1 \leq j < N_\text{A}-1\} \cup \{J J_\text{Heis} Y_j Y_{j+1} | 1 \leq j < N_\text{A}-1\} \cup \{J J_\text{Heis} Z_j Z_{j+1} | 1 \leq j < N_\text{A}-1\} \cup \{ J B_x X_j | 1 \leq j < N_\text{A}\} \cup \{ J U_{zz} Z_j^2 | 1 \leq j < N_\text{A}\}
```
s.t. the BW violating Ansatz reads
 ```math
\hat{H}_\text{A}^\text{BWV} = J \left ( J_\text{Heis} \sum_{i=1}^{N_\text{A}-1} ( J_{i,i+1}^x X_i X_{i+1} + 
    J_{i,i+1}^y Y_i Y_{i+1} + J_{i,i+1}^z Z_{i} Z_{i+1} ) 
    + B_x \sum_{i=1}^{N_\text{A}} B_i^x X_i  
     +  U_{zz} \sum_{i=1}^{N_\text{A}} U_i^{zz} Z_i^2 \right ) 
```

# Utils 

```@docs 
rho_A(H::AbstractBlock, A::AbstractVector{Int}, N::Int; S::Union{Rational, Int} = 1//2) 
```