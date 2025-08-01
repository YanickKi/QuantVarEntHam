# Models

For each model, the Hamiltonian and an type, where the model settings are saved, are implemented. 

Instantiating a subtype of [`AbstractModel`](@ref) is fairly simple.
As an example, if one wants find the EH of the TFIM with 8 sites, Γ=1
and the provided default values (J=-1, OBC), one can instantiate the [`TFIM`](@ref) object
```jldoctest Models
julia> using QuantVarEntHam

julia> model = TFIM(8,4,1)
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix
TFIM

Spin 1//2
Number of spins in the composite system N =  8
Number of spins in the subsystem N_A = 4
Boundary conditions: open
Global prefactor J = -1
Transverse field strength Γ = 1

```
This object can be passed to the cost functions together with the variational Ansatz.

If you want to obtain the corresponding Hamiltonian of the composite system, you can 
use  [`H_TFIM`](@ref)
```jldoctest Models
julia> H = H_TFIM(8,1)
Block
Spin 1//2
Number of spins: 8

-Z₁⊗ Z₂ - Z₂⊗ Z₃ - Z₃⊗ Z₄ - Z₄⊗ Z₅ - Z₅⊗ Z₆ - Z₆⊗ Z₇ - Z₇⊗ Z₈ - X₁ - X₂ - X₃ - X₄ - X₅ - X₆ - X₇ - X₈

```


```@docs 
AbstractModel
```
## Transverse field Ising model (TFIM)
```@docs 
H_TFIM(N::Int, Γ::Real; J::Real = -1, periodic::Bool=false, S::Union{Int64, Rational} = 1//2)
```
```@docs
TFIM
```
### Ansätze
The BW Ansatz [`H_A_BW`](@ref) has the blocks 
```math
\hat{h}_i = J \left ( \frac{1}{2} \sum_{j \in \langle j,i \rangle \cap \text{A}} Z_j Z_i - \Gamma X_i \right ) 
```
The BW violating Ansatz [`H_A_BWV`](@ref) has the blocks 
```math
\hat{h}_i \in  \{ J Z_j Z_{j+1} |  1 \leq j < N_\text{A}-1 \} \cup \{ J\Gamma X_j | 1 \leq j < N_\text{A}\}
```
s.t. the BW violating Ansatz reads
 ```math
\hat{H}_\text{A}^\text{BWV} = J \left (\sum_{i=1}^{N_\text{A} -1} J_{i,i+1}Z_i Z_{i+1} + \Gamma \sum_{i=1}^{N_\text{A}} \Gamma_i X_i \right )
```
## XXZ model 
```@docs 
H_XXZ(N::Int, Δ::Real; periodic::Bool=false, J::Real = +1, S::Union{Int64, Rational} = 1//2)
```
```@docs
XXZ
```
### Ansätze
The BW Ansatz [`H_A_BW`](@ref) has the blocks 
```math
\hat{h}_i = \frac{J}{2} \sum_{j \in \langle j,i \rangle \cap \text{A}} ( X_i X_j + Y_i Y_j + J \Delta Z_i Z_j)
```
The BW violating Ansatz [`H_A_BWV`](@ref) has the blocks 
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
### Ansätze
The BW Ansatz [`H_A_BW`](@ref) has the blocks 
```math
\hat{h}_i = J \left ( \frac{J_\text{Heis}}{2} \sum_{j \in \langle j,i \rangle \cap \text{A}} \vec{S}_i \cdot \vec{S}_{j} + B_x X_i + U_{zz} Z_i^2 \right )
```
The BW violating Ansatz [`H_A_BWV`](@ref) has the blocks 
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

# Getter 

```@docs 
getrho_A(model::AbstractModel) 
```