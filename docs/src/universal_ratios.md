# Universal ratios 

The eigenspectrum ``\{ \xi_{\alpha} \}`` ("Entanglement spectrum", short: ES) of the Entanglement Hamiltonian (EH) ``H_{\text{A}}`` is the target quantity of this package. 
However, the [`Commutator`](@ref) and the [`QCFL`](@ref) cannot determine the ES uniquely.
These cost functions determine the EH up to a scale factor and an additive constant ``H_{\text{A}}' = \beta H_{\text{A}} + c``,
s.t. the ES is not uniquely determined aswell ``\xi_{\alpha}' = \beta \xi_{\alpha} + c``.
To compare the spectra, the universal ratios are introduced 
```math
\kappa_{\alpha} = \frac{\xi_{\alpha} - \xi_{\alpha_0}}{\xi_{\alpha_1} - \xi_{\alpha_0}},
```
s.t. the scale factor and additive constant are eleminated by division and subtraction, respectively.
```@docs 
universal_ratios(A::AbstractMatrix; α0::Integer=1, α1::Integer=5)
```
```@docs 
var_exact_universal_ratios(cost::AbstractCostFunction, g::Vector{<:AbstractFloat}, α0::Integer=1, α1::Integer=5)
```