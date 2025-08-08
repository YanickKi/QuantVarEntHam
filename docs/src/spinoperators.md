# Spin Operators

Pauli strings and their linear combination (blocks) are implemented s.t the user 
deals with these objects instead of the matrices.
This allows for nice printing and offers an easy overview over the ansätze, hamiltonians and observables. 

```@docs
AbstractBlock
```
```@docs
PauliString

```
```@docs
Block
```
```@docs
mat(block::AbstractBlock)
```

# Algebra

Adding or subtracting two pauli strings and blocks is possible 


```jldoctest algebra
julia> using QuantVarEntHam

julia> ps1 = z(4, (1,2));

julia> ps2 = x(4, (3,4));

julia> ps3 = z(4,(1,2));

julia> ps1+ps2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + X₃⊗ X₄

julia> block1= ps1+ps2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + X₃⊗ X₄

julia> block1 + ps3
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + X₃⊗ X₄ + Z₁⊗ Z₂

julia> block2 = ps1-ps2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ - X₃⊗ X₄

julia> block1-block2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + X₃⊗ X₄ - Z₁⊗ Z₂ + X₃⊗ X₄

```

Multiplication by scalars of a pauli string or block is defined aswell
```jldoctest algebra
julia> 2*ps1
Block
Spin 1//2
Number of spins: 4

2*Z₁⊗ Z₂

julia> 2*block1
Block
Spin 1//2
Number of spins: 4

2*Z₁⊗ Z₂ + 2*X₃⊗ X₄

```

Multiplication of pauli strings or blocks is not implemented, since the current models do not need it. 
However, taking the power of a pauli string is allowed 
```jldoctest algebra 
julia> ps1^3
Pauli string
Spin 1//2
Number of spins: 4

(Z₁⊗ Z₂)³

```