# Spinoperators

Pauli strings and their linear combination (blocks) are implemented s.t the user 
deals with these objects instead of the matrices.
This allows nice printing and offers an easy overview over the ansätze, hamiltonians and observables. 

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


```jlcon 
julia> ps1 = PauliString(4,"Z", (1,2));

julia> ps2 = PauliString(4,"X", (1,2));

julia> ps3 = PauliString(4,"Y", (1,2));

julia> ps1+ps2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + X₁⊗ X₂

julia> block1 = ps1+ps2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + X₁⊗ X₂

julia> block1 + ps3
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + X₁⊗ X₂ + Y₁⊗ Y₂

julia> block2 = ps1-ps2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ - X₁⊗ X₂

julia> block1-block2
Block
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂ + X₁⊗ X₂ - Z₁⊗ Z₂ + X₁⊗ X₂

```

Multiplication by scalars of the a pauli string or block is defined aswell
```jlcon
julia> 2*ps1
Block
Spin 1//2
Number of spins: 4

2.0Z₁⊗ Z₂

julia> 2*block1
Block
Spin 1//2
Number of spins: 4

2.0Z₁⊗ Z₂ + 2*X₁⊗ X₂

```

Multiplication of pauli strings or blocks is not implemented, since the current models do not need it. 
However, taking the power of a pauli string is allowed 
```jlcon 
julia> ps1^3
Pauli string
Spin 1//2
Number of spins: 4

(Z₁⊗ Z₂)³

```