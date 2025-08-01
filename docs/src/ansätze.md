# Ansätze

In General, all variational ansätze are linear combinations of [`Block`](@ref)s.
On the one hand, there is an ansatz following the BW theorem - the BW like Ansatz - [`H_A_BW`](@ref), which 
is a discretized version of the BW theorem. 
On the other hand, allowing multiple independent, parameters per lattice site often times leads to improved convergence, which
is the motivation for the BW violating Ansatz [`H_A_BWV`](@ref).  

In the case of the BW violating Ansatz, it may happen that some blocks are essentially the zero matrix,
which is why the number of blocks also depends on the hamilton parameters.
As an example, consider the TFIM with `N=8` sites in the composite system, `N_A=4` sites in the subsystem and a tranvserse field strength of
`Γ=1`. 
```jldoctest Ansätze 
julia> using QuantVarEntHam

julia> model1 = TFIM(8,4,1);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix

julia> bwv1 = H_A_BWV(model1)
H_A_BWV
Number of blocks: 7

Spin 1//2
Number of spins: 4

Block 1: 
	-X₁

Block 2: 
	-Z₁⊗ Z₂

Block 3: 
	-X₂

Block 4: 
	-Z₂⊗ Z₃

Block 5: 
	-X₃

Block 6: 
	-Z₃⊗ Z₄

Block 7: 
	-X₄
```

Since `Γ` is non-zero, the pauli X terms appear.
However, if one sets `Γ` to zero, the pauli X term does not appear in the composite system hamiltonian.
Including in the variational Ansatz would not make sense, which is why the number of blocks reduces if one sets `Γ=0`
```jldoctest Ansätze  
julia> model0 = TFIM(8,4,0);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix

julia> ansatz0 = H_A_BWV(model0)
H_A_BWV
Number of blocks: 3

Spin 1//2
Number of spins: 4

Block 1: 
	-Z₁⊗ Z₂

Block 2: 
	-Z₂⊗ Z₃

Block 3: 
	-Z₃⊗ Z₄
```
and only the ising coupling part is included in the variational Ansatz.

Including long range corrections can be achived by increasing `r_max`.
For example including next nearest neighbour corrections in the TFIM
```jldoctest Ansätze 
julia> ansatz1_nnn = H_A_BWV(model1, 2)
H_A_BWV
r_max = 2
Number of blocks: 9

Spin 1//2
Number of spins: 4

Block 1: 
	-X₁

Block 2: 
	-Z₁⊗ Z₂

Block 3: 
	-X₂

Block 4: 
	-Z₂⊗ Z₃

Block 5: 
	-X₃

Block 6: 
	-Z₃⊗ Z₄

Block 7: 
	-X₄

Block 8: 
	Z₁⊗ Z₃

Block 9: 
	Z₂⊗ Z₄
```
Note that only the already existing nearest neighbour couplings are extended to the long range corrections.
```@docs
AbstractAnsatz
```
## BW like Ansatz  
```@docs
H_A_BW
```

## BW violating Ansatz 
```@docs
H_A_BWV
```

## Utils 
```@docs
H_A(ansatz::AbstractAnsatz, g::Vector{<:Real}; digits::Int = 16)
```

## Getter

```@docs 
getr_max(ansatz::AbstractAnsatz)
```
```@docs 
getblocks(ansatz::AbstractAnsatz)
```

