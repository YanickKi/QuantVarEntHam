# Cost

Cost functions, subtypes of [`AbstractCostFunction`](@ref), are objects which can passed to the [`optimize`](@ref) function. 
They contain all cost function specific settings (such as the model, the ansatz, buffers, etc...).
The cost objects are callable, allowing easy syntax to evaluate the cost function at a given parameter set. 
Consider the following example: 

```jlcon 
julia> model = TFIM(8,4,1);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix

julia> ansatz = H_A_BW(model)
H_A_BW
Number of blocks: 4

Spin 1//2
Number of spins: 4

Block 1: 
	-X₁ - 0.5*Z₁⊗ Z₂

Block 2: 
	-X₂ - 0.5*Z₁⊗ Z₂ - 0.5*Z₂⊗ Z₃

Block 3: 
	-X₃ - 0.5*Z₂⊗ Z₃ - 0.5*Z₃⊗ Z₄

Block 4: 
	-X₄ - 0.5*Z₃⊗ Z₄
```
As you can see, the ansatz has four blocks s.t. four parameters are needed to call the object. 

```jlcon
julia> cost = QCFL(model, ansatz, 1)
QCFL

Model: TFIM (S=1//2, N=8, OBC, N_A = 4, J=-1, Γ=1)
Ansatz: H_A_BW
Integration method: Tanh-sinh quadrature (atol=0, rtol=1.4901161193847656e-8, h0=1, maxlevel=12)
T_max=1
```
Here, we set up the [`QCFL`](@ref) for the [`TFIM`](@ref) with `N=8`, `N_A=4`, `Γ=1`, the BW like Ansatz [`H_A_BW`](@ref) and a maximum integration
time of `T_max=1`. 
See [`QCFL`](@ref) for the default values therein.
The most important settings are printed in the REPL besides the observables, since that would blow up the terminal.
In case you want to see which observables are used, use [`print_observables`](@ref).
```jlcon
julia> print_observables(cost)
Number of Pauli strings: 3

Spin 1//2
Number of spins: 4

Pauli string 1: 
	 	Z₁⊗ Z₂

Pauli string 2: 
	 	Z₂⊗ Z₃

Pauli string 3: 
	 	Z₃⊗ Z₄
```
Now to obtain the cost function value, simply call the `cost` with some parameters
```
julia> g = [1,2,3,4];

julia> cost(g)
0.00031157017647431444
```


```@docs 
AbstractCostFunction
```
```@docs 
AbstractFreeCostFunction
```
## Quantum Classical Feedback Loop (QCFL) 

```@docs 
QCFL
```
```@docs
QCFLBuffer
```

## Commutator

```@docs 
Commutator
```
```@docs
CommutatorBuffer
```

## Relative Entropy

```@docs 
RelativeEntropy
```
```@docs 
RelativeEntropyBuffer
```

## Wrapper for fixing parameters 

```@docs
FixedCost
```
```@docs
fill_full_g(fc::FixedCost, g::Vector{<:Real})
```
## Gradient
```@docs
gradient(cost::AbstractCostFunction, g::Vector{<:Real})
```

## Getter
```@docs
getmodel(cost::AbstractCostFunction)
```
```@docs 
getansatz(cost::AbstractCostFunction)
```
```@docs 
getobservables(cost::AbstractCostFunction)
```
## Printing 
```@docs 
print_model(cost::AbstractCostFunction)
```
```@docs 
print_ansatz(cost::AbstractCostFunction)
```
```@docs 
print_observables(cost::AbstractCostFunction)
```

## Defining and passing buffers

This package allows the user to construct buffers. 
The reason behind this is the following: Imagine one wants to check the convergence of the parameters with respect to the maximum integration time 
`T_max` of the [`QCFL`](@ref) (or testing different ansätze, integration schemes, etc...). 
Instantiating a buffer object for each iteration in e.g. a loop, would allocate huge amounts of memory. 
```jl
using QuantVarEntHam

model = TFIM(8,4,1)

ansatz = H_A_BW(model)

max_times = 0.1:0.1:10 

g_init = [1,2,3,4]

for T_max in max_times 
    cost = QCFL(model, ansatz, T_max) # This allotes a new buffer object in each iteration internally
    optimize(cost, g_init)
end
```
One can instantiate the buffer once and reuse it for each cost.
```jl
using QuantVarEntHam

model = TFIM(8,4,1)

ansatz = H_A_BW(model)

max_times = 0.1:0.1:10 

g_init = [1,2,3,4]

buffer = QCFLBuffer(ansatz) # predefine buffer once

for T_max in max_times 
    cost = QCFL(model, ansatz, T_max, buffer=buffer) # reuse buffer in each iteration
    optimize(cost, g_init)
end
``` 

