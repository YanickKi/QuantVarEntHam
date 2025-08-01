# Cost

Cost functions, subtypes of [`AbstractCostFunction`](@ref), are objects which can passed to the [`optimize`](@ref) function. 
They contain all cost function specific settings (such as the model, the ansatz, buffers, etc...).
The cost objects are callable, allowing easy syntax to evaluate the cost function at a given parameter set. 
Consider the following example: 

```jldoctest Cost 
julia> using QuantVarEntHam

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

```jldoctest Cost
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
```jldoctest Cost
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
```jldoctest Cost
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
## Integration 

Integration is required for the [`QCFL`](@ref).
Often times, the results depend on the accuracy of the integration, which is why different integration methods are offered. 
Constructing and handing over integration methods and its corresponding settings is fairly simple.

Consider an example where one wants to use the [`MidPoint`](@ref) rule with `dt=1e-2`.

```jldoctest Integration 
julia> using QuantVarEntHam

julia> model = TFIM(8,4,1);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix

julia> ansatz = H_A_BW(model);

julia> integrator = MidPoint(1e-2)
MidPoint(0.01)

julia> cost_mp = QCFL(model, ansatz, 1, integrator = integrator)
QCFL

Model: TFIM (S=1//2, N=8, OBC, N_A = 4, J=-1, Γ=1)
Ansatz: H_A_BW
Integration method: Midpoint method (dt=0.01)
T_max=1
```

If no integrator is provided, it will be automatically set to the [`TanhSinh`](@ref) quadrature, with its recommended default values.

```jldoctest Integration  
julia> cost_ts = QCFL(model, ansatz, 1)
QCFL

Model: TFIM (S=1//2, N=8, OBC, N_A = 4, J=-1, Γ=1)
Ansatz: H_A_BW
Integration method: Tanh-sinh quadrature (atol=0, rtol=1.4901161193847656e-8, h0=1, maxlevel=12)
T_max=1
```

As mentioned, The optimal parameters will often times differ for different integrator  
```jldoctest Integration 
julia> g_init = [1,2,3,4];

julia> optimize(cost_mp, g_init, print_result=false, show_trace=false)
([1.5442375810179547, 4.346266068051659, 5.8102829999823715, 6.145496386707021], 2.2719995750666052e-5)

julia> optimize(cost_ts, g_init, print_result=false, show_trace=false)
([1.5440537768894778, 4.345710296202192, 5.80964334649578, 6.14473725960163], 2.2723132440614417e-5)
```

In this case, the differences are small, since the `dt` is chosen small for the mid point rule.
However, varying integration settings can lead to significantly different results.

```@docs
AbstractIntegrator
```

### Tanh-sinh quadrature 

```@docs
TanhSinh
```

### Midpoint rule 
```@docs
MidPoint
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