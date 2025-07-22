# Quick Start 

Let's say you want to learn the Entanglement Hamiltonian (EH) of the [`TFIM`](@ref) with 
a transverse field strength of `Γ = 1`, `N=8` spins in the composite system, open boundary conditions (default) and `N_A = 4` sites in the subsystem.
So first of all let's define the model
```jl
using QuantVarEntHam

model = TFIM(8, 4, 1)
```
All model spefic settings and the reduced density matrix are saved in `model`.
Once defining the model, you need to obtain the variational Ansatz, which you want to you use.
Let's decide for the Ansatz, which follows the Bisognano Wichmman theorem (BW-theorem, see [`H_A_BW`](@ref))
```jl
ansatz = H_A_BW(model)
```
After defining the model and the variational Ansatz you want to learn the EH with, choose the cost function.
This package spefically focuses on the Quantum Classical Feedback Loop ([`QCFL`](@ref)).
So let's define an object, which will represent the cost function with the recommended default values.
```jl
T_max = 1
cost = QCFL(model, ansatz, T_max)
```
The (default) observables $\{ Z_i Z_{i+1} | 1 \leq i < N_\text{A} \}$ are monitored up to a time of `T_max = 1`
and with the Tanh-sinh quadrature as an (default) integration method (see [`TanhSinh`](@ref)).
All model and cost function specific settings are set and now just simply pass the 
`cost` object to the optimizer (see [`optimize`](@ref)) together with some initial parameters.
Since we chose the BW Ansatz, we have one block for each lattice site in the subsystem (here: four sites) s.t. we need four initial parameters
```jlcon
julia> g_init = [1,2,3,4];

julia> g_opt, c = optimize(cost, g_init, print_result=false, show_trace = false)
([1.5440537768894778, 4.345710296202192, 5.80964334649578, 6.14473725960163], 2.2723132440614417e-5)

```
A tuple is returned with an vector containing the optimal parameters and the minimum of the cost function.
Note that `print_result=false` and `show_trace=false` avoids that the result and the trace of the minimization (info about iterations, etc...) is printed.
You can obtain the Entanglement Hamiltonian via [`H_A`](@ref)
```jlcon
julia> H_A_opt = H_A(ansatz, g_opt, digits = 2)
Block
Spin 1//2
Number of spins: 4

-1.54*X₁ - 2.945*Z₁⊗ Z₂ - 4.35*X₂ - 5.08*Z₂⊗ Z₃ - 5.81*X₃ - 5.975*Z₃⊗ Z₄ - 6.14*X₄
```
where we set the number of digits in the coefficients to `digits=2` for readability.
If you want to get the cost function value for a given parameter set you can simply call the cost function object
```jlcon
julia> cost(g_opt)
2.2723132440614417e-5
```

# Fixing parameters
Let's say you want to fix some parameters during the minimization procedure.
For this purpose you can use the [`FixedCost`](@ref) wrapper.
Consider the example above with the already defined cost function `cost`.
Let's say you want to fix the first and third parameter to a value of `1.5` and `3.2` respectively.
```jl
fixed_indices = [1,3]
fixed_values = [1.5,3.2]
``` 
Now wrap the existing cost function `cost` into a `FixedCost` with the given indices and values 
```jl
fixed_cost = FixedCost(cost, fixed_indices, fixed_values)
```
and run the optimization with some initial parameters
```jlcon
julia> g_init_fixed = [2,4];

julia> g_opt_fixed, c_fixed = optimize(fixed_cost, g_init_fixed, print_result=false, show_trace = false)
([3.372309853793615, 3.9788217387878415], 0.0001512490878026316)
```
Note that only two initial parameters must be provided since we fixed two out of four parameters.
To obtain the full parameter set you can use [`fill_full_g`](@ref)
```jlcon
julia> full_g = fill_full_g(fixed_cost, g_opt_fixed)
4-element Vector{Float64}:
 1.5
 3.372309853793615
 3.2
 3.9788217387878415
```