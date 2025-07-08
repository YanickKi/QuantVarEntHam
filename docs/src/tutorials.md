# Quick Start 

Let's say you want to learn the Entanglement Hamiltonian (EH) of the [`TFIM`](@ref) with 
a transversal field strength of `Γ = 1`, `N=8` spins in the composite system, open boundary conditions (default) and `N_A = 4` sites in the subsystem.
So first of all let's define the model
```jl
using QuantVarEntHam

model = TFIM(8, 4, 1)
```
All model spefic settings and the reduced density matrix are saved in `model`.
Once defining the model, you need to obtain the blocks for the variational Ansatz, which you want to you use.
Let's decide for the Ansatz, which follows the Bisognano Wichmman theorem (BW-theorem, see [`H_A_BW`](@ref))
```jl
blocks = H_A_BW(model)
```
After defining the model and the variational Ansatz you want to learn the EH with, choose the cost function.
This package spefically focuses on the Quantum Classical Feedback Loop ([`QCFL`](@ref)).
So let's define an object, which will represent the cost function with the recommended default values.
```jl
T_max = 1
cost = QCFL(model, blocks, T_max)
```
The (default) observables $\{ Z_i Z_{i+1} | 1 \leq i < N_\text{A} \}$ are monitored up to a time of `T_max = 1`
and with the Tanh-sinh quadrature as an (default) integration method (see [`TanhSinh`](@ref)).
All model and cost function specific settings are set and now just simply pass the 
`cost` object to the optimizer (see [`optimize`](@ref)) together with some initial parameters.
Since we chose the BW Ansatz, we have one block for each lattice site in the subsystem (here: four sites) s.t. we need four initial parameters
```jlcon
julia> g_init = [1,2,3,4];

julia> g_opt, c = optimize(cost, g_init, print_result=false, show_trace = false)
([1.5440537768894713, 4.345710296202176, 5.809643346495764, 6.1447372596016185], 2.27231324406153e-5)
```
A tuple is returned with an vector containing the optimal parameters and the minimum of the cost function.
Note that `print_result=false` and `show_trace=false` avoids that the result and the trace of the minimization (info about iterations, etc...) is printed.
If you want to get the cost function value for a given parameter set you can simply call the cost function object
```jlcon
julia> cost(g_opt)
2.27231324406153e-5
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
([3.3723098537936185, 3.978821738787843], 0.0001512490878026331)
```
Note that only two initial parameters must be provided since we fixed two out of four parameters.
To obtain the full parameter set you can use [`fill_full_g`](@ref)
```jlcon
julia> full_g = fill_full_g(fixed_cost, g_opt_fixed)
4-element Vector{Float64}:
 1.5
 3.3723098537936185
 3.2
 3.978821738787843
```