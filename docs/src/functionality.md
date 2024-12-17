# Functionality

This section shows how all components come together to a running optimzation and serves as a walkthrough.

## Walkthrough 

The backbones of this package are the concrete types of [`Settings`](@ref) i.e. the structs containing 
the settings e.g. number of sites in the composite system `N` or subsystem `N_A`, observables `observables`, reduced density matrix `ρ_A`, etc... 
For each lattice model there is an own subtype e.g. [`Settings_XXZ`](@ref) for the XXZ-model, which can 
be conveniently constructed with a constructor e.g. [`XXZ`](@ref).
Second of all after choosing the lattice model, one needs to choose the variational Ansatz e.g [`H_A_BW`](@ref), which returns
a struct [`H_A_Var`](@ref) containing the Hamiltonian Blocks.
To keep the code as short as possible for the user, this is passed to a function [`initialize`](@ref) which returns 
a struct of type [`Init`](@ref) containing the settings, buffer, integration tables for the double-exponential formula and the variational Ansatz s.t.
everything needed for the optimization is saved in one variable, which can then be passed to the optimizer.

## In short

1. Choose a lattice model and use its corresponding constructor (e.g. [`XXZ`](@ref) for XXZ-model).
2. Choose a variational Ansatz (e.g. [`H_A_BW`](@ref) for an Ansatz corresponding to the Bisognano-Wichmann-theorem).
3. Wrap it all up in the initializer function [`initialize`](@ref).
4. Pass the returned struct form step 3 together with initial parameters to the optimizer and run the optimizer.

## Example 

Finding the optimal parameters for the TFIM with `N=10`, `N_A =5`, `Γ=1` with an maximum integration time of `T_max=1.5` with OBC (default)
with [`H_A_BW`](@ref) as an variational ansatz.

```jldoctest
using QuantVarEntHam

init = initialize(TFIM(10, 5, 1., 1.5), H_A_BW)
g_init = [1.,2.,3.,4.,5.]
optimize_LBFGS(g_init, init, gtol = 1e-16, maxiter = 100, print_result = false, show_trace = false)

# output
Diagonalizing the Hamitlonian via exact diagonalization for constructing the ground state density matrix
g_init: [1.0, 2.0, 3.0, 4.0, 5.0]
N: 10
N_A: 5
T_max: 1.5
([0.7998845505676324, 2.339898908834232, 3.4440561575399675, 3.9916008402215475, 4.146669334126589], 1.0913856563696261e-5)
```
The output 
```shell
Diagonalizing the Hamitlonian via exact diagonalization for constructing the ground state density matrix
```
stems from the function [`get_rhoA`](@ref) giving a hint which method is used to extract the ground state.
The optimizer prints the most important underlying settings once and the  real output 
is 
```shell
([0.7998845505676324, 2.339898908834232, 3.4440561575399675, 3.9916008402215475, 4.146669334126589], 1.0913856563696261e-5)
```
The vector corresponds to the optimal parameters and the second element in the tuple is the cost function value at its minimum