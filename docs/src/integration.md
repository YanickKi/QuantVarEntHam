# Integration 

Integration is required for the [`QCFL`](@ref).
Often times, the results depend on the accuracy of the integration, which is why different integration methods are offered. 
Constructing and handing over integration methods and its corresponding settings is fairly simple.

Consider an example where one wants to use the [`MidPoint`](@ref) rule with `dt=1e-2`.

```jlcon 
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

```jlcon 
julia> cost_ts = QCFL(model, ansatz, 1)
QCFL

Model: TFIM (S=1//2, N=8, OBC, N_A = 4, J=-1, Γ=1)
Ansatz: H_A_BW
Integration method: Tanh-sinh quadrature (atol=0, rtol=1.4901161193847656e-8, h0=1, maxlevel=12)
T_max=1
```

As mentioned, The optimal parameters will often times differ for different integrator  
```jlcon
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

## Tanh-sinh quadrature 

```@docs
TanhSinh
```

## Midpoint rule 
```@docs
MidPoint
```