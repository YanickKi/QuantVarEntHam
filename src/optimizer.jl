using Optim

export optimize

"""
    optimize(cost::AbstractCostFunction, g_init::Vector{<:Real}; ∇_tol::Real=1e-16, maxiter::Integer=1000, show_trace::Bool=true, print_result::Bool= true)

Minimize cost function using the [LBFGS-optimizer from `Optim.jl`](https://julianlsolvers.github.io/Optim.jl/stable/algo/lbfgs/) and return 
the optimal parameters and the cost function value at the minimum.

# Required arguments
- `cost`: cost function.
- `g_init`: initial parameters.

# Keyword arguments 
- `∇_tol`: maximum infinity norm of the gradient in order to stop the minimization.
- `maxiter`: maximum number of iterations in order to stop minimizing.
- `show_trace`: true for showing the trace of the minimizing procedure, false otherwise.
- `print_result` true to print optimal parameters and result, false otherwise.

# Example 

Finding the Entanglement Hamiltonian of the [`TFIM`](@ref) with `N=8`, OBC, `Γ=1`, `N_A=4` and the BW like Ansatz [`H_A_BW`](@ref), using 
the [`QCFL`](@ref).


```jlcon 
julia> model = TFIM(8,4,1);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix

julia> ansatz = H_A_BW(model);

julia> cost = QCFL(model, ansatz, 1);

julia> g_init = [1,2,3,4];

julia> optimize(cost, g_init);
Iter     Function value   Gradient norm 
     0     3.115702e-04     6.983472e-04
 * time: 0.020919084548950195
     1     8.449873e-05     1.344302e-04
 * time: 0.7815248966217041
     2     2.800500e-05     5.200157e-05
 * time: 0.8070659637451172
     3     2.641797e-05     5.919818e-06
 * time: 0.8340210914611816
     4     2.634565e-05     3.705769e-06
 * time: 0.8612968921661377
     5     2.632947e-05     4.628599e-07
 * time: 0.8869969844818115
     6     2.632558e-05     2.331509e-06
 * time: 0.922123908996582
     7     2.529544e-05     2.372577e-05
 * time: 0.983971118927002
     8     2.527188e-05     2.625788e-05
 * time: 1.0214660167694092
     9     2.346181e-05     7.991659e-06
 * time: 1.1042149066925049
    10     2.311557e-05     7.170567e-06
 * time: 1.1501638889312744
    11     2.302321e-05     6.947690e-06
 * time: 1.1689889430999756
    12     2.272929e-05     1.160222e-06
 * time: 1.1969239711761475
    13     2.272481e-05     4.668009e-07
 * time: 1.2245891094207764
    14     2.272317e-05     1.700687e-07
 * time: 1.2523829936981201
    15     2.272313e-05     1.039140e-08
 * time: 1.2805180549621582
    16     2.272313e-05     6.533077e-11
 * time: 1.3085510730743408
    17     2.272313e-05     1.704691e-14
 * time: 1.326706886291504
    18     2.272313e-05     2.371184e-17
 * time: 1.3455560207366943
 * Status: success

 * Candidate solution
    Final objective value:     2.272313e-05

 * Found with
    Algorithm:     L-BFGS

 * Convergence measures
    |x - x'|               = 4.26e-10 ≰ 0.0e+00
    |x - x'|/|x'|          = 6.94e-11 ≰ 0.0e+00
    |f(x) - f(x')|         = 1.56e-19 ≰ 0.0e+00
    |f(x) - f(x')|/|f(x')| = 6.86e-15 ≰ 0.0e+00
    |g(x)|                 = 2.37e-17 ≤ 1.0e-16

 * Work counters
    Seconds run:   1  (vs limit Inf)
    Iterations:    18
    f(x) calls:    66
    ∇f(x) calls:   66

[1.5440537768894778, 4.345710296202192, 5.80964334649578, 6.14473725960163]

```

"""
function optimize(
    cost::AbstractCostFunction,
    g_init::Vector{<:Real};
    ∇_tol::Real=1e-16,
    maxiter::Integer=1000,
    show_trace::Bool=true,
    print_result::Bool=true,
)
    free_indices = get_free_indices(cost)

    @assert length(g_init) == length(free_indices) "You provided $(length(g_init)) initial parameters but there are $(length(free_indices)) free parameters!"

    result = Optim.optimize(
        Optim.only_fg!((F, G, g) -> fg!(F, G, cost, g)),
        Float64.(g_init),
        LBFGS(),
        Optim.Options(;
            g_tol=∇_tol,
            store_trace=false,
            show_trace=show_trace,
            show_warnings=true,
            iterations=maxiter,
        ),
    )
    g_opt = Optim.minimizer(result)
    if print_result == true
        println(result)
        println(g_opt)
    end
    return g_opt, cost(g_opt)
end
