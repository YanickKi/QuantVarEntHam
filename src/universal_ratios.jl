export universal_ratios, var_exact_universal_ratios

"""
    var_exact_universal_ratios(cost::AbstractCostFunction, g::Vector{<:Real}, α0::Integer=1, α1::Integer=5)
    var_exact_universal_ratios(model::AbstractModel, ansatz::AbstractAnsatz, g::Vector{<:Real},  α0::Integer=1, α1::Integer=5)
    var_exact_universal_ratios(ρ_A::AbstractMatrix, H_A_var::AbstractMatrix,  α0::Integer=1, α1::Integer=5)


Return the exact and variational universal ratios.

The exact density matrix `ρ_A` is contained in `cost` and `model`, where as the variational Entanglement Hamiltonian `H_A_var`
is obtained with the `ansatz` (saved in `cost` aswell) with a parameter set `g`. 

Handing over the matrices is allowed, too.

# Example 

Variational and exact univeral ratios after a minimization.

```jldoctest  
julia> using QuantVarEntHam

julia> model = TFIM(8,4,1);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix
Regularizing aka making the exact reduced density matrix positive semidefinite with ϵ_reg=1.0e-8

julia> ansatz = H_A_BW(model);

julia> cost = QCFL(model, ansatz, 1);

julia> g_init = [1,2,3,4];

julia> g_opt, _ = optimize(cost, g_init, print_result = false, show_trace = false);

julia> var, exact = var_exact_universal_ratios(cost, g_opt);

julia> var
16-element Vector{Float64}:
 0.0
 0.14870038217678444
 0.5179643526932122
 0.6666647348699999
 1.0
 1.1487003821767865
 1.5179643526932114
 1.5210033796528426
 1.6666647348699977
 1.66970376182963
 2.0389677323460527
 2.1876681145228387
 2.521003379652842
 2.6697037618296284
 3.0389677323460518
 3.1876681145228383

julia> exact
16-element Vector{Float64}:
 0.0
 0.13108370696369803
 0.49093656075557146
 0.6219973957553208
 1.0
 1.1308577295601454
 1.4896925449737326
 1.6203103569003199
 1.6983480175752588
 1.8293816131481835
 2.190224913729735
 2.3023631380571357
 2.419714640996394
 2.5195012809117405
 2.7216788895026087
 2.778789825378197
```
"""
function var_exact_universal_ratios(
    cost::AbstractCostFunction, g::Vector{<:AbstractFloat}, α0::Integer=1, α1::Integer=5
)
    ρ_A = cost.model.ρ_A
    H_A_var = get_H_A!(cost, g)
    return var_exact_universal_ratios(ρ_A, H_A_var, α0, α1)
end

function var_exact_universal_ratios(
    model::AbstractModel,
    ansatz::AbstractAnsatz,
    g::Vector{<:Real},
    α0::Integer=1,
    α1::Integer=5,
)
    ρ_A = model.ρ_A
    H_A_var = Matrix(mat(H_A(ansatz, g)))
    return var_exact_universal_ratios(ρ_A, H_A_var, α0, α1)
end

function var_exact_universal_ratios(
    ρ_A::AbstractMatrix, H_A_var::AbstractMatrix, α0::Integer=1, α1::Integer=5
)
    H_A_exact = -log(Hermitian(Matrix(ρ_A)))
    return universal_ratios(Matrix(H_A_var), α0, α1), universal_ratios(H_A_exact, α0, α1)
end

"""
    universal_ratios(A::AbstractMatrix, α0::Integer=1, α1::Integer=5)

Return the universal ratios of a matrix `A`. 

# Example 

Universal ratios of the variational Ansatz with some arbitrary parameters. 

```jldoctest
julia> using QuantVarEntHam

julia> model = TFIM(8,4,1);
Diagonalizing the Hamiltonian via exact diagonalization for constructing the ground state density matrix
Regularizing aka making the exact reduced density matrix positive semidefinite with ϵ_reg=1.0e-8

julia> ansatz = H_A_BW(model);

julia> g = [1,2,3,4];

julia> H_A_var = mat(H_A(ansatz, g));

julia> universal_ratios(H_A_var)
16-element Vector{Float64}:
 0.0
 0.1731182188567562
 0.5391313151820285
 0.7122495340387878
 1.0
 1.173118218856759
 1.5391313151820274
 1.6606468448909688
 1.7122495340387864
 1.8337650637477279
 2.199778160072996
 2.372896378929754
 2.660646844890969
 2.833765063747727
 3.199778160072994
 3.3728963789297537
```
"""
function universal_ratios(A::AbstractMatrix, α0::Integer=1, α1::Integer=5)
    #ishermitian(A) || throw(ArgumentError("A must be hermitian!"))
    ξ, _ = eigen(Matrix(A))
    return (ξ .- ξ[α0])/(ξ[α1] - ξ[α0])
end
