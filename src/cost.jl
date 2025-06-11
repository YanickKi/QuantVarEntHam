using LinearAlgebra
using QuadGK


function get_H_A!(g::Vector{<:AbstractFloat}, init::Init)
    init.buff.H_A .= g[1].*init.blocks[1]

    @fastmath @inbounds @simd for i in 2:length(g)
        init.buff.H_A .+= g[i].*init.blocks[i]
    end 
end 
```
    function QCFL(;integration_method = tanh_sinh(), g1::Real=NaN)

Constructs the cost function for the QCFL with the given integration method and fixed first parameter (if value provided)

# Keyowrd arguments 
- `integration_method=tanh_sinh()`: integration scheme to evaluate the cost function 
- `g1::Real=NaN`: fix first paramter to given value if proved, will be optimized otherwise

# Returns 
- `(F, G, g, init) -> fg!(F, G, g, init, integration_method, fix_first_index = !isnan(g1))`: computes cost and gradient provided `F` and `G` are not `nothing` and returns the cost.
    - `F`: cost function is not evaluated if `F = nothing`, is evaluated otherwise
    - `G`: vector to save the gradient in, if `G=nothing` gradient is not computed 
    - `g`: parameters to be optimized
    - `init`: see [`Init`](@ref)
- `g1`: value of first parameter if provided, NaN otherwise

This function allows one to conveniently hand the cost function to be minimized to [optimize_LBFGS](@ref).
```
function QCFL(;integration_method = tanh_sinh(), g1::Real=NaN)
    return (F, G, g, init) -> fg!(F, G, g, init, integration_method, fix_first_index = !isnan(g1)), Float64(g1)
end 



function fg!(F, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, init::Init, integration_method; fix_first_index::Bool = false)
    @unpack observables, T_max = init.set
    get_H_A!(g, init)
    init.buff.H_A .*= -1im
    method_scalar, method_vector, need_buffer = integration_method
    if G !== nothing
        if need_buffer == true
            method_vector(init.buff.C_G_result, t -> integrand_gradient(t, init, fix_first_index), T_max, init.buff.Σ)
        else 
            method_vector(init.buff.C_G_result, t -> integrand_gradient(t, init, fix_first_index), T_max)
        end 
        init.buff.C_G_result ./= length(observables)*T_max
        G[:] .= @view init.buff.C_G_result[2:end]
    end
    if F !== nothing 
        if G === nothing
            In = method_scalar(t -> integrand_cost(t, init), T_max)/(length(observables)*T_max)
            return In
        else 
            return init.buff.C_G_result[1]
        end 
    end
end 


function integrand_gradient(t::AbstractFloat, init::Init, fix_first_index::Bool)
    @unpack ρ_A, meas0, observables = init.set
    buff = init.buff
    
    buff.H_A_forexp .=  t .* buff.H_A 
    pullback =  own_rrule(exp, buff.H_A_forexp, init.exp_buf)
    
    evolve(init)  
    
    δ = dev(init, 1)
    buff.sumobs .= δ.*observables[1]
    buff.C_G[1] = δ^2

    @fastmath @inbounds @simd for i in 2:lastindex(observables)
        δ = dev(init, i)
        buff.sumobs .+= δ.*observables[i]
        buff.C_G[1] += δ^2
    end 

    pullback(mul!(buff.dAforpb, buff.ρ_A_right, buff.sumobs))

    @fastmath @inbounds @simd for i in 1+fix_first_index:lastindex(buff.C_G)-1+fix_first_index
        buff.C_G[i+1-fix_first_index] = gradient_component(init, i, t)
    end    
    
    return buff.C_G
end

function integrand_cost(t::AbstractFloat,  init::Init)
    @unpack ρ_A, meas0, observables = init.set 

    buff = init.buff

    buff.H_A_forexp .=  t .* buff.H_A 
    exp_bufered!(buff.H_A_forexp, init.exp_buf)
    evolve(init)

    c = 0.
    @fastmath @inbounds @simd for i in eachindex(observables)
        c += dev(init, i)^2
    end 
    return c
end

function evolve(init::Init)
    @unpack ρ_A = init.set
    buff = init.buff
    U = init.exp_buf.X
    mul!(buff.ρ_A_evolved, U, mul!(buff.ρ_A_right, ρ_A, U'))
end 

function dev(init::Init, index_observable::Integer)
    @unpack meas0, observables = init.set
    @inbounds mul!(init.buff.evolob, observables[index_observable], init.buff.ρ_A_evolved)
    @inbounds δ = real(tr(init.buff.evolob))- meas0[index_observable]
    return δ
end 

function gradient_component(init::Init, index_parameter, t::Real)
    buff = init.buff
    frech = init.exp_buf.∂X
    mul!(buff.frechTimesBlock, frech, init.blocks[index_parameter])
    return 4*t*imag(tr(buff.frechTimesBlock))
end


#########################################################################################
#                                                                                       #
#                   Commutator cost function                                            #
#                                                                                       #
#########################################################################################

```
    commutator()

Constructs the commutator as cost function.

# Returns 
- `(F, G, g, init) -> comm_fg!(F, G, g, init)
    - `F`: cost function is not evaluated if `F = nothing`, is evaluated otherwise
    - `G`: vector to save the gradient in, if `G=nothing` gradient is not computed 
    - `g`: parameters to be optimized
    - `init`: see [`Init`](@ref)
- `NaN`: no parameter to be fixed. 

This function allows one to conveniently hand the cost function to be minimized to [optimize_LBFGS](@ref).
```
function commutator()
    return (F, G, g, init) -> comm_fg!(F, G, g, init), NaN 
end 

function comm_fg!(F::Union{AbstractFloat, Nothing}, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, init::Init)
 
    get_H_A!(g, init)
    C = 0.
    if G !== nothing
        C = comm_grad!(G, init)
    end
    if F !== nothing 
        if G === nothing
            C = comm_cost(init)
            return C
        else 
            return C
        end 
    end
end 


function comm_cost(init::Init)

    comm = init.buff.dAforpb
    H_A = init.buff.H_A
    ρ_A = init.set.ρ_A

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm
    
    return norm(comm)/(2*norm(H_A)*norm(ρ_A))
end 

function comm_grad!(G::Vector{<:AbstractFloat}, init::Init)

    comm = init.buff.dAforpb
    H_A = init.buff.H_A
    ρ_A = init.set.ρ_A    
    temp = init.buff.H_A_forexp

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm


    Λ_comm = norm(comm)
    Λ_H_A = norm(H_A)
    Λ_ρ_A = norm(ρ_A)

    @fastmath @inbounds @simd for index in eachindex(G)
        mul!(mul!(temp, ρ_A, init.blocks[index]), init.blocks[index], ρ_A, 1, -1) # compute [h_index, ρ_A] and save it in temp
        ∂Λ_comm = 1/Λ_comm * real(sum(had!(temp, temp, comm)))
        ∂Λ_H_A = 1/Λ_H_A * real(sum(had!(temp, H_A, init.blocks[index])))
        G[index] = 1/(2*Λ_ρ_A*Λ_H_A^2)*(∂Λ_comm * Λ_H_A - Λ_comm * ∂Λ_H_A) 
    end 
    return Λ_comm/(2*Λ_H_A * Λ_ρ_A) 
end 

function had!(buf::AbstractMatrix, A::AbstractMatrix,B::AbstractMatrix)
    n = size(A)[1]
    for j in 1:n
       for i in 1:n
         @inbounds buf[i,j] = A[i,j] *B[i,j]
       end
    end
    return buf
end


```
    relative_entropy()

Constructs the commutator as cost function.

# Returns 
- `(F, G, g, init) -> relative_entropy_fg!
    - `F`: cost function is not evaluated if `F = nothing`, is evaluated otherwise
    - `G`: vector to save the gradient in, if `G=nothing` gradient is not computed 
    - `g`: parameters to be optimized
    - `init`: see [`Init`](@ref)
- `NaN`: no parameter to be fixed. 

This function allows one to conveniently hand the cost function to be minimized to [optimize_LBFGS](@ref).
```

#########################################################################################
#                                                                                       #
#                   Relative entropy                                                    #
#                                                                                       #
#########################################################################################

function relative_entropy()
    return (F, G, g, init) -> relative_entropy_fg!(F, G, g, init), NaN
end 

function relative_entropy_fg!(F::Union{AbstractFloat, Nothing}, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, init::Init)
 
    get_H_A!(g, init)
    C = 0
    if G !== nothing
        C = relative_entropy_grad!(G, length(g), init)
    end
    if F !== nothing 
        if G === nothing
            C = relative_entropy_cost(init)
            return C
        else 
            return C
        end 
    end
end 

function relative_entropy_cost(init::Init)
    @unpack ρ_A  = init.set
    buff = init.buff
    
    buff.H_A_forexp .= -1 .*init.buff.H_A
    exp_bufered!(buff.H_A_forexp, init.exp_buf)
    exp_H_A = init.exp_buf.X

    S1 = tr(mul!(buff.dAforpb, ρ_A, buff.H_A))
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 

function relative_entropy_grad!(G::Union{AbstractVector, Nothing}, num_params::Integer, init::Init)
    @unpack ρ_A  = init.set
    buff = init.buff
    buff.H_A_forexp .= -1 .*buff.H_A

    exp_bufered!(buff.H_A_forexp, init.exp_buf)
    exp_H_A = init.exp_buf.X
    γ = -1/(tr(exp_H_A))

    @fastmath @inbounds @simd for i in 1:num_params
        ∂S1 = tr(mul!(buff.dAforpb, ρ_A, init.blocks[i]))
        ∂S2 = γ * tr(mul!(buff.dAforpb, exp_H_A, init.blocks[i]))
        G[i] = ∂S1+∂S2
    end

    S1 = tr(mul!(buff.dAforpb, ρ_A, buff.H_A))
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 

```
    tanh_sinh_integrand(g::Vector{<:AbstractFloat}, u1::Real, u2::Real, num_points::Int, init::Init)

Return a vector containing the integrand after the  Tanh-sinh transformation at different points.
The integrand is sampled uniformly between `u1` and `u2` `num_points ` times.

# Arguments
- `g::Vector{<:AbstractFloat}`: parameters
- `u1::Real`: lowest sampling point 
- `u2::Real`: highest sampling point
- `num_points::Int`: number of points
- `init`: see [Init](@ref)
```

#########################################################################################
#                                                                                       #
#                   INTEGRAND AFTER TANH-SINH TRANSFORMATION                            #
#                                                                                       #
#########################################################################################

function tanh_sinh_integrand(g::Vector{<:AbstractFloat}, u1::Real, u2::Real, num_points::Int, init::Init)
    @unpack T_max, observables = init.set
    u1_ = Float64(u1)
    u2_ = Float64(u2)
    u = range(start = u1_, stop = u2_, length = num_points)
    sinhu = sinh.(u)
    Φ = tanh.(sinhu .* π/2)
    ϕ′ = (cosh.(u) .*π/2) ./ (cosh.(sinhu .* π/2)).^2
    c = Vector{Float64}(undef, length(u))
    get_H_A!(g, init)
    for i in eachindex(u) 
        c[i] = integrand_onlycost(T_max/2*(Φ[i]+1), init) * ϕ′[i]
    end 

    return c.*T_max/(2*length(observables))
end 


```
    cost_count(g::Vector{<:AbstractFloat}, init::Init; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Int = 12, h0::Real = 1.)

Return the cost function value and the number of samples needed for integration with the Tanh-sinh quadrature as a tuple.

# Required Arguments
- `g::Vector{<:AbstractFloat}`: parameters
- `init`: see [Init](@ref)

# Keyword arguments 
- `atol::Real=0.0`: absolute tolerance for integration 
- `rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64))`: relative tolerance for integration
- `maxlevel::Int=12` maximum number of levels, i.e. halving of the integration step width 
- `h0::Real=1`: initial integration step size  
```
#########################################################################################
#                                                                                       #
#                                COUNT INTEGRAND EVALUTAIONS                            #
#                                                                                       #
#########################################################################################

function cost_count(g::Vector{<:AbstractFloat}, init::Init; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Int = 12, h0::Float64 = 1.)
    @unpack observables, T_max = init.set
    q = integration_tables(maxlevel = maxlevel, h0 = h0)
    count = [0]
    get_H_A!(g, init)
    init.buff.H_A .*= -1im 
    I = tanh_sinh(t -> integrand_onlycost_count(t, init, count),0., T_max, q, atol = atol, rtol = rtol)/(length(observables)*T_max)
    return I, count[1]
end 

function integrand_count(t::AbstractFloat,  init::Init, count::Vector{Int})
    @unpack ρ_A, meas0, observables = init.set
    count[1] += 1
    return integrand_onlycost(t, init)
end

```
    quadgk_count!(g::Vector{<:AbstractFloat}, init::Init; rtol=sqrt(eps), atol=0, maxevals=10^7, order=7)

Return the cost function value and the number of samples needed for integration with the Gauss-Kronrod method as a tuple.

# Required Arguments
- `g::Vector{<:AbstractFloat}`: parameters
- `init`: see [Init](@ref)

# Keyword arguments 
- `atol=0.0`: absolute tolerance for integration 
- `rtol=sqrt(eps)`: relative tolerance for integration
- `maxevals=10^7` maximum number of samples 
- `order`: order for the quadrature  
```
#########################################################################################
#                                                                                       #
#                                QUADGK                                                 #
#                                                                                       #
#########################################################################################

function quadgk_count!(g::Vector{<:AbstractFloat}, init::Init; rtol=sqrt(eps), atol=0, maxevals=10^7, order=7)
    @unpack observables, T_max, atol, rtol = init.set
    get_H_A!(g, init)
    init.buff.H_A .*= -1im
    I, E, count =  quadgk_count(t -> integrand_cost(t, init),0., T_max, rtol = rtol, atol = atol, maxevals=maxevals, order=order)
    return I/(length(observables)*T_max), count 
end