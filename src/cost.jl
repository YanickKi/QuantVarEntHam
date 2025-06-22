using LinearAlgebra
using QuadGK


function get_H_A!(H_A::AbstractMatrix, g::Vector{<:AbstractFloat}, blocks)
    H_A .= g[1].*blocks[1]

    @fastmath @inbounds @simd for i in 2:length(g)
        H_A .+= g[i].*blocks[i]
    end 
    return H_A 
end 
"""
    QCFL(T_max::Real ;integration_method = tanh_sinh(), g1::Real=NaN)

Constructs the cost function for the QCFL with the given integration method and fixed first parameter (if value provided)

# Required arguments
- `T_max::Real`: maximum integration time

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

This function allows one to conveniently hand the cost function to be minimized to [`optimize_LBFGS`](@ref).
"""
function QCFL(set::Settings, T_max::Real, H_A::Function; 
    observables::Vector{<:AbstractMatrix}=[repeat(set.N_A, Z, (i,i+1), S = set.S) for i in 1:set.N_A-1], integration_method = tanh_sinh(), g1::Real=NaN,
    buffer::Union{Nothing, QCFL_buffer} = nothing)
    @unpack ρ_A = set
    
    meas0 = [expect(observable, ρ_A) for observable in observables]

    blocks = H_A(set)

    buff = something(buffer, make_QCFL_buffer(set, length(blocks)-!isnan(g1), observables))
 
    return (F, G, g) -> fg!(F, G, g, integration_method, T_max, ρ_A, observables, meas0, blocks, buff, fix_first_index = !isnan(g1)), Float64(g1), length(blocks)
end 



function fg!(F, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, integration_method, T_max::Real, ρ_A::AbstractMatrix, 
    observables::Vector{<:AbstractMatrix}, meas0::Vector{<:AbstractFloat}, blocks::Vector{<:AbstractMatrix}, buff; fix_first_index::Bool = false)

    get_H_A!(buff.qcfl_buff.H_A, g, blocks)
    buff.qcfl_buff.H_A .*= -1im
    method_scalar, method_vector, need_buffer = integration_method
    if G !== nothing
        if need_buffer == true
            method_vector(buff.qcfl_buff.C_G_result, t -> integrand_gradient(t, ρ_A ,observables, meas0, blocks, buff, fix_first_index), T_max, buff.qcfl_buff.Σ)
        else 
            method_vector(buff.qcfl_buff.C_G_result, t -> integrand_gradient(t, ρ_A, observables, meas0, blocks, buff, fix_first_index), T_max)
        end 
        buff.qcfl_buff.C_G_result ./= length(observables)*T_max
        G[:] .= @view buff.qcfl_buff.C_G_result[2:end]
    end
    if F !== nothing 
        if G === nothing
            In = method_scalar(t -> integrand_cost(t, ρ_A, observables, meas0, buff), T_max)/(length(observables)*T_max)
            return In
        else 
            return buff.qcfl_buff.C_G_result[1]
        end 
    end
end 


function integrand_gradient(t::AbstractFloat, ρ_A::AbstractMatrix, observables::Vector{<:AbstractMatrix}, meas0::Vector{<:AbstractFloat}, blocks, buff, fix_first_index::Bool)
    
    qcfl_buff = buff.qcfl_buff

    qcfl_buff.H_A_forexp .=  t .* qcfl_buff.H_A 
    pullback =  own_rrule(exp, qcfl_buff.H_A_forexp, buff.exp_buff)
    
    evolve(ρ_A, buff.exp_buff.X, qcfl_buff)  
    
    δ = dev(observables[1], meas0[1], qcfl_buff)
    qcfl_buff.sumobs .= δ.*observables[1]
    qcfl_buff.C_G[1] = δ^2

    @fastmath @inbounds @simd for i in 2:lastindex(observables)
        δ = dev(observables[i], meas0[i], qcfl_buff)
        qcfl_buff.sumobs .+= δ.*observables[i]
        qcfl_buff.C_G[1] += δ^2
    end 

    pullback(mul!(qcfl_buff.dAforpb, qcfl_buff.ρ_A_right, qcfl_buff.sumobs))

    @fastmath @inbounds @simd for i in 1+fix_first_index:lastindex(qcfl_buff.C_G)-1+fix_first_index
        qcfl_buff.C_G[i+1-fix_first_index] = gradient_component(buff, blocks[i], t)
    end    
    
    return qcfl_buff.C_G
end

function integrand_cost(t::AbstractFloat,  ρ_A::AbstractMatrix, observables::Vector{<:AbstractMatrix}, meas0::Vector{<:AbstractFloat}, buff)

    qcfl_buff = buff.qcfl_buff

    qcfl_buff.H_A_forexp .=  t .* qcfl_buff.H_A 
    exp_only_buffered!(qcfl_buff.H_A_forexp, buff.exp_buff)

    evolve(ρ_A, buff.exp_buff.X, qcfl_buff)  

    c = 0.
    @fastmath @inbounds @simd for i in eachindex(observables)
        c += dev(observables[i], meas0[i], qcfl_buff)^2
    end 
    return c
end

function evolve(ρ_A::AbstractMatrix, U::AbstractMatrix, qcfl_buff)
    mul!(qcfl_buff.ρ_A_evolved, U, mul!(qcfl_buff.ρ_A_right, ρ_A, U'))
end 

function dev(observable::AbstractMatrix, expect_at_0::AbstractFloat, qcfl_buff)
    @inbounds mul!(qcfl_buff.evolob, observable, qcfl_buff.ρ_A_evolved)
    @inbounds δ = real(tr(qcfl_buff.evolob))- expect_at_0
    return δ
end 

function gradient_component(buff, block::AbstractMatrix, t::Real)
    qcfl_buff = buff.qcfl_buff
    frech = buff.exp_buff.∂X
    mul!(qcfl_buff.frechTimesBlock, frech, block)
    return 4*t*imag(tr(qcfl_buff.frechTimesBlock))
end


#########################################################################################
#                                                                                       #
#                   Commutator cost function                                            #
#                                                                                       #
#########################################################################################

"""
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
"""
function commutator(set::Settings, H_A::Function; buff::Union{Nothing, Commutator_buffer}=nothing)
    @unpack ρ_A = set
    
    blocks = H_A(set)

    buff = something(buff, make_commutator_buffer(set))

    return (F, G, g) -> comm_fg!(F, G, g, ρ_A, blocks, buff), NaN, length(blocks)
end 

function comm_fg!(F::Union{AbstractFloat, Nothing}, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, ρ_A::Matrix{ComplexF64}, blocks::Vector{<:AbstractMatrix}, buff)
 
    get_H_A!(buff.H_A, g, blocks)

    C = 0.
    if G !== nothing
        C = comm_grad!(G, ρ_A, blocks, buff)
    end
    if F !== nothing 
        if G === nothing
            C = comm_cost(ρ_A, buff)
            return C
        else 
            return C
        end 
    end
end 


function comm_cost(ρ_A::Matrix{ComplexF64}, buff)

    comm = buff.comm
    H_A = buff.H_A

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm
    
    return norm(comm)/(2*norm(H_A)*norm(ρ_A))
end 

function comm_grad!(G::Vector{<:AbstractFloat}, ρ_A::Matrix{ComplexF64}, blocks::Vector{<:AbstractMatrix}, buff)

    comm = buff.comm
    H_A = buff.H_A
    temp = buff.temp

    mul!(mul!(comm, ρ_A, H_A), H_A, ρ_A, 1, -1) # compute [H_A, ρ_A] and save it in comm


    Λ_comm = norm(comm)
    Λ_H_A = norm(H_A)
    Λ_ρ_A = norm(ρ_A)

    @fastmath @inbounds @simd for index in eachindex(G)
        mul!(mul!(temp, ρ_A, blocks[index]), blocks[index], ρ_A, 1, -1) # compute [h_index, ρ_A] and save it in temp
        ∂Λ_comm = 1/Λ_comm * real(sum(had!(temp, temp, comm)))
        ∂Λ_H_A = 1/Λ_H_A * real(sum(had!(temp, H_A, blocks[index])))
        G[index] = 1/(2*Λ_ρ_A*Λ_H_A^2)*(∂Λ_comm * Λ_H_A - Λ_comm * ∂Λ_H_A) 
    end 
    return Λ_comm/(2*Λ_H_A * Λ_ρ_A) 
end 

# compute hadamard product of two square matrices A and B and save it in buf
function had!(buf::AbstractMatrix, A::AbstractMatrix,B::AbstractMatrix)
    n = size(A)[1]
    for j in 1:n
       for i in 1:n
         @inbounds buf[i,j] = A[i,j] *B[i,j]
       end
    end
    return buf
end



#########################################################################################
#                                                                                       #
#                   Relative entropy                                                    #
#                                                                                       #
#########################################################################################



"""
    relative_entropy()

Constructs the commutator as cost function.

# Returns 
- `(F, G, g, init) -> relative_entropy_fg!(F, G, g, init), NaN`
    - `F`: cost function is not evaluated if `F = nothing`, is evaluated otherwise
    - `G`: vector to save the gradient in, if `G=nothing` gradient is not computed 
    - `g`: parameters to be optimized
    - `init`: see [`Init`](@ref)
- `NaN`: no parameter to be fixed. 

This function allows one to conveniently hand the cost function to be minimized to [`optimize_LBFGS`](@ref).
"""
function relative_entropy(set::Settings, H_A::Function; buff::Union{Nothing, Relative_entropy_buffer} = nothing)
    @unpack ρ_A = set

    blocks = H_A(set)

    buff = something(buff, make_relative_entropy_buffer(set))
    
    return (F, G, g) -> relative_entropy_fg!(F, G, g, ρ_A, blocks, buff), NaN, length(blocks)
end 

function relative_entropy_fg!(F::Union{AbstractFloat, Nothing}, G::Union{Vector{<:AbstractFloat}, Nothing}, g::Vector{<:Real}, 
    ρ_A::AbstractMatrix, blocks::Vector{<:AbstractMatrix}, buff)
 
    get_H_A!(buff.H_A, g, blocks)

    C = 0
    if G !== nothing
        C = relative_entropy_grad!(G, ρ_A, blocks, buff)
    end
    if F !== nothing 
        if G === nothing
            C = relative_entropy_cost(ρ_A, buff)
            return C
        else 
            return C
        end 
    end
end 

function relative_entropy_cost(ρ_A::AbstractMatrix, buff)
    
    buff.H_A_forexp  .= -1 .* buff.H_A
    exp_only_buffered!(buff.H_A_forexp, buff.exp_buff)
    exp_H_A = buff.exp_buff.X

    S1 = tr(mul!(buff.ρ_A_times_H_A, ρ_A, buff.H_A))
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 

function relative_entropy_grad!(G::AbstractVector, ρ_A::AbstractMatrix, blocks::Vector{<:AbstractMatrix}, buff)

    buff.H_A_forexp .= -1 .*buff.H_A


    ρ_A_times_block = buff.ρ_A_times_H_A # just a renaming so it's more readable

    exp_only_buffered!(buff.H_A_forexp, buff.exp_buff)
    exp_H_A = buff.exp_buff.X
    γ = -1/(tr(exp_H_A))

    @fastmath @inbounds @simd for i in eachindex(G)
        ∂S1 = tr(mul!(ρ_A_times_block, ρ_A, blocks[i]))
        ∂S2 = γ * tr(mul!(ρ_A_times_block, exp_H_A, blocks[i]))
        G[i] = ∂S1+∂S2
    end

    S1 = tr(mul!(ρ_A_times_block, ρ_A, buff.H_A))
    S2 = log(tr(exp_H_A))

    return S1 + S2
end 


#########################################################################################
#                                                                                       #
#                   INTEGRAND AFTER TANH-SINH TRANSFORMATION                            #
#                                                                                       #
#########################################################################################

#=

"""
    tanh_sinh_integrand(g::Vector{<:AbstractFloat}, u1::Real, u2::Real, num_points::Int, init::Init)

Return a vector containing the integrand after the  Tanh-sinh transformation at different points.
The integrand is sampled uniformly between `u1` and `u2` `num_points ` times.

# Arguments
- `g::Vector{<:AbstractFloat}`: parameters
- `u1::Real`: lowest sampling point 
- `u2::Real`: highest sampling point
- `num_points::Int`: number of points
- `init`: see [Init](@ref)
"""
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


#########################################################################################
#                                                                                       #
#                                COUNT INTEGRAND EVALUTAIONS                            #
#                                                                                       #
#########################################################################################

"""
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
"""
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



#########################################################################################
#                                                                                       #
#                                QUADGK                                                 #
#                                                                                       #
#########################################################################################

"""
    quadgk_count!(g::Vector{<:AbstractFloat}, init::Init; rtol=sqrt(eps), atol=0, maxevals=10^7, order=7)

Return the cost function value and the number of samples needed for integration with the [Gauss-Kronrod method
implemented in `quadgk.jl`](https://juliamath.github.io/QuadGK.jl/stable/api/#QuadGK.quadgk) as a tuple.

# Required Arguments
- `g::Vector{<:AbstractFloat}`: parameters
- `init`: see [Init](@ref)

# Keyword arguments 
- `atol=0.0`: absolute tolerance for integration 
- `rtol=sqrt(eps)`: relative tolerance for integration
- `maxevals=10^7` maximum number of samples 
- `order`: order for the quadrature  
"""
function quadgk_count!(g::Vector{<:AbstractFloat}, init::Init; rtol=sqrt(eps), atol=0, maxevals=10^7, order=7)
    @unpack observables, T_max, atol, rtol = init.set
    get_H_A!(g, init)
    init.buff.H_A .*= -1im
    I, E, count =  quadgk_count(t -> integrand_cost(t, init),0., T_max, rtol = rtol, atol = atol, maxevals=maxevals, order=order)
    return I/(length(observables)*T_max), count 
end

=#