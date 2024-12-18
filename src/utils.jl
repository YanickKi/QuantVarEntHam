using LinearAlgebra, CairoMakie, LaTeXStrings

function calc_universal_ratios(ξ::Vector{<:AbstractFloat}, α0::Integer, α1::Integer)
    κ = (ξ .- ξ[α0])/(ξ[α1] - ξ[α0])
    return κ
end

function universal_ratios(g::Vector{<:AbstractFloat}, init::Init; α0::Integer = 1, α1::Integer = 5)
    @unpack ρ_A = init.set
    H_A = H_A = @inbounds sum(g[i].*init.blks.matrices[i] for i in eachindex(g))
    H_A_exact = -log(Hermitian(ρ_A.state))
    ξ_var, v = eigen(Hermitian(H_A))
    ξ_exact, v =  eigen(Hermitian(H_A_exact))
    return calc_universal_ratios(ξ_var, α0, α1), calc_universal_ratios(ξ_exact, α0, α1)
end

function universal_ratios_svd(g::Vector{<:AbstractFloat}, init::Init; α0::Integer = 1, α1::Integer = 5)
    @unpack ρ_A = init.set
    H_A = H_A = @inbounds sum(g[i].*init.blks.matrices[i] for i in eachindex(g))
    F = svd(ρ_A.state)
    probs = F.S 
    ξ_var, v = eigen(Hermitian(H_A))
    ξ_exact = -log.(probs.^2)
    return calc_universal_ratios(ξ_var, α0, α1), calc_universal_ratios(ξ_exact, α0, α1)
end

function print_H_A(g::Vector{<:AbstractFloat}, init::Init)
    print(sum(g[i]*init.blks.blocks[i] for i in eachindex(g)))
end

function get_H_A(g::Vector{<:AbstractFloat}, init::Init)
    return sum(g[i].*init.blks.matrices[i] for i in eachindex(g))
end 

function plot_universal_ratios(filename::String, universal_ratios...)
    F = Figure()
    ax = Axis(F[1,1], xlabel = L"α", ylabel = L"κ_α", xgridvisible = false, ygridvisible = false)
    for i in eachindex(universal_ratios)
        scatter!(ax, universal_ratios[i])
    end 
    save(filename, F)
end 

function plot_universal_ratios(filename::String, labels::Union{Vector{LaTeXStrings.LaTeXString}, Vector{String}}, universal_ratios...)
    F = Figure()
    ax = Axis(F[1,1], xlabel = L"α", ylabel = L"κ_α", xgridvisible = false, ygridvisible = false,)
    for i in eachindex(universal_ratios)
        scatter!(ax, universal_ratios[i], label = labels[i])
    end 
    axislegend(ax, position = :lt)
    save(filename, F)
end 