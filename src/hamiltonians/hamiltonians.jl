abstract type Settings end 
using KrylovKit: eigsolve



function get_rhoA(H::AbstractBlock, A::AbstractRange) 
    values, vectors = eigsolve(mat(H) ,1 ,:SR, ishermitian=true)
    rhoA = density_matrix(ArrayReg(vectors[1]), A)
    return rhoA
end 

mutable struct H_A_Var
    blocks::Vector{AbstractBlock}
    matrices::Vector{Matrix{ComplexF64}}
end 


function H_A_BW(set::Settings) 
    @unpack N, N_A, r_max, periodic = set
    
    if 2*N_A != N && periodic == false 
        @warn "Be aware: The Bisognano-Wichmann theorem for the case of open boundary conditions is only valid for N = 2*N_A i.e. for a half plane!" 
    end 
    
    blks = Vector{AbstractBlock}(undef, 0)
    
    H_A_BW_wo_corrections!(blks, set)
    
    if r_max > 1
        corrections!(blks, set)
    end 

    return H_A_Var(blks, Matrix.(blks))

end 

function H_A_not_BW(set::Settings) 
    @unpack N_A, r_max = set
    
    blks = Vector{AbstractBlock}(undef, 0)
    
    H_A_notBW_wo_corrections!(blks, set)
    
    if r_max > 1
        corrections!(blks, set)
    end 

    return H_A_Var(blks, Matrix.(blks))

end 



function H_A_BW_wo_corrections!(blks::Vector{AbstractBlock}, set::Settings)
    @unpack N_A = set
    for i in 1:N_A 
        push!(blks, hi(i, set))
    end     
end


function corrections!(blks::Vector{AbstractBlock}, set::Settings)
    @unpack N_A, r_max = set
    for r in 2:r_max
        for i in 1:N_A-r
            correction!(blks, i, r, set)
        end
    end 
end 

function print_H_A_Var(g::Vector{<:AbstractFloat}, blks::H_A_Var)
    println(sum(g[i] * blks.blocks[i] for i in eachindex(g)))
end 

include("xxz.jl")
include("tfim.jl")




#####################################################
#                                                   #
#           OLD STUPID STUFF BELOW HERE             #
#                                                   #
#####################################################



#=
function H_A_BW(g::Vector{<:AbstractFloat}, set::Settings, hi::Function, correction::Function, params_per_order::Integer)
    @unpack N_A, r_max = set
    
    H_A::Add{2} = map(1:N_A) do i 
        g[i]*hi(i,set)
    end |> sum

    if r_max > 1
        H_A+=corrections(g, set, correction, params_per_order, N_A+1)
    end 
    return H_A 
end


function H_A_notBW(g::Vector{<:AbstractFloat}, set::Settings, H_A_notBW_wo_corretions::Function, correction::Function, params_per_order::Integer, start_count::Integer)
    @unpack N_A, r_max, = set
    H_A  = H_A_notBW_wo_corretions(g, set)
    if r_max > 1
        H_A+=corrections(g, set, correction, params_per_order, start_count)
    end
    return H_A
end

function corrections(g::Vector{<:AbstractFloat}, set::Settings, correction::Function, params_per_order::Integer, start_count::Integer)
    @unpack N_A, r_max = set
    count = start_count
    for r in 2:r_max
        for i in 1:N_A-r
            for par in 1:params_per_order
                push!(corr,g[count]*correction(i, r, par, set))
                count += 1
            end
        end
    end
    return sum(corr)
end 
=#