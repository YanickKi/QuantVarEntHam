function H_A_BW_wo_corrections!(blocks::Vector{<:Block}, model::AbstractModel{S,N_A}) where {S,N_A}
    for i in 1:N_A 
        push!(blocks, hi(model, i))
    end     
end


function corrections!(blocks::Vector{<:Block}, model::AbstractModel{S,N_A}) where {S,N_A}
    r_max = model.r_max
    for r in 2:r_max
        for i in 1:N_A-r
            correction!(blocks, model, i, r)
        end
    end 
end 