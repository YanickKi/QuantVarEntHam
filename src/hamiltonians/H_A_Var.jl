using Yao 

mutable struct H_A_Var{T<:Real, D} <: CompositeBlock{D}
    n::Int 
    coeffs::Vector{T}  # Koeffizienten x1, x2, ...
    blocks::Vector{AbstractBlock{D}}
end



#Yao.mat(::Type{T}, c::H_A_Var{D})  where {T, D} = sum(c.coeffs.*mat.(T, c.blocks))
Yao.mat(::Type{T}, c::H_A_Var{D})  where {T, D} = mat(T, c.coeffs'*c.blocks)

niparams(::Type{<:H_A_Var{D}})  where {D} = n
getiparams(x::H_A_Var{D}) where {D} = x.coeffs
setiparams!(r::H_A_Var{D}, params::Vector{<:Real})  where {D} = (r.coeffs = params; r)
Base.adjoint(x::H_A_Var{D})  where {D} = x
Yao.nqudits(r::H_A_Var{D}) where {D}  = r.n
Yao.YaoAPI.subblocks(c::H_A_Var{T,D}) where {T, D} = c.coeffs.*c.blocks
function Yao.YaoAPI.unsafe_apply!(r::AbstractRegister, c::H_A_Var{T,D}) where {T,D}
    blocks = subblocks(c)
    res = mapreduce(blk -> YaoAPI.unsafe_apply!(copy(r), blk), regadd!, blocks[1:end-1])
    YaoAPI.unsafe_apply!(r, blocks[end])
    regadd!(r, res)
    r
    return r
end