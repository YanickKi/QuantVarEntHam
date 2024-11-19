using Yao 

mutable struct LinearCombinationBlock{T<:Real, D} <: CompositeBlock{D}
    n::Int 
    coeffs::Vector{T}  # Koeffizienten x1, x2, ...
    blocks::Vector{AbstractBlock{D}}
end



Yao.mat(::Type{T}, c::LinearCombinationBlock{D})  where {T, D} = sum(c.coeffs.*mat.(T, c.blocks))

niparams(::Type{<:LinearCombinationBlock{D}})  where {D} = 5
getiparams(x::LinearCombinationBlock{D}) where {D} = x.coeffs
setiparams!(r::LinearCombinationBlock{D}, params::Vector{<:Real})  where {D} = (r.coeffs = params; r)
Base.adjoint(x::LinearCombinationBlock{D})  where {D} = x
Yao.nqudits(r::LinearCombinationBlock{D}) where {D}  = r.n
Yao.YaoAPI.subblocks(c::LinearCombinationBlock{T,D}) where {T, D} = c.coeffs.*c.blocks
test = LinearCombinationBlock{Float64, 2}(2, [0.,1.], [repeat(2,Z,(1,2)), put(2, 2=>Z)])


#reg = ArrayReg([0, 1, -1+0.0im, 0])
#
#apply!(reg, test)
#
#
#print_table(reg)
#
#
#TE = time_evolve(test, 1.)
#print(mat(TE))