using QuantumOpticsBase


function repeat_cust(N, sig::Function, loc::AbstractVector{<:Int}; S::Rational=1//1)
    b = SpinBasis(S) 

    Ops = []
    for qudit in 1:N
        if qudit ∉ loc 
            push!(Ops, sparse(identityoperator(b)))
        end
        if qudit ∈ loc 
            push!(Ops, sparse(sig(b)))
        end
    end
    TensOp = tensor(Ops...)
    return TensOp
end 



function main()

    M = repeat_cust(2, sigmay, [1,2])
    println(M.data)
end 

main()