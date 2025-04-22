function repeat_cust(N, sig::Function, loc::AbstractVector{<:Int}; S::Rational=1//1)
    b = SpinBasis(S) 

    Ops = []
    for qudit in 1:N
        if qudit ∉ loc 
            push!(Ops, identityoperator(b))
        end
        if qudit ∈ loc 
            push!(Ops, sig(b))
        end
    end
    TensOp = tensor(Ops...)
    return TensOp
end 

function repeat_diag(N, loc::AbstractVector{<:Int}; S::Rational=1//1)
    b = SpinBasis(S) 

    Ops = []
    for qudit in 1:N
        if qudit ∉ loc 
            push!(Ops, identityoperator(b))
        end
        if qudit ∈ loc 
            push!(Ops, sig(b))
        end
    end
    TensOp = tensor(Ops...)
    return TensOp
end 
