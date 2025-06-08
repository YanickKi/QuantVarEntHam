using LinearAlgebra
import LinearAlgebra.BlasFloat
using ChainRulesCore: add!!

mutable struct ExpBuffer{T}
    Inn::Matrix{T}
    A::Matrix{T}
    A2::Matrix{T}
    P::Matrix{T}
    tempP::Matrix{T}
    W::Matrix{T}
    V::Matrix{T}
    tempV::Matrix{T}
    U::Matrix{T}
    X::Matrix{T}
    tempX::Matrix{T}
    temp2X::Matrix{T}
    VminU::Matrix{T}
    ∂A2::Matrix{T}
    ΔAA::Matrix{T}
    AΔA::Matrix{T}
    ∂temp::Matrix{T}
    ∂X::Matrix{T}
    ∂W::Matrix{T}
    ∂U::Matrix{T}
    Apows::Vector{Matrix{T}}
    Xpows::Vector{Matrix{T}}
    C::Vector{Vector{T}}
end

function ExpBuffer(T::Type, n::Int)
    C1 = T[17643225600.0,8821612800.0,2075673600.0,302702400.0,30270240.0,2162160.0,110880.0,3960.0,90.0,1.0,]
    C2 = T[17297280.0, 8648640.0, 1995840.0, 277200.0, 25200.0, 1512.0, 56.0, 1.0]
    C3 = T[30240.0, 15120.0, 3360.0, 420.0, 30.0, 1.0]
    C4 = T[120.0, 60.0, 12.0, 1.0]
    C5 = T[64764752532480000.0,32382376266240000.0,7771770303897600.0,1187353796428800.0,129060195264000.0,10559470521600.0,670442572800.0,33522128640.0,1323241920.0,40840800.0,960960.0,16380.0,182.0,1.0,]
    return ExpBuffer(
        Matrix{T}(I, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        zeros(T, n, n),
        Matrix{T}[zeros(T,n,n) for i in 1:6],  # Apows
        Matrix{T}[zeros(T,n,n) for i in 1:5],  # Xpows
        [C1,C2,C3,C4,C5]
    )
end

function exp_bufered!(A::StridedMatrix{T}, buf::ExpBuffer{T}) where {T<:BlasFloat}    
    n = LinearAlgebra.checksquare(A)
    ilo, ihi, scale = LAPACK.gebal!('B', A)  # modifies A
    nA = opnorm(A, 1)
    ## For sufficiently small nA, use lower order Padé-Approximations
    if (nA <= 2.1)
        if nA > 0.95
            C = buf.C[1]
        elseif nA > 0.25
            C = buf.C[2]
        elseif nA > 0.015
            C = buf.C[3]
        else
            C = buf.C[4]
        end
        si = 0
    else
        C = buf.C[5]
        s = log2(nA / 5.4)  # power of 2 later reversed by squaring
        si = ceil(Int, s)
    end
    if si > 0
        A ./= convert(T, 2^si)
    end

    mul!(buf.A2, A, A)  # A2 = A * A
    buf.P .= buf.Inn  # P = I
    buf.W .= C[2] .* buf.P
    buf.V .= C[1] .* buf.P
    sizeApows = (div(size(C, 1), 2) - 1)
    for k in 1:sizeApows
        k2 = 2 * k
        mul!(buf.tempP, buf.P, buf.A2)
        buf.P .= buf.tempP
        buf.Apows[k] .= buf.P # CARE, APOWS MIGHT BE TOO SHORT
        buf.W .+= C[k2 + 2] .* buf.P
        buf.V .+= C[k2 + 1] .* buf.P
    end
    mul!(buf.U, A, buf.W)
    buf.X .= buf.V .+ buf.U
    buf.V .= buf.V .- buf.U
    buf.tempV .= buf.V
    F = lu!(buf.tempV)  # NOTE: use lu! instead of LAPACK.gesv! so we can reuse factorization
    ldiv!(F, buf.X)
    buf.Xpows[1] .= buf.X
    if si > 0  # squaring to reverse dividing by power of 2
        for t in 1:si
            mul!(buf.tempX, buf.X, buf.X)
            buf.X .= buf.tempX
            buf.Xpows[t+1] .= buf.X # CARE, XPOWS MIGHT BE TOO SHORT
        end
    end
    sizeXpows = si+1
    _unbalance!(buf.X, ilo, ihi, scale, n)
    return (ilo, ihi, scale, C, si, F, sizeApows, sizeXpows)
end

function _balance!(X, ilo, ihi, scale, n)
    n = size(X, 1)
    if ihi < n
        for j in (ihi + 1):n
            LinearAlgebra.rcswap!(j, Int(scale[j]), X)
        end
    end
    if ilo > 1
        for j in (ilo - 1):-1:1
            LinearAlgebra.rcswap!(j, Int(scale[j]), X)
        end
    end

    for j in ilo:ihi
        scj = scale[j]
        for i in 1:n
            X[j, i] /= scj
        end
        for i in 1:n
            X[i, j] *= scj
        end
    end
    return X
end

function _unbalance!(X, ilo, ihi, scale, n)
    for j in ilo:ihi
        scj = scale[j]
        for i in 1:n
            X[j, i] *= scj
        end
        for i in 1:n
            X[i, j] /= scj
        end
    end

    if ilo > 1
        for j in (ilo - 1):-1:1
            LinearAlgebra.rcswap!(j, Int(scale[j]), X)
        end
    end
    if ihi < n
        for j in (ihi + 1):n
            LinearAlgebra.rcswap!(j, Int(scale[j]), X)
        end
    end
    return X
end


function own_rrule(::typeof(exp), A0::StridedMatrix{<:BlasFloat}, buf::ExpBuffer)
    # TODO: try to make this more type-stable
    buf.A .= A0
    intermediates = exp_bufered!(buf.A, buf)
    function exp_pullback(ΔX)
        _matfun_frechet!(exp, ΔX, buf.A, intermediates, buf)
    end
    return exp_pullback
end


function _matfun_frechet!(
    ::typeof(exp), ΔA, A::StridedMatrix{T}, (ilo, ihi, scale, C, si, F, sizeApows, sizeXpows), buf::ExpBuffer
) where {T<:BlasFloat}
    n = LinearAlgebra.checksquare(A)
    _balance!(ΔA, ilo, ihi, scale, n)
    if si > 0
        ΔA ./= convert(T, 2^si)
    end
    
    ∂A2 = buf.∂A2
    mul!(buf.ΔAA, ΔA ,A)
    mul!(buf.AΔA, A  ,ΔA)
    ∂A2 .= buf.ΔAA .+ buf.AΔA 
    A2 = first(buf.Apows)
    # we will repeatedly overwrite ∂temp and ∂P below
    
    ∂temp = buf.∂temp
    ∂P = buf.P 
    ∂W = buf.∂W
    ∂V = buf.V
    ∂P .= ∂A2
    ∂W .= C[4] .* ∂P
    ∂V .= C[3] .* ∂P
    for k in 2:(sizeApows-1)
        k2 = 2 * k
        P = buf.Apows[k - 1]
        #∂P, ∂temp = mul!(mul!(∂temp, ∂P, A2), P, ∂A2, true, true), ∂P
        buf.ΔAA .=  ∂P
        mul!(∂temp, ∂P, A2)
        mul!(∂temp, P, ∂A2, true, true)
        ∂P .= ∂temp
        ∂temp .= buf.ΔAA
        axpy!(C[k2 + 2], ∂P, ∂W)
        axpy!(C[k2 + 1], ∂P, ∂V)
    end
 
    #∂U, ∂temp = mul!(mul!(∂temp, A, ∂W), ΔA, buf.W, true, true), ∂W
    ∂U = buf.∂U
    buf.ΔAA .= ∂W
    mul!(∂temp, A, ∂W)
    mul!(∂temp, ΔA, buf.W, true, true)
    ∂U .= ∂temp
    ∂temp .= ∂U .- ∂V
    ∂X = buf.∂X
    ∂X .= add!!(∂U, ∂V)
    mul!(∂X, ∂temp, first(buf.Xpows), true, true)
    ldiv!(F, ∂X)
    X = buf.X
    if si > 0
        for t in 1:sizeXpows-1
            X = buf.Xpows[t]
            #∂X, ∂temp = mul!(mul!(∂temp, X, ∂X), ∂X, X, true, true), ∂X
            buf.ΔAA .= ∂X
            mul!(∂temp, X, ∂X)
            mul!(∂temp, ∂X, X, true, true)
            ∂X .= ∂temp
            ∂temp .=  buf.ΔAA 
        end
    end
    _unbalance!(∂X, ilo, ihi, scale, n)
end