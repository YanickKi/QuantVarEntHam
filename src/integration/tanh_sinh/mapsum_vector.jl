@noinline function mapsum!(
    Σ, f, A::AbstractVector, ifirst::Integer, ilast::Integer, blksize::Int
)
    if ifirst == ilast
        @inbounds a1 = A[ifirst]
        return f(Σ, a1)
    elseif ifirst + blksize > ilast
        # sequential portion
        @inbounds a1 = A[ifirst]
        @inbounds a2 = A[ifirst + 1]
        f(Σ, a1)
        f(Σ, a2)
        @simd for i in (ifirst + 2):ilast
            @inbounds ai = A[i]
            f(Σ, ai)
        end
    else
        # pairwise portion
        imid = (ifirst + ilast) >> 1
        mapsum!(Σ, f, A, ifirst, imid, blksize)
        mapsum!(Σ, f, A, imid+1, ilast, blksize)
    end
end

@noinline function mapsum_withmax!(
    Σ::AbstractVector,
    f,
    A::AbstractVector,
    ifirst::Integer,
    ilast::Integer,
    blksize::Int,
    maximum::Float64,
)
    if ifirst == ilast
        @inbounds a1 = A[ifirst]
        max1 = f(Σ, a1)
        return max(maximum, max1)

    elseif ifirst + blksize > ilast
        # sequential portion
        @inbounds a1 = A[ifirst]
        @inbounds a2 = A[ifirst + 1]
        max1 = f(Σ, a1)
        max2 = f(Σ, a2)
        maximum = max(max1, max2, maximum)
        @simd for i in (ifirst + 2):ilast
            @inbounds ai = A[i]
            maxi = f(Σ, ai)
            maximum = max(maxi, maximum)
        end
        return maximum
    else
        # pairwise portion
        imid = (ifirst + ilast) >> 1
        max1 = mapsum_withmax!(Σ, f, A, ifirst, imid, blksize, maximum)
        max2 = mapsum_withmax!(Σ, f, A, imid+1, ilast, blksize, maximum)
        return max(max1, max2)
    end
end
#
#"""
#mapsum!(f, A::AbstractVector, ifirst::Integer=1, ilast::Integer=length(A))
#
#Return the total summation of items in `A` from `istart`-th through `iend`-th
#with applying a function `f`.
#
#NOTE: This function doesn't check `ifirst` and `ilast`. Be careful to use.
#"""
function mapsum!(
    Σ::AbstractVector, f, A::AbstractVector, ifirst::Integer=1, ilast::Integer=length(A)
)
    return mapsum!(Σ, f, A, ifirst, ilast, 512)
end

function mapsum_withmax!(
    Σ::AbstractVector, f, A::AbstractVector, ifirst::Integer=1, ilast::Integer=length(A)
)
    maximum::Float64 = 0
    return mapsum_withmax!(Σ, f, A, ifirst, ilast, 512, maximum)
end
