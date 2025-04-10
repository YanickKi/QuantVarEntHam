@noinline function mapsum(f, A::AbstractVector, Σ,
    ifirst::Integer, ilast::Integer, blksize::Int,)
    if ifirst == ilast
        @inbounds a1 = A[ifirst]
        return f(a1,Σ)
    elseif ifirst + blksize > ilast
        # sequential portion
        @inbounds a1 = A[ifirst]
        @inbounds a2 = A[ifirst+1]
        f(a1,Σ) 
        f(a2,Σ)
        @simd for i in ifirst + 2 : ilast
        @inbounds ai = A[i]
        f(ai,Σ)
    end
    else
        # pairwise portion
        imid = (ifirst + ilast) >> 1
        mapsum(f, A, Σ, ifirst, imid, blksize)
        mapsum(f, A, Σ, imid+1, ilast, blksize)
    end
end

@noinline function mapsum_withmax(f, A::AbstractVector, Σ::AbstractVector,
    ifirst::Integer, ilast::Integer, blksize::Int, maximum::Float64)
    if ifirst == ilast
        @inbounds a1 = A[ifirst]
        max1 = f(a1, Σ)
        return max(maximum, max1)

    elseif ifirst + blksize > ilast
        # sequential portion
        @inbounds a1 = A[ifirst]
        @inbounds a2 = A[ifirst+1]
        max1 = f(a1, Σ) 
        max2 = f(a2, Σ)
        maximum = max(max1, max2, maximum)
        @simd for i in ifirst + 2 : ilast
        @inbounds ai = A[i]
        maxi = f(ai, Σ)
        maximum = max(maxi, maximum)
    end
    return maximum
    else
        # pairwise portion
        imid = (ifirst + ilast) >> 1
        max1 = mapsum_withmax(f, A, Σ, ifirst, imid, blksize, maximum)
        max2 = mapsum_withmax(f, A, Σ, imid+1, ilast, blksize, maximum)
    return max(max1, max2)
    end
end
#
#"""
#mapsum(f, A::AbstractVector, ifirst::Integer=1, ilast::Integer=length(A))
#
#Return the total summation of items in `A` from `istart`-th through `iend`-th
#with applying a function `f`.
#
#NOTE: This function doesn't check `ifirst` and `ilast`. Be careful to use.
#"""
function mapsum(f, A::AbstractVector, Σ::AbstractVector, ifirst::Integer=1, ilast::Integer=length(A))
    return mapsum(f, A, Σ, ifirst, ilast, 512)

end
 
function mapsum_withmax(f, Σ::AbstractVector, A::AbstractVector, ifirst::Integer=1, ilast::Integer=length(A))
    maximum::Float64 = 0 
    return mapsum_withmax(f, A, Σ, ifirst, ilast, 512, maximum)
end