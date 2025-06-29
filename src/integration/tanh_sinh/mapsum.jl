@noinline function mapsum(f, A::AbstractVector,
    ifirst::Integer, ilast::Integer, blksize::Int)
    if ifirst == ilast
        @inbounds a1 = A[ifirst]
        return f(a1)
    elseif ifirst + blksize > ilast
        # sequential portion
        @inbounds a1 = A[ifirst]
        @inbounds a2 = A[ifirst+1]
        v = f(a1) + f(a2)
        @simd for i in ifirst + 2 : ilast
        @inbounds ai = A[i]
        v = v + f(ai)
    end
    return v
    else
        # pairwise portion
        imid = (ifirst + ilast) >> 1
        v1 = mapsum(f, A, ifirst, imid, blksize)
        v2 = mapsum(f, A, imid+1, ilast, blksize)
    return v1 + v2
    end
end

@noinline function mapsum_withmax(f, A::AbstractVector,
    ifirst::Integer, ilast::Integer, blksize::Int, maximum::Float64)
    if ifirst == ilast
        @inbounds a1 = A[ifirst]
        f1, max1 = f(a1)
        return f1, max(maximum, max1)
    elseif ifirst + blksize > ilast
        # sequential portion
        @inbounds a1 = A[ifirst]
        @inbounds a2 = A[ifirst+1]
        v1, max1 = f(a1) 
        v2, max2 = f(a2)
        maximum = max(max1, max2, maximum)
        v = v1+v2
        @simd for i in ifirst + 2 : ilast
        @inbounds ai = A[i]
        fi, maxi = f(ai)
        maximum = max(maxi, maximum)
        v = v + fi
    end
    return v, maximum
    else
        # pairwise portion
        imid = (ifirst + ilast) >> 1
        v1, max1 = mapsum_withmax(f, A, ifirst, imid, blksize, maximum)
        v2, max2 = mapsum_withmax(f, A, imid+1, ilast, blksize, maximum)

    return v1 + v2, max(max1, max2)
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
function mapsum(f, A::AbstractVector, ifirst::Integer=1, ilast::Integer=length(A))
    return mapsum(f, A, ifirst, ilast, 512)
    #if Threads.nthreads() == 1 || multithreading == false
    #    return mapsum(f, A, ifirst, ilast, 512)
    #end
    #tasks = map(index_chunks(ifirst:ilast; n=1*Threads.nthreads())) do chunk
    #    Threads.@spawn mapsum(f, A, chunk[1], chunk[end], 512)
    #end::Vector{Task}
    #chunk_sums::Vector{T} = fetch.(tasks)
    #return sum(chunk_sums)
end
 
function mapsum_withmax(f, A::AbstractVector, ifirst::Integer=1, ilast::Integer=length(A))
    maximum::Float64 = 0 
    return mapsum_withmax(f, A, ifirst, ilast, 512, maximum)
    #if Threads.nthreads() == 1 || multithreading == false
    #    return mapsum(f, A, ifirst, ilast, 512)
    #end
    #tasks = map(index_chunks(ifirst:ilast; n=1*Threads.nthreads())) do chunk
    #    Threads.@spawn mapsum(f, A, chunk[1], chunk[end], 512)
    #end::Vector{Task}
    #chunk_sums::Vector{T} = fetch.(tasks)
    #return sum(chunk_sums)
end