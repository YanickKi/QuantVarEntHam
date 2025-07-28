abstract type BufferTrait end
struct NeedBuffer <: BufferTrait end
struct NoNeedBuffer <: BufferTrait end

shorten_buffer!(::NoNeedBuffer, ::AbstractVectorIntegrator, ::Integer) = nothing

function shorten_buffer!(
    ::NeedBuffer, vector_integrate::AbstractVectorIntegrator, how_often::Integer
)
    for _ in 1:how_often
        pop!(vector_integrate.buffer)
    end
end
