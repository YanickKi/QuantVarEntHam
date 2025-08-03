abstract type BufferTrait end
struct NeedBuffer <: BufferTrait end
struct NoNeedBuffer <: BufferTrait end

shorten_buffer!(::NoNeedBuffer, ::AbstractVectorIntegrator, ::Integer) = nothing

function shorten_buffer!(
    ::NeedBuffer, vector_integrator::AbstractVectorIntegrator, how_often::Integer
)
    v = vector_integrator.buffer
    deleteat!(v, (lastindex(v) - how_often + 1):lastindex(v))
end
