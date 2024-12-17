# Optimizer 

```@docs
optimize_LBFGS(g_init::Vector{<:AbstractFloat}, init::Init; g1::AbstractFloat=NaN, gtol::AbstractFloat=1e-12, maxiter::Integer = 200, show_trace::Bool=true)
```