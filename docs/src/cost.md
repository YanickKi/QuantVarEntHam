# Cost

```@docs 
QCFL(;integration_method = tanh_sinh(), g1::Real=NaN)
```
```@docs 
commutator()
```
```@docs 
relative_entropy()
```
```@docs 
tanh_sinh_integrand(g::Vector{<:AbstractFloat}, u1::Real, u2::Real, num_points::Int, init::Init)
```
```@docs 
cost_count(g::Vector{<:AbstractFloat}, init::Init; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Int = 12, h0::Real = 1.)
```
```@docs 
quadgk_count!(g::Vector{<:AbstractFloat}, init::Init; rtol=sqrt(eps), atol=0, maxevals=10^7, order=7)
```