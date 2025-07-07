# Integration 

 

```@docs
AbstractIntegrator
```

## Tanh-sinh quadrature 

```@docs
TanhSinh
```
```@docs
TanhSinh(; atol::Real=0.0, rtol::Real=atol > 0 ? 0. : sqrt(eps(Float64)), maxlevel::Integer=12, h0::Real=1)
```

## midpoint rule 
```@docs
MidPoint
```