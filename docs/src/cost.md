# Cost

```@docs 
AbstractCostFunction
```
```@docs 
AbstractFreeCostFunction
```
## Quantum Classical Feedback Loop (QCFL) 

```@docs 
QCFL
```
```@docs 
QCFL(model::AbstractModel, blocks::Vector{<:AbstractMatrix}, T_max::Real; integrator::Union{Nothing,AbstractIntegrator} = nothing, observables::Union{Nothing, Vector{<:AbstractMatrix}} = nothing,
buffer::Union{Nothing, QCFL_buffer} = nothing) 
```
```@docs
QCFLBuffer
```
```@docs
QCFLBuffer(model::AbstractModel, blocks::Vector{<:AbstractMatrix}, observables::Vector{<:AbstractMatrix})
```

## Commutator

```@docs 
Commutator
```
```@docs 
Commutator(model::AbstractModel, blocks::Vector{<:AbstractMatrix})
```
```@docs
CommutatorBuffer
```
```@docs
CommutatorBuffer(model::AbstractModel)
```

## Relative Entropy

```@docs 
RelativeEntropy
```
```@docs 
RelativeEntropy(model::AbstractModel, blocks::Vector{<:AbstractMatrix})
```
```@docs
RelativeEntropyBuffer
```
```@docs
RelativeEntropyBuffer(model::AbstractModel)
```

## Wrapper for fixing parameters 

```@docs
FixedCost
```
```@docs
FixedCost(c::AbstractFreeCostFunction, fixed_indices::Vector{<:Integer}, fixed_values::Vector{<:Real})
```
```@docs
fill_full_g(fc::FixedCost, g::Vector{Float64}) 
```
## Gradient

```@docs
gradient!(G::Vector{<:Real}, c::AbstractCostFunction, g::Vector{<:Real})
```

