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
QCFLBuffer
```

## Commutator

```@docs 
Commutator
```
```@docs
CommutatorBuffer
```

## Relative Entropy

```@docs 
RelativeEntropy
```
```@docs 
RelativeEntropyBuffer
```

## Wrapper for fixing parameters 

```@docs
FixedCost
```
```@docs
fill_full_g(fc::FixedCost, g::Vector{Float64}) 
```
## Gradient
```@docs
gradient(cost::AbstractCostFunction, g::Vector{<:Real})
```
```@docs
gradient!(G::Vector{<:Real}, cost::AbstractCostFunction, g::Vector{<:Real})
```

## Getter
```@docs
getmodel(cost::AbstractCostFunction)
```
```@docs 
getblocks(cost::AbstractCostFunction)
```
```@docs 
getobservables(cost::AbstractCostFunction)
```
## Printing 
```@docs 
print_model(cost::AbstractCostFunction)
```
```@docs 
print_blocks(cost::AbstractCostFunction)
```
```@docs 
print_observables(cost::AbstractCostFunction)
```