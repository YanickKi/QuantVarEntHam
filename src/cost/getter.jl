getblocks(cost::AbstractFreeCostFunction) = cost.blocks
getobservables(cost::AbstractFreeCostFunction) = cost.observables

getblocks(fc::FixedCost) = unwrap(getblocks, fc)
getobservables(fc::FixedCost) = unwrap(getobservables, fc)
    
