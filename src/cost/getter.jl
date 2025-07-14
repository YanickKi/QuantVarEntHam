getblocks(cost::AbstractFreeCostFunction) = cost.str_blocks
getobservables(cost::AbstractFreeCostFunction) = cost.str_observables

getblocks(fc::FixedCost) = unwrap(getblocks, fc)
getobservables(fc::FixedCost) = unwrap(getobservables, fc)
    
