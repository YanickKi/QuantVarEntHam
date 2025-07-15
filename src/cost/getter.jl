"""
    getblocks(cost::AbstractCostFunction)

Return the blocks of a given cost function.
"""
getblocks(cost::AbstractCostFunction) = getblocks(cost)
getblocks(cost::AbstractFreeCostFunction) = cost.blocks
getblocks(fc::FixedCost) = unwrap(getblocks, fc)

"""
    getobservables(cost::AbstractCostFunction)

Return the observables of a given cost function.
"""
getobservables(cost::AbstractCostFunction) = getobservables(cost)
getobservables(cost::AbstractFreeCostFunction) = cost.observables
getobservables(fc::FixedCost) = unwrap(getobservables, fc)


"""
    getmodel(cost::AbstractCostFunction)

Return the model of a given cost function.
"""
getmodel(cost::AbstractCostFunction) = getmodel(cost)
getmodel(cost::AbstractFreeCostFunction) = cost.model
getmodel(fc::FixedCost) = unwrap(getmodel, fc)