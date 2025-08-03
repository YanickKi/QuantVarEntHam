# COV_EXCL_START
"""
    getansatz(cost::AbstractCostFunction)

Return the ansatz of a given cost function.
"""
getansatz(cost::AbstractCostFunction) = getansatz(cost)
getansatz(cost::AbstractFreeCostFunction) = cost.ansatz
getansatz(fc::FixedCost) = unwrap(getansatz, fc)

"""
    getobservables(cost::AbstractCostFunction)

Return the observables of a given cost function.
"""
getobservables(cost::AbstractCostFunction) = getobservables(hasobservables(cost), cost)
getobservables(::HasObservables, cost::AbstractFreeCostFunction) = cost.observables
function getobservables(::HasNoObservables, cost::AbstractFreeCostFunction)
    error("The cost function $(nameof(typeof(cost))) has no observables!")
end
getobservables(fc::FixedCost) = unwrap(getobservables, fc)

"""
    getmodel(cost::AbstractCostFunction)

Return the model of a given cost function.
"""
getmodel(cost::AbstractCostFunction) = getmodel(cost)
getmodel(cost::AbstractFreeCostFunction) = cost.model
getmodel(fc::FixedCost) = unwrap(getmodel, fc)
# COV_EXCL_STOP
