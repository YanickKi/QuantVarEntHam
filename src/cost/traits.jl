struct HasObservables end
struct HasNoObservables end

hasobservables(::QCFL) = HasObservables()
hasobservables(::AbstractCostFunction) = HasNoObservables()
