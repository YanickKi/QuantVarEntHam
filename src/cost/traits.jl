# COV_EXCL_START
struct HasObservables end
struct HasNoObservables end

hasobservables(::QCFL) = HasObservables()
hasobservables(::AbstractCostFunction) = HasNoObservables()
# COV_EXCL_STOP
