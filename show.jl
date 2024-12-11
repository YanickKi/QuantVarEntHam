using QuantVarEntHam

function sh()
    init = initialize(XXZ(10, 5, 1., 1.), H_A_BW)
    g = [1., 2., 3., 4., 5.]
    optimize_LBFGS(g, init, gtol = 1e-16)    
end 

sh()