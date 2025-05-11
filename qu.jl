using QuantVarEntHam


function main()
    init = initialize(TFIM(8, 4, 1, 2), H_A_not_BW)
    giniit = [1,1,2,2,2,3,3]
    g,c = optimize_LBFGS(giniit, init)
    println(g/g[1])
end 

main()