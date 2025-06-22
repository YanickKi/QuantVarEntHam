using QuantVarEntHam


function main()
    init = initialize(TFIM(8, 4, 1, 1), H_A_BW)
    ginit = [1,2,3,4]
    g,c = QuantVarEntHam.optimize_LBFGS(ginit, init)
    #I, count = QuantVarEntHam.cost_count!(g, init)
    #println(I)
    #println(count)
    #g = 1e-3*[0.5786116515453814, 1.137519377265569, 1.6576902448971031, 2.12141045941895, 2.5128885850966975, 2.8187933038866, 3.02870739748262]
    #println(init.buff.count)
end

main()