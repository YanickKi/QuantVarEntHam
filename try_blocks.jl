using QuantVarEntHam
using Cthulhu

function main()
    N = 5
    S = 1//1
    ps = PauliString(N, "Z", 1)
    @descend mat(ps)
end 



main()