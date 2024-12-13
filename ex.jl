using LinearAlgebra
function whole()
    A = rand(ComplexF64, 2,2)
    B = (A+A')/2
    println(typeof(Hermitian(B)))

end 

whole()