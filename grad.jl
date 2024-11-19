using LinearAlgebra
using Zygote
using Yao 

# Definition der Matrix A als Funktion von g und h
function A(g, h)
    return [g h; 0.0 g * h]
end

# Matrixexponential als Funktion von g und h
function expA(g, h)
    return exp(A(g, h))
end

# Berechnung der Jacobi-Matrix in Bezug auf (g, h)
function jacobian_expA(g, h)
    jacobian(x -> expA(x[1], x[2]), [g, h])
end

# Beispiel für g = 1.0, h = 2.0
#g = 1.0
#h = 2.0
#expA_result = expA(g, h)
#jacobian_result = jacobian_expA(g, h)

#println("Matrixexponential exp(A(g, h)) für g = $g, h = $h:\n", expA_result)
#println("Jacobi-Matrix der Ableitung von exp(A(g, h)) in Bezug auf (g, h):\n", jacobian_result)

N = 2
h1 = (put(4,1=>Z))
h2 = (put(4,2=>Z))
H(g1,g2) = [g1*h1, g2*h2]
reg =  (density_matrix(product_state(bit"1010"))).state

function expH(g1, g2)
    return exp(2.0*1im*Matrix(sum(Matrix.(H(g1, g2)))))
end


function ew(g1, g2, reg)
    U = expH(g1, g2)
    return tr(reg'*U*reg)
end

g1 = 1.
g2 = 2. 
function grad(g1, g2, reg)
    r = gradient(x -> real(ew(x[1],x[2], reg)), [g1,g2])
    i = gradient(x -> imag(ew(x[1],x[2], reg)), [g1,g2])
    return r, i
end

#expH_result = expH(g1, g2)
#r, i = jacobian_expH(g1,g2)
#@time jacobian_expH(g1,g2)

res = ew(g1,g2,reg)
println(typeof(real(res)))
println(typeof(imag(res)))


r, i = grad(g1, g2, reg)
@time grad(g1, g2, reg)
println(r)
println(i)
#println(expH_result)
#println(r)
#println(i)