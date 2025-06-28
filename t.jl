using QuantVarEntHam

mp = MidPoint(1e-2)

mp.scalar_integrate(t -> sin(t), pi)

I = zeros(2)

mp.vector_integrate(I, t -> [sin(t), cos(t)], pi)

I