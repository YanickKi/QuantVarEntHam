using QuantVarEntHam

model = TFIM(8,4,1)

ansatz = H_A_BW(model)

max_times = 0.1:0.1:10 

g_init = [1,2,3,4]

buffer = QCFLBuffer(ansatz) # predefine buffer once

for T_max in max_times 
    cost = QCFL(model, ansatz, T_max, buffer=buffer) # reuse buffer in each iteration
    optimize(cost, g_init)
end