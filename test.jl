using Distributed
using Yao 
using Parameters
using Cthulhu
include("src/QuantVarEntHam.jl")
using .QuantVarEntHam
using BenchmarkTools
using LinearAlgebra
using QuadGK


function mul_test()
    init = initialize(TFIM(14, 7, 1., 1.), H_A_BW)
    ginit = [1., 2., 3., 4., 5., 6., 7.]

    @btime QuantVarEntHam.cost_grad!(1., $ginit, zeros(length($ginit)), $init)

end 
#34.137 ms (13846 allocations: 64.84 MiB)
#= init = initialize(TFIM(14, 7, 1., 1.), H_A_BW), ginit = [1., 2., 3., 4., 5., 6., 7.]
vor der nutzung von der struktur der matrizen: 746.306 ms (14631 allocations: 1.05 GiB)
668.967 ms (20287 allocations: 1.05 GiB) matrobs mit mat aber ohne TS
661.708 ms (22105 allocations: 1.05 GiB) matrobs mit mat und buffer diagonal ohne TS
680.973 ms (14631 allocations: 1.05 GiB) matrobs mit mat und buffer diagonal mit TS
663.747 ms (16163 allocations: 1.05 GiB) auch für die matrices in blocks die mat funktion benutzt ohne TS
647.123 ms (14695 allocations: 1.05 GiB) done everything and everyting is type stable
679.109 ms (14695 allocations: 1.05 GiB) everyting type stable and parametric buffers, why slower???????
683.174 ms (14695 allocations: 1.05 GiB) before buffer for -1im*t*H_A
680.699 ms (14392 allocations: 1.03 GiB) after buffer for -1im*t*H_A
671.513 ms (14190 allocations: 1.03 GiB) after summing to c without sum and broadcasting 
650.395 ms (14129 allocations: 1.03 GiB) after making blocks dense again 
649.219 ms (14126 allocations: 1.03 GiB) no Hermitian Matrix constructor in cost_grad anymore 
665.716 ms (14108 allocations: 1.03 GiB) loop to make H_A not the sum Function
660.618 ms (14087 allocations: 1.03 GiB) working in place for skalar times matrix
=#
#Everything with T_max = 2. here 
#before optimizing 40.244 ms (17738 allocations: 66.21 MiB) (type instability)
#fixing bug with static arrays but type instable: 39.218 ms (14308 allocations: 66.14 MiB)
#using another buffer 39.284 ms (14294 allocations: 66.14 MiB)
# for T_max = 1 : 18.798 ms (7047 allocations: 32.49 MiB)
mul_test()