# QuantVarEntHam

[![Build Status](https://gitlab.com/quantum-computing-software/QuantVarEntHam/-/badges/release.svg)](https://gitlab.com/quantum-computing-software/QuantVarEntHam)
[![pipeline status](https://gitlab.com/quantum-computing-software/QuantVarEntHam/badges/main/pipeline.svg)](https://gitlab.com/quantum-computing-software/QuantVarEntHam/-/pipelines)
[![coverage report](https://gitlab.com/quantum-computing-software/QuantVarEntHam/badges/main/coverage.svg)](https://gitlab.com/quantum-computing-software/QuantVarEntHam/-/graphs/main/charts)
[![documentation (placeholder)](https://img.shields.io/badge/docs-latest-blue.svg)](https://quantum-computing-software.gitlab.io/QuantVarEntHam/)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)


This package aims to find the Entanglement Hamiltonian (EH) of one dimensional spin lattice models with arbitrary spin. 
Kokail et al. presented the algorithm -a Quantum Classical Feedback Loop (QCFL)- , which is the foundation of this package.
It learns the EH by minimizing the time variation of observables. 
The time evolution and the monitoring of observables is the quantum part, while the classical minimization of a cost function is the classical part.
For benchmarking, there are more cost functions implemented besides the QCFL.
See the [documentation](https://quantum-computing-software.gitlab.io/QuantVarEntHam/) for more information.


# Installation

Add the package via the package manager 
```julia
] add https://gitlab.com/quantum-computing-software/QuantVarEntHam
```
Test whether the installation was succesfull in the REPL
```julia
julia> using QuantVarEntHam

julia> z(4,(1,2))
Pauli string
Spin 1//2
Number of spins: 4

Z₁⊗ Z₂
```