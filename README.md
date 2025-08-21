# QuantVarEntHam

[![Build Status](https://gitlab.dlr.de/ma-kind/QuantVarEntHam/-/badges/release.svg)](https://gitlab.dlr.de/ma-kind/QuantVarEntHam)
[![pipeline status](https://gitlab.dlr.de/ma-kind/QuantVarEntHam/badges/main/pipeline.svg)](https://gitlab.dlr.de/ma-kind/QuantVarEntHam/-/pipelines)
[![coverage report](https://gitlab.dlr.de/ma-kind/QuantVarEntHam/badges/main/coverage.svg)](https://gitlab.dlr.de/ma-kind/QuantVarEntHam/-/graphs/main/charts)
[![documentation (placeholder)](https://img.shields.io/badge/docs-latest-blue.svg)](https://quantvarentham-ma-kind-28334cbd97e4f53ca6d853b2420a209527d76ff3.pages.gitlab.dlr.de/)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)


This package aims to find the Entanglement Hamiltonian (EH) of one dimensional spin lattice models with arbitrary spin. 
Kokail et al. presented the algorithm -a Quantum Classical Feedback Loop (QCFL)- , which is the foundation of this package.
It learns the EH by minimizing the time variation of observables. 
The time evolution and the monitoring of observables is the quantum part, while the classical minimization of a cost function is the classical part.
For benchmarking, there are more cost functions implemented besides the QCFL.


# Installation

Clone the repo 
```
git clone git@gitlab.dlr.de:ma-kind/QuantVarEntHam.git
```
and navigate into the folder `QuantVarEntHam`.
Add the package via the package manager 
```julia
] add .
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