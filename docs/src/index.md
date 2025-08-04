# QuantVarEntHam.jl Documentation

# What? 

This package aims to find the Entanglement Hamiltonian (EH) of one dimensional spin lattice models with arbitry spin. 
Kokail et al. presented the algorithm -a Quantum Classical Feedback Loop (QCFL)- , which is the foundation for this package.
It learns the EH by minimizing the time evolution of observables. 
The time evolution and the monitoring of observables is the quantum part, while the classical minimization of a cost function is the classical part.
For benchmarking, there are more cost functions implemented besides the [`QCFL`](@ref) (see [`Commutator`](@ref), [`RelativeEntropy`](@ref)).


# Theory in short

Given a composite Hilbert space ``\mathscr{H} = \mathscr{H}_\text{A} \otimes \mathscr{H}_\text{B}``, a general quantum state 
after a Schmidt decomposition can be written as 
```math 
| \Psi \rangle  = \sum_{\alpha = 1}^{\chi} e^{-\xi_{\alpha}/2} | \Phi^{\alpha}_\text{A} \rangle  \otimes | \Phi^{\alpha}_\text{B} \rangle,
```
where ``\chi`` is called the Schmidt rank.
Tracing out the degrees of freedom related to subsystem B in the quantum state 
``\rho = | \Psi \rangle \langle \Psi |``, 
the reduced density matrix on subsystem A reads 
```math 
    \rho_\text{A} = \text{Tr}_\text{B} \left [ \rho \right ] = \sum_{\alpha} e^{-\xi_{\alpha}} | \Phi_\text{A}^{\alpha} \rangle \langle  \Phi_\text{A}^{\alpha} | 
   =\mathrel{\mathop:} e^{-H_\text{A}},
``` 
which defines the EH ``H_\text{A}``on subsystem A.
Its spectrum - the Entanglemen spectrum - ``\{ \xi_{\alpha} \}``is the target quantity of this algorithm and gives insights about topological phases.
 