This package implements the Density Matrix Renormalization Group algorithm which is a fast eigenvalue and eigenvector solver. The method was applied to the Heisenberg XXX model, or sometimes known as the Ising model or a spin chain in Condensed Matter Physics.

The algorithm proceeds as follows:

```
  1.
  2.
  3.
  4.
  5.
```
The simulation is self-contained, and users can simply change the variables in the file DMRG_functions.py:

```
J = 1.0k = 300
LancStates = 6
stateSize = [10] 
```
where J is the interaction strength of the spin chain
k is th elength of the system in number of lattice sites
LancStates is the number of eigenvalues and eigenvectors for a fast sparse matrix solver called the Lanczos algorithm
stateSize is the number of energy states to save in the intermediate steps of the DMRG algorithm

Plots of the energy convergence, the time taken for the simulations, the interaction strength (or called the correlation between the different lattice sites as the distance is increased) and also the entanglement entropy is plotted as a final step of the script.
