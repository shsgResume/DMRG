This package implements the Density Matrix Renormalization Group algorithm [See_this_paper](https://web.archive.org/web/20070721172908/http://hedrock.ps.uci.edu/dmrgpaper/dmrgpap.pdf) which is a fast eigenvalue and eigenvector solver. The method was applied to the Heisenberg XXX model, or sometimes known as the Ising model or a spin chain in Condensed Matter Physics.

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

where J is the interaction strength of the spin chain;
k is th elength of the system in number of lattice sites;
LancStates is the number of eigenvalues and eigenvectors for a fast sparse matrix solver called the Lanczos algorithm;
stateSize is the number of energy states to save in the intermediate steps of the DMRG algorithm;

Plots of the energy convergence, the time taken for the simulations, the interaction strength (or called the correlation between the different lattice sites as the distance is increased) and also the entanglement entropy is plotted as a final step of the script.

The entanglement entropy is also a method to measure the correlations between the different sites. The entropy used here is the Von Neumann entripy of entanglement [See_This_paper](https://arxiv.org/abs/hep-th/9303048). It is given by the formula 


$$S(\rho_A) = -Tr \left[\rho_{A} log \rho_A \right] = -Tr \left[\rho_{B} log \rho_B \right] = S(\rho_B) = - \Sigma_i \alpha^2_i log (\alpha^2_i)$$

