<h1>About</h1>
This package implements the Density Matrix Renormalization Group algorithm [1](https://web.archive.org/web/20070721172908/http://hedrock.ps.uci.edu/dmrgpaper/dmrgpap.pdf) which is a fast eigenvalue and eigenvector solver. The method was applied to the Heisenberg XXX model, or sometimes known as the Ising model or a spin chain in Condensed Matter Physics.

The DMRG algorithm iteratively increases the system length $N$ by finding the ground state of the system at each step, by solving a matrix of size $N \times N$ . The algorithm proceeds as follows:

```
  1. The Hamiltonian of the superblock of size N, defined as the new system length is obtained by a Kronecker product
  2. The Hamiltonian is solved to get the ground state of the system
  3. The reduced density matrix is computed for each half system
  4. The ground state energy is used to truncate the operators
  5. One site is added onto the system for each half block
  6. The operators (including the Hamiltonian) for the reduced system + 1 becomes the new half block
  7. The new superblock of size N + 2 is created and the operations are repeated
```
The simulation is self-contained, and users can simply change the variables in the file DMRG_functions.py:

```
J = 1.0
k = 300
LancStates = 6
stateSize = [10] 
```

where J is the interaction strength of the spin chain;
k is th elength of the system in number of lattice sites;
LancStates is the number of eigenvalues and eigenvectors for a fast sparse matrix solver called the Lanczos algorithm;
stateSize is the number of energy states to save in the intermediate steps of the DMRG algorithm;

Plots of the energy convergence, the time taken for the simulations, the interaction strength (or called the correlation between the different lattice sites as the distance is increased) and also the entanglement entropy is plotted as a final step of the script.

The entanglement entropy is also a method to measure the correlations between the different sites. The entropy used here is the Von Neumann entripy of entanglement [2](https://arxiv.org/abs/hep-th/9303048). It is given by the formula 


$$S(\rho_A) = -Tr \left[\rho_{A} log \rho_A \right] = -Tr \left[\rho_{B} log \rho_B \right] = S(\rho_B) = - \Sigma_i \alpha^2_i log (\alpha^2_i)$$

<h1>Requirements</h1>

> [!NOTE] 
> Python<br>
> Numpy<br>
> Matplotlib<br>
