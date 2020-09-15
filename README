# Bfpt – what it is

Bfpt (Brute Force Perturbation Theory) is a numeric solver of quantum Hamiltonian eigenproblem for discrete one dimensional periodic systems.

The solver perform exact denationalization (ED) in a subspace of the system Hilbert space. The subspace is constructed in a spirit of perturbation calculus. The construction starts with the 0-th order subspace filled exclusively with quantum states provded by the solver user. Then the next order subspaces are constructed iteratively: The (n+1)-th order subspace is the n-th order subspace extended by the Hamiltonian image of the n-th order subspace. The solver user defines the target subspace order. Having the desired subspace constructed, Hamiltonian restricted to the subspace is diagonalized using a spare matrices numeric method.

All calculations (both at the subspace generation stage and the diagonalization stage) are performed in the ‘inverse space’ -- in the ream of quantum states with given Bloch's theorem pseudo-momentum.

The archetype problem fitting the solver domain is one dimensional Heisenberg antiferromagnet. Eventhough the Heisenberg model is a lighthouse for development, the author original motivation is to address spin-orbital Kugel–Khomskii models (like the one described for copper fluoride KCuF3).

The solver is conceived to be general and easily extensible so to encore further experimenting. The project tries to achieve the objectives by involving easy to play modern C++ with balanced static and dynamic polymorphism, limited usage of external libraries focusing only on battle-tested ones:  Armadillo, BLAS, LAPACK, ARPACK libraries stack for linear algebra, and boost library for standard general purposes.

## Bfpt – what it is not

- The solver is not meant to be a competitor for Bethe ansatz based solutions.
- Bfpt is not meant to be the fastest numeric solver created.
- Bfpt is not meant to target continuous systems.
- Bfpt is not meant to target  two and more dimensional systems.

# Build

The project requires build tools (cmake3.14 and gcc9) and 3p libraries (armadillo9 boost1.71) to be installed. On Ubuntu 20.04:

```
apt install libboost-dev
apt install libboost-program-options-dev
apt install libarmadillo-dev
cd build_release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j 4
```

# Example

The project contains an executable (model_monostar) showing how to use the solver for 1D Heisenberg antiferromagnet.
