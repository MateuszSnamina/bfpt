# Bfpt

## What it is
Bfpt (Brute Force Perturbation Theory) is a numeric solver of quantum Hamiltonian eigenproblem for finite discrete one dimensional periodic systems.

The solver performs exact denationalization (ED) in a subspace of the system Hilbert space. The subspace is constructed in a spirit of perturbation calculus. The construction starts with the `0`-th order subspace filled exclusively with quantum states provded (as an input) by the solver user. Then the next order subspaces are constructed iteratively: The `(n+1)`-th order subspace is the `n`-th order subspace extended by the Hamiltonian image of the `n`-th order subspace. The solver user defines the target subspace order. Having the desired subspace constructed, Hamiltonian restricted to the subspace is diagonalized using a spare matrices numeric method.

Note that for sufficiently high value of the order parameter the relevant subspace is an invariant subspace of Hamiltonian and
the solver performs ED.

All calculations (both at the subspace generation stage and the diagonalization stage) are performed in the 'inverse space' (or 'momentum space') -- in the realm of quantum states with given Bloch's theorem pseudo-momentum.

The archetype problem fitting the solver domain is one dimensional Heisenberg antiferromagnet. Even though the Heisenberg model is a lighthouse for development, the author's original motivation is to address spin-orbital Kugel–Khomskii models (like the one described for copper fluoride KCuF3).

The solver is conceived to be general and easily extensible so to encore further experimenting. The project tries to achieve the objectives by involving easy to play modern C++ with balanced static and dynamic polymorphism, limited usage of external libraries focusing only on battle-tested ones: `boost C++` libraries for general purposes and `Armadillo`, `LAPACK`, `ARPACK`, `BLAS` libraries stack for linear algebra.

## What it is not

- The solver is not meant to be a competitor for Bethe ansatz based solutions.
- Bfpt is not meant to be the fastest numeric solver created.
- Bfpt is not meant to target continuous systems.
- Bfpt is not meant to target  two and more dimensional systems.

# Compile/build/install

Compiling/linking the project requires basic dev tools (`cmake3.16` and `gcc9`) and 3p libraries (`armadillo9`, `boost1.71`) to be installed.
On Ubuntu 20.04:

```
apt install gcc-9
apt install cmake
apt install libboost-dev
apt install libboost-program-options-dev
apt install libarmadillo-dev
mkdir build_release
pushd build_release
cmake -DCMAKE_BUILD_TYPE=Release ..
make "-j$(nproc)"
popd
build_release/bin/model_monostar # try it!
```

# Simple example -- model monostar

The project contains an executable (`model_monostar`) showing
how to use the solver for models like 1D Heisenberg (anti)ferromagnet.
 
## Model defintion

**Hilbert space** Monostar model is defined for quantum systems made of replicas of a two-level subsystem.
The replicas are arranged in a way they form a chain (and here will be called chain nodes); 
The nodes interactions are assumed to be restricted to the nearest neighbors.

**Interactions pattern** A chain node quantum levels are denoted: `gs` and `es` (like ground state and excited state).
The system Hamiltonian `H` is translationally invariant (with periodic boundary conditions)
and is prescribed by two-nodes kernel Hamiltionian `H_12` (`H = \sum_{<ij>} H_12(i,j)` with the summation over pairs of adjacent nodes).
`H_12` is further decomposed into `H_12_diag` and `H_12_off_diag` parts;
`H_12_diag` defines diagonal energies for `(gs, gs)`, `(gs, es)`, `(es, gs)`, `(es, es)` nodes pairs,
whereas `H_12_off_diag` describes couplings between the pairs.

**Considered Hamiltonians** Model monostar is a framework that 

TODO

## Model monostat with `af` and `fm` Hamiltonians

The model goes in two flavors: `fm` and `af`, each parameterized with two real values `J_classical` and `J_quantum`.
The `H_12_diag` does not depend on the `fm`/`af` variant and is governed by `J_classical`:
| states pair | energy contribution |
|-------------|---------------------|
| `(gs, gs)`  | `-J_classical/4`    |
| `(gs, es)`  | `+J_classical/4`    |
| `(es, gs)`  | `+J_classical/4`    |
| `(es, es)`  | `-J_classical/4`    |

`H_12_off_diag` in `fm` model variant is given by:
| coupled states pairs    | coupling contribution |
|-------------------------|-----------------------|
| `(es, gs)` ⟷ `(gs, es)` | `-J_quantum/2`        |

`H_12_off_diag` in `af` model variant is given by:
| coupled states pairs    | coupling contribution |
|-------------------------|-----------------------|
| `(gs, gs)` ⟷ `(es, es)` | `+J_quantum/2`        |

`fm` monostar model is trivially equivalent to Heisenberg ferromagnet model with `gs` translated into `spin down` and `es` translated into `spin up`. `af` monostar model is equivalent to Heisenberg antiferromagnet with the monostar states to spin states association given by `gs`≡`down`, `es`≡`up` for nodes on one magnetic sub-lattice, and `gs`≡`up`, `es`≡`down` on the other.
To keep the equivalence explicit, here term "monostar model" is used instead of "Heisenberg (anti)ferromagnet mode" to avoid confusion.

## Results from bftp solver

The plot below presents the system Hamiltonian eigenenergies for `af` monostar system of 20 nodes calculated with `bfpt` at different levels of approximation. The considered system states are: the system ground state (horizontal lines on the plot) and system states from the first excited band (the curves on the plot). Result obtained in calculations with the Hilbert space subspaces of orders from 1 to 8 (inclusive) are represented by lines with colors going from red to blue. In addition exact results for infinite chain was included for reference; the reference energies are represented by gray lines. (The ground state energy was scaled so to preserve the quantum correlation energy per node.)

!["Monostar model -- absolute energies"](img/model_monostar_20sites_absolute_energy.png "Monostar model -- absolute energies")

The model excitation energies are obtained as a difference of the relevant excited states energy and the ground state energy.

!["Monostar model -- excitation energies"](img/model_monostar_20sites_excitation_enery.png "Monostar model -- excitation energies")

The content of the eigenstates may be printed in a pretty form in which the state configuration contributions for are sorted and displayed in human friendly way. The screenshot below shows the decomposition of the groundstate of 16 nodes `af` monostar model.

!["Monostar model -- excitation energies"](img/monostar_16sites_gs_eigenvector.png )

## Executable CLI synopsis

[TODO]

# Example -- a spinorbital
