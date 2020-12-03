#pragma once

#include<monostar_hamiltonians/hamiltonian_params_agile_affo.hpp>
#include<monostar_hamiltonians/hamiltonian_params_jkl01.hpp>

namespace monostar_hamiltonians {

HamiltonianParamsJkl01 dacay_hamiltonian_params_agile_affo(
        HamiltonianParamsAgileAffo hamiltonian_params_agile_affo,
        double orbital_theta);

}  // end of namespace monostar_hamiltonians
