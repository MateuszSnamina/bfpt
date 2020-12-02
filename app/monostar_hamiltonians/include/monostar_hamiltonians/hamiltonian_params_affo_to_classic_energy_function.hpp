#pragma once

#include <monostar_hamiltonians/hamiltonian_params_affo.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>

#include <sstream>
#include <iomanip>

// #######################################################################
// ## hamiltonian_params_affo_to_classic_energy_function                ##
// #######################################################################

namespace monostar_hamiltonians {

AcosPlusBsinPlusCsqcosPlusZ
hamiltonian_params_affo_to_classic_energy_function_of_theta(
    HamiltonianParamsAffo params,
    double SS_average);

}  // end of namespace monostar_hamiltonians
