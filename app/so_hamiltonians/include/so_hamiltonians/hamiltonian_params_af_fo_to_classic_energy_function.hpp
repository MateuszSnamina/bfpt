#pragma once

#include <so_hamiltonians/hamiltonian_params_af_fo.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>

#include <sstream>
#include <iomanip>

// #######################################################################
// ## hamiltonian_params_af_fo_to_classic_energy_function               ##
// #######################################################################

namespace so_hamiltonians {

monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ
hamiltonian_params_af_fo_to_classic_energy_function_of_theta(
        HamiltonianParamsAfFo params,
        double SS_average);

}  // end of namespace so_hamiltonians
