#pragma once

#include <monostar_hamiltonians/hamiltonian_params_fo.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>

// #######################################################################
// ## hamiltonian_params_fo_to_classic_energy_function                  ##
// #######################################################################

namespace monostar_hamiltonians {

monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ hamiltonian_params_fo_to_classic_energy_function(HamiltonianParamsFo params);

}  // end of namespace monostar_hamiltonians
