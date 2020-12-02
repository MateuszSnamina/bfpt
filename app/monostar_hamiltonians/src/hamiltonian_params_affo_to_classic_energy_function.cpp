#include <monostar_hamiltonians/hamiltonian_params_affo_to_classic_energy_function.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_to_classic_energy_function.hpp>

#include <sstream>
#include <iomanip>

// #######################################################################
// ## hamiltonian_params_affo_to_classic_energy_function                ##
// #######################################################################

namespace monostar_hamiltonians {

monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ hamiltonian_params_affo_to_classic_energy_function_of_theta(
    HamiltonianParamsAffo params, double SS_average) {
    const AcosPlusBsinPlusCsqcosPlusZ fun1 =
        hamiltonian_params_fo_to_classic_energy_function(params.get_hamiltonian_params_fo());
    const AcosPlusBsinPlusCsqcosPlusZ fun2 =
        hamiltonian_params_fo_to_classic_energy_function(params.get_hamiltonian_params_ss_fo());
    return fun1 + (SS_average * fun2);
}

}  // end of namespace monostar_hamiltonians
