#include <so_hamiltonians/hamiltonian_params_af_fo_to_classic_energy_function.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo_helpers.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_to_classic_energy_function.hpp>

#include <sstream>
#include <iomanip>

// #######################################################################
// ## hamiltonian_params_af_fo_to_classic_energy_function               ##
// #######################################################################

namespace so_hamiltonians {

monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ hamiltonian_params_af_fo_to_classic_energy_function(
        HamiltonianParamsAfFo params, double SS_average) {
    const monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ fun1 =
            monostar_hamiltonians::hamiltonian_params_fo_to_classic_energy_function(params.get_hamiltonian_params_fo());
    const monostar_hamiltonians::AcosPlusBsinPlusCsqcosPlusZ fun2 =
            monostar_hamiltonians::hamiltonian_params_fo_to_classic_energy_function(params.get_hamiltonian_params_ss_fo());
    return fun1 + (SS_average * fun2);//TODO finish
}

}  // end of namespace so_hamiltonians
