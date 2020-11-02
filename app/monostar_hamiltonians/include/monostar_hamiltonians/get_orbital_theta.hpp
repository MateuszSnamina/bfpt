#pragma once

#include <monostar_hamiltonians/hamiltonian_params_fo.hpp>

#include <optional>

// #######################################################################
// ## get_orbital_theta                                                 ##
// #######################################################################

namespace monostar_hamiltonians {

double get_orbital_theta(
    const HamiltonianParamsFo& hamiltonian_fo_params,
    std::optional<double> user_defined_overrule);

}  //end of namespace monostar_hamiltonians
