#pragma once

#include <so_hamiltonians/hamiltonian_params_af_fo.hpp>

#include <optional>

// #######################################################################
// ## get_orbital_theta                                                 ##
// #######################################################################

namespace so_hamiltonians {

double get_orbital_theta(
    const HamiltonianParamsAfFo& hamiltonian_af_fo_params,
    std::optional<double> user_defined_overrule,
    double average_ss);

}  //end of namespace so_hamiltonians
