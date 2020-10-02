#pragma once

#include<model_monostar/hamiltonian_fo_params.hpp>

#include<optional>

// #######################################################################
// ## get_orbital_theta                                                 ##
// #######################################################################

double get_orbital_theta(
        const HamiltonianFoParams& hamiltonian_fo_params,
        std::optional<double> user_defined_overrule);
