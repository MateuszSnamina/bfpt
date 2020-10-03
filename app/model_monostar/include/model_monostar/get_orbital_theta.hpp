#pragma once

#include<model_monostar/hamiltonian_params_fo.hpp>

#include<optional>

// #######################################################################
// ## get_orbital_theta                                                 ##
// #######################################################################

double get_orbital_theta(
        const HamiltonianParamsFo& hamiltonian_fo_params,
        std::optional<double> user_defined_overrule);
