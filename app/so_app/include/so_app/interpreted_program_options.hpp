#pragma once

#include <so_app/raw_program_options.hpp>
#include <so_app/enum_model_type.hpp>
#include <so_app/enum_run_type.hpp>
#include <so_app/enum_es_momentum_domain.hpp>

#include <so_hamiltonians/hamiltonian_params_af_fo.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo.hpp>

#include <bfpt_common/common_recipe_print_flags.hpp>

#include <string>
#include <optional>

namespace so_app {

struct InterpretedProgramOptions {
    unsigned n_sites;
    unsigned n_pt;
    std::optional<unsigned> n_max_site_spin_excitations;
    std::optional<unsigned> n_max_site_orbit_excitations;
    bool accept_orbit_site_excitations_only_if_near_domain_wall;
    ModelType model_type;
    so_hamiltonians::HamiltonianParamsAfFo hamiltonian_params_af_fo = so_hamiltonians::HamiltonianParamsAfFo::Builder().build();
    std::optional<double> orbital_theta;
    double average_ss;
    RunType run_type;
    EsMomentumDomainVariant es_momentum_domain;
    bfpt_common::CommonRecipePrintFlags print_flags;
    unsigned n_threads;
};

InterpretedProgramOptions interpret_program_options(const RawProgramOptions&);

}  // end of namespace so_app
