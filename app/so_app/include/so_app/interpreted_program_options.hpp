#pragma once

#include <so_app/raw_program_options.hpp>
#include <so_app/enum_model_type.hpp>
#include <so_app/enum_run_type.hpp>
#include <so_app/enum_es_momentum_domain.hpp>

#include <monostar_hamiltonians/hamiltonian_params_af_fm.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo.hpp>

#include <bfpt_common/common_recipe_print_flags.hpp>

#include <string>
#include <optional>

namespace so_app {

struct InterpretedProgramOptions {
    unsigned n_sites;
    unsigned n_pt;
    ModelType model_type;
    monostar_hamiltonians::HamiltonianParamsAfFm hamiltonian_params_af_fm;
    monostar_hamiltonians::HamiltonianParamsFo hamiltonian_params_fo;
    std::optional<double> orbital_theta;
    RunType run_type;
    EsMomentumDomainVariant es_momentum_domain;
    bfpt_common::CommonRecipePrintFlags print_flags;
    unsigned n_threads;
};

InterpretedProgramOptions interpret_program_options(const RawProgramOptions&);

}  // end of namespace so_app
