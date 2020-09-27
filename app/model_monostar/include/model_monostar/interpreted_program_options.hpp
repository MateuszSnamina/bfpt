#pragma once

#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/enum_model_type.hpp>
#include <model_monostar/enum_run_type.hpp>
#include <model_monostar/enum_es_momentum_domain.hpp>
#include <model_monostar/hamiltonian_params.hpp>

#include <bfpt_common/common_recipe_print_flags.hpp>

#include <string>

struct InterpretedProgramOptions {
  unsigned n_sites;
  unsigned n_pt;
  ModelType model_type;
  HamiltonianParams hamiltonian_params;
  RunType run_type;
  EsMomentumDomainVariant es_momentum_domain;
  bfpt_common::CommonRecipePrintFlags print_flags;
  unsigned n_threads;
};

InterpretedProgramOptions interpret_program_options(const RawProgramOptions&);
