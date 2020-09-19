#pragma once

#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/enum_model_type.hpp>
#include <model_monostar/enum_run_type.hpp>
#include <model_monostar/enum_es_momentum_domain.hpp>

#include <bfpt_common/common_recipe_print_flags.hpp>

#include <string>

struct InterpretedProgramOptions {
  unsigned n_sites;
  unsigned n_pt;
  ModelType model_type;
  double J_classical;
  double J_quantum;
  RunType run_type;
  EsMomentumDomainVariant es_momentum_domain;
  bfpt_common::CommonRecipePrintFlags print_flags;
};

InterpretedProgramOptions interpret_program_options(const RawProgramOptions&);
