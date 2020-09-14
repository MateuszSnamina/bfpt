#pragma once

#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/enum_model_type.hpp>

#include <string>

struct InterpretedProgramOptions {
  unsigned n_sites;
  unsigned n_pt;
  ModelType model_type;
  double J_classical;
  double J_quantum;
};

InterpretedProgramOptions interpret_program_options(const RawProgramOptions&);
