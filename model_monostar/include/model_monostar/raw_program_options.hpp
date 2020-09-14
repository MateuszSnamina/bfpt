#pragma once

#include <string>

struct RawProgramOptions {
  unsigned n_sites;
  unsigned n_pt;
  std::string model_type_string;
  double J_classical;
  double J_quantum;
};

RawProgramOptions grep_program_options(int argc, char** argv);
