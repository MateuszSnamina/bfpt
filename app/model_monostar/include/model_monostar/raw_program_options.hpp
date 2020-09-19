#pragma once

#include <string>

struct RawProgramOptions {
  unsigned n_sites;
  unsigned n_pt;
  std::string model_type_string;
  double J_classical;
  double J_quantum;
  std::string run_type_string;
  bool print_unpopulated_basis_flag;
  bool print_unpopulated_basis_size_flag;
  bool print_populated_basis_flag;
  bool print_populated_basis_size_flag;
  bool print_sp_hamiltonian_flag;
  bool print_hamiltonian_flag;
  bool print_eigen_values_flag;
  bool print_eigen_vectors_flag;
  std::string es_momentum_domain_string;
  unsigned es_n_k;
};

RawProgramOptions grep_program_options(int argc, char** argv);
