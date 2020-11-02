#pragma once

#include <utility>

// #######################################################################
// ## CommonRecipePrintFlags                                            ##
// #######################################################################

namespace bfpt_common {

struct CommonRecipePrintFlags {
    bool print_unpopulated_basis_flag = false;
    bool print_unpopulated_basis_size_flag = true;
    bool print_populated_basis_flag = false;
    bool print_populated_basis_size_flag = true;
    bool print_hamiltonian_stats = true;
    bool print_sp_hamiltonian_flag = false;
    bool print_hamiltonian_flag = false;
    bool print_eigen_values_flag = true;
    bool print_eigen_vectors_flag = false;
    bool print_pretty_vectors_flag = false;
    bool print_density_operator_flag = false;
    std::pair<unsigned, unsigned> print_pretty_min_max_n_kstates{10, 50};
    double print_pretty_probability_treshold = 0.05;
};

}  // namespace bfpt_common
