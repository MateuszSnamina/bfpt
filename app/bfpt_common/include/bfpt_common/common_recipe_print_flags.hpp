#ifndef BFPT_COMMON_COMMON_RECIPIE_PRINT_FLAGS_HPP
#define BFPT_COMMON_COMMON_RECIPIE_PRINT_FLAGS_HPP

// #######################################################################
// ## CommonRecipePrintFlags                                            ##
// #######################################################################

namespace bfpt_common {

struct CommonRecipePrintFlags {
    bool print_unpopulated_basis_flag = false;
    bool print_unpopulated_basis_size_flag = true;
    bool print_populated_basis_flag = false;
    bool print_populated_basis_size_flag = true;
    bool print_sp_hamiltonian_flag = false;
    bool print_hamiltonian_flag = false;
    bool print_eigen_values_flag = true;
    bool print_eigen_vectors_flag = false;
};

}  // namespace bfpt_common

#endif