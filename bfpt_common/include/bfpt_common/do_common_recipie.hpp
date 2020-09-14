#ifndef BFPT_COMMON_DO_COMMON_RECIPIE_HPP
#define BFPT_COMMON_DO_COMMON_RECIPIE_HPP

#include <bfpt_common/i_kstate_hamiltonian.hpp>
#include <bfpt_common/populate_pt_basis.hpp>
#include <bfpt_common/common_recipe_print_flags.hpp>

#include <linear_algebra/linear_algebra.hpp>

#include <kstate/basis.hpp>
#include <kstate/kstate_concrete.hpp>

#include <armadillo>

#include <iostream>
#include <string>

#include <cassert>

// #######################################################################
// ## do_common_recipe                                                  ##
// #######################################################################

namespace bfpt_common {

//template <typename Element> //TODO remove!
template<typename KstateT>
inline double do_common_recipe(const IKstatePopulator<KstateT>& bais_populator,
                               const IKstateHamiltonian<KstateT>& hamiltonian,
                               kstate::Basis<KstateT>& basis,
                               const unsigned max_pt_order, const unsigned k_n,
                               CommonRecipePrintFlags print_flags) {
    arma::wall_clock timer;
    __attribute__((unused)) const size_t n_sites = basis.n_sites();
    const std::string message_prefix = "[common-recipe] ";
    const std::string progress_tag = "[progress] ";
    const std::string data_tag = "[data    ] ";
    const std::string time_tag = "[time    ] ";
    assert(k_n < n_sites);
    // --------------------------------------------------
    if (print_flags.print_unpopulated_basis_flag) {
        std::cout << message_prefix << data_tag << "Unpopulated basis (0'th pt-order basis):";
        std::cout << basis;
    }
    if (print_flags.print_unpopulated_basis_size_flag) {
        std::cout << message_prefix << data_tag << "Unpopulated basis (0'th pt-order basis) size: "
                  << basis.size() << "." << std::endl;
    }
    // --------------------------------------------------
    // Define hamiltonian and basis:
    std::cout << message_prefix << progress_tag << "About to populate pt-basis." << std::endl;
    // Generate higher pt-orders subspace basis:
    timer.tic();
    populate_pt_basis(bais_populator, max_pt_order, basis);
    const double time_populating_pt_basis = timer.toc();
    std::cout << message_prefix << time_tag << "Populating pt-basis took " << time_populating_pt_basis << "s." << std::endl;
    std::cout << message_prefix << progress_tag << "Has populated pt-basis." << std::endl;
    // --------------------------------------------------
    if (print_flags.print_populated_basis_flag) {
        std::cout << message_prefix << data_tag << "Populated basis (" << max_pt_order << "'th pt-order basis):";
        std::cout << basis;
    }
    if (print_flags.print_populated_basis_size_flag) {
        std::cout << message_prefix << data_tag << "Populated basis (" << max_pt_order << "'th pt-order basis) size: "
                  << basis.size() << "." << std::endl;
    }
    // --------------------------------------------------
    // Generate hamiltoniam matrix:
    std::cout << message_prefix << progress_tag << "About to generate hamiltoniam." << std::endl;
    timer.tic();
    const auto kn_hamiltonian_matrix = hamiltonian.make_kn_hamiltonian_matrix(basis, k_n);
    const double time_generating_kn_hamiltonian_matrix = timer.toc();
    std::cout << message_prefix << time_tag << "Generating kn-hamiltoniam matrix took " << time_generating_kn_hamiltonian_matrix << "s." << std::endl;
    std::cout << message_prefix << progress_tag << "Has generated kn-hamiltoniam." << std::endl;
    // --------------------------------------------------
    if (print_flags.print_sp_hamiltonian_flag) {
        std::cout << message_prefix << data_tag << "kn_hamiltonian_matrix:";
        std::cout << kn_hamiltonian_matrix;
    }
    if (print_flags.print_hamiltonian_flag) {
        std::cout << message_prefix << data_tag << "kn_hamiltonian_matrix:";
        std::cout << arma::cx_mat(kn_hamiltonian_matrix);
    }
    // --------------------------------------------------
    std::cout << message_prefix << progress_tag << "About to solve eigen problem." << std::endl;
    timer.tic();
    const auto& eigs_sym_result = lin_alg::fallbacked_eigs_sym(kn_hamiltonian_matrix, 1, 1e-6);
    const double time_solving_eigen_problem = timer.toc();
    if (eigs_sym_result.is_err()) {
        std::cout << message_prefix << time_tag << "Solving eigen problem took: " << time_solving_eigen_problem << "s." << std::endl;
        std::cout << message_prefix << progress_tag << "Failed to solve eigen problem." << std::endl;
        std::cout << message_prefix << progress_tag << "The reported error message:" << std::endl;
        std::cout << eigs_sym_result.unwrap_err().what() << std::endl;
        return arma::datum::nan;
    }
    const auto eigen_info = eigs_sym_result.unwrap();
    const arma::vec& eigen_values = eigen_info.eigen_values;
    const arma::cx_mat& eigen_vectors = eigen_info.eigen_vectors;
    std::cout << message_prefix << time_tag << "Solving eigen problem took: " << time_solving_eigen_problem << "s." << std::endl;
    std::cout << message_prefix << progress_tag << "Has solved eigen problem." << std::endl;
    // --------------------------------------------------
    if (print_flags.print_eigen_values_flag) {
        std::cout << message_prefix << data_tag << "eigen_values:" << std::endl;
        std::cout << message_prefix << eigen_values;
        // std::cout << "min eigen_value: " << eigen_values(0) << std::endl;
    }
    if (print_flags.print_eigen_vectors_flag) {
        std::cout << message_prefix << data_tag << "eigen_vectors:" << std::endl;
        std::cout << message_prefix << eigen_vectors;
        // std::cout << "eigen_vectors.col(0): " << std::endl
        //           << eigen_vectors.col(0) / eigen_vectors(0, 0) << std::endl;
    }
    // --------------------------------------------------
    return eigen_values(0);
}

}  // namespace bfpt_common

#endif
