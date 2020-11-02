#pragma once

#include <koperator_trait/trait_koperator.hpp>

#include <kbasis/basis.hpp>

#include <armadillo>

#include <omp.h>

#include <chrono>

namespace bfpt_common {

template <typename KoperatorTraitT>
arma::sp_cx_mat
generate_operator_matrix(
    const typename KoperatorTraitT::KoperatorT& koperator,
    const typename KoperatorTraitT::BasisT& basis,
    const unsigned k_n,
    unsigned n_threads) {
    // *********** asserts ****************************************************************
    static_assert(koperator_trait::IsTraitKoperator<KoperatorTraitT>::value);
    static_assert(KoperatorTraitT::is_koperator_trait);
    // *********** prepare ****************************************************************
    std::vector<arma::sp_cx_mat> kn_operator_builder_matrix_all(n_threads);
    for (unsigned i = 0; i < n_threads; i++) {
        kn_operator_builder_matrix_all[i] = arma::sp_cx_mat(basis.size(), basis.size());
    }
    // *********** filling ****************************************************************
    //const auto tp_fill_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
#pragma omp parallel for num_threads(n_threads) schedule(guided)
    for (arma::uword ket_kstate_idx = 0; ket_kstate_idx < basis.size(); ket_kstate_idx++) {
        const auto tid = omp_get_thread_num();
        KoperatorTraitT::fill_kn_operator_builder_matrix_coll(koperator, basis, ket_kstate_idx, kn_operator_builder_matrix_all[tid], k_n);
    }
    // const auto tp_fill_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //std::cout << "fill took     : " << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_fill_2 - tp_fill_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
    // *********** reduction **************************************************************
    //const auto tp_reduce_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
    for (unsigned d = 1; d < n_threads; d *= 2) {
#pragma omp parallel num_threads(n_threads)
        {
            const auto tid = omp_get_thread_num();
            if (tid % (2 * d) == 0) {
                const auto idx1 = tid;
                const auto idx2 = tid + d;
                if (idx2 < n_threads) {
                    kn_operator_builder_matrix_all[idx1] += kn_operator_builder_matrix_all[idx2];
                }
            }
        }
    }
    //const auto tp_reduce_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //std::cout << "reduce took   : " << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_reduce_2 - tp_reduce_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
    // *********** add transpose **********************************************************
    //const auto tp_trans_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
    kn_operator_builder_matrix_all[0] += kn_operator_builder_matrix_all[0].t();
    //const auto tp_trans_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //std::cout << "transpose took: " << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_trans_2 - tp_trans_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
    // *********** return *****************************************************************
    return kn_operator_builder_matrix_all[0];
}

}  // namespace bfpt_common
