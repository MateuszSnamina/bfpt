#pragma once

#include <kbasis/basis.hpp>

#include <kstate_view_amend_spec/amend_spec.hpp>

#include <kstate_trait/trait_kstate.hpp>

#include <omp.h>

#include <armadillo>
#include <map>
#include <complex>
#include <vector>

// #######################################################################
// ## DensityOperator12                                                 ##
// #######################################################################

namespace bfpt_common {

template<typename KstateTraitT>
void
calculate_reduced_density_operator_12_impl(
        const kbasis::Basis<KstateTraitT>& basis,
        const arma::cx_vec& eigen_vector,
        const arma::uword ket_kstate_idx,
        arma::cx_mat& result_accumulator) {
    // static asserts:
    static_assert(kstate_trait::IsTraitKstate<KstateTraitT>::value);
    static_assert(KstateTraitT::is_kstate_trait);
    // helper types:
    using KstateT = typename KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    // helper values:
    const auto n_sites = basis.n_sites();
    const auto site_basis_dim = SiteStateTraitT::site_basis_dim();
    // [[expect]]:
    assert(result_accumulator.n_rows == site_basis_dim * site_basis_dim);
    assert(result_accumulator.n_cols == site_basis_dim * site_basis_dim);
    assert(eigen_vector.n_rows == basis.size());
    // the job:
    const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
    assert(ket_kstate_ptr);
    const auto ket_kstate_view = KstateTraitT::to_view(*ket_kstate_ptr);
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % n_sites) {
        const auto ket_kernel_site_1 = KstateTraitT::view_n_th_site_state(ket_kstate_view, n_delta);
        const auto ket_kernel_site_2 = KstateTraitT::view_n_th_site_state(ket_kstate_view, n_delta_p1);
        const auto ket_kernel_site_1_idx = SiteStateTraitT::get_index(ket_kernel_site_1);
        const auto ket_kernel_site_2_idx = SiteStateTraitT::get_index(ket_kernel_site_2);
        for (unsigned bra_kernel_site_1_idx = 0; bra_kernel_site_1_idx < site_basis_dim; bra_kernel_site_1_idx++) {
            for (unsigned bra_kernel_site_2_idx = 0; bra_kernel_site_2_idx < site_basis_dim; bra_kernel_site_2_idx++) {
                const auto& bra_kernel_site_1 = SiteStateTraitT::from_index(bra_kernel_site_1_idx);
                const auto& bra_kernel_site_2 = SiteStateTraitT::from_index(bra_kernel_site_2_idx);
                const auto refined_holder_1 = kstate_view_amend_spec::refined(n_delta, bra_kernel_site_1);// Must outlive bra_kstate_view
                const auto refined_holder_2 = kstate_view_amend_spec::refined(n_delta_p1, bra_kernel_site_2);// Must outlive bra_kstate_view.
                const auto bra_kstate_view_preproduct = KstateTraitT::refined_view(ket_kstate_view, refined_holder_1);// Must outlive bra_kstate_view
                const auto bra_kstate_view = KstateTraitT::refined_view(bra_kstate_view_preproduct, refined_holder_2);
                const size_t bra_kstate_n_unique_shift = KstateTraitT::view_n_unique_shift(bra_kstate_view);
                const auto roration_spec = kstate_view_amend_spec::rotated(bra_kstate_n_unique_shift);
                const auto bra_kstate_view_unique_shifted = KstateTraitT::rotated_view(bra_kstate_view, roration_spec); // equivalent to `kstate::make_unique_shift(bra_kstate)`
                if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(bra_kstate_view_unique_shifted)) {
                    const auto bra_kstate_idx = *bra_kstate_optional_idx;
                    const double pre_norm_1 = KstateTraitT::norm_factor(*basis.vec_index()[bra_kstate_idx]) * KstateTraitT::norm_factor(*basis.vec_index()[ket_kstate_idx]);
                    const size_t bra_n_least_replication_shift = KstateTraitT::n_least_replication_shift(*basis.vec_index()[bra_kstate_idx]);
                    const size_t bra_n_replicas = n_sites / bra_n_least_replication_shift;
                    const std::complex<double> value = std::conj(eigen_vector(bra_kstate_idx)) * eigen_vector(ket_kstate_idx);
                    const auto ket_kstate_matrix_idx = site_basis_dim * ket_kernel_site_1_idx + ket_kernel_site_2_idx;
                    const auto bra_kstate_matrix_idx = site_basis_dim * bra_kernel_site_1_idx + bra_kernel_site_2_idx;
                    result_accumulator(bra_kstate_matrix_idx, ket_kstate_matrix_idx) += pre_norm_1 * bra_n_replicas * value;
                }
            } // end of bra_kernel_site_1_idx loop
        } // end of bra_kernel_site_1_idx loop
    } // end of n_delta loop
}

template<typename KstateTraitT>
arma::cx_mat
calculate_reduced_density_operator_12(
        const kbasis::Basis<KstateTraitT>& basis,
        const arma::cx_vec& eigen_vector,
        unsigned n_threads) {
    // *********** asserts ****************************************************************
    static_assert(kstate_trait::IsTraitKstate<KstateTraitT>::value);
    static_assert(KstateTraitT::is_kstate_trait);
    // *********** helper types ***********************************************************
    using KstateT = typename KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    // *********** helper values **********************************************************
    const auto site_basis_dim = SiteStateTraitT::site_basis_dim();
    // *********** [[expect]] *************************************************************
    assert(eigen_vector.n_rows == basis.size());
    // *********** prepare ****************************************************************
    std::vector<arma::cx_mat> density_matrix_all(n_threads);
    for (unsigned i = 0; i < n_threads; i++) {
        density_matrix_all[i] = arma::cx_mat(site_basis_dim * site_basis_dim, site_basis_dim * site_basis_dim, arma::fill::zeros);
    }
    // *********** filling ****************************************************************
    //const auto tp_fill_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
#pragma omp parallel for num_threads(n_threads) schedule(guided)
    for (arma::uword ket_kstate_idx = 0; ket_kstate_idx < eigen_vector.n_rows; ket_kstate_idx++) {
        const auto tid = omp_get_thread_num();
        calculate_reduced_density_operator_12_impl(basis, eigen_vector, ket_kstate_idx, density_matrix_all[tid]);
    }
    // const auto tp_fill_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //std::cout << "fill took     : " << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_fill_2 - tp_fill_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
    // *********** reduction **************************************************************
    //const auto tp_reduce_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
    for (unsigned d = 1; d < n_threads; d *=2) {
#pragma omp parallel num_threads(n_threads)
        {
            const auto tid = omp_get_thread_num();
            if (tid % (2 * d) == 0) {
                const auto idx1 = tid;
                const auto idx2 = tid + d;
                if (idx2 < n_threads) {
                    density_matrix_all[idx1] += density_matrix_all[idx2];
                }
            }
        }
    }
    //const auto tp_reduce_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //std::cout << "reduce took   : " << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_reduce_2 - tp_reduce_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
    // *********** return *****************************************************************
    return density_matrix_all[0];
}

} // end of namespace bfpt_common;

// #######################################################################
// ## DensityOperator1                                                  ##
// #######################################################################

namespace bfpt_common {

template<typename KstateTraitT>
void
calculate_reduced_density_operator_1_impl(
        const kbasis::Basis<KstateTraitT>& basis,
        const arma::cx_vec& eigen_vector,
        const arma::uword ket_kstate_idx,
        arma::cx_mat& result_accumulator) {
    // static asserts:
    static_assert(kstate_trait::IsTraitKstate<KstateTraitT>::value);
    static_assert(KstateTraitT::is_kstate_trait);
    // helper types:
    using KstateT = typename KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    // helper values:
    const auto n_sites = basis.n_sites();
    const auto site_basis_dim = SiteStateTraitT::site_basis_dim();
    // [[expect]]:
    assert(result_accumulator.n_rows == site_basis_dim);
    assert(result_accumulator.n_cols == site_basis_dim);
    assert(eigen_vector.n_rows == basis.size());
    // the job:
    const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
    assert(ket_kstate_ptr);
    const auto ket_kstate_view = KstateTraitT::to_view(*ket_kstate_ptr);
    for (size_t n_delta = 0; n_delta < n_sites; n_delta++) {
        const auto ket_kernel_site_1 = KstateTraitT::view_n_th_site_state(ket_kstate_view, n_delta);
        const auto ket_kernel_site_1_idx = SiteStateTraitT::get_index(ket_kernel_site_1);
        for (unsigned bra_kernel_site_1_idx = 0; bra_kernel_site_1_idx < site_basis_dim; bra_kernel_site_1_idx++) {
            const auto bra_kernel_site_1 = SiteStateTraitT::from_index(bra_kernel_site_1_idx);
            const auto refined_holder_1 = kstate_view_amend_spec::refined(n_delta, bra_kernel_site_1); // Must outlive bra_kstate_view.
            const auto bra_kstate_view = KstateTraitT::refined_view(ket_kstate_view, refined_holder_1);
            const size_t bra_kstate_n_unique_shift = KstateTraitT::view_n_unique_shift(bra_kstate_view);
            const auto roration_spec = kstate_view_amend_spec::rotated(bra_kstate_n_unique_shift);
            const auto bra_kstate_view_unique_shifted = KstateTraitT::rotated_view(bra_kstate_view, roration_spec); // equivalent to `kstate::make_unique_shift(bra_kstate)`
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(bra_kstate_view_unique_shifted)) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                const double pre_norm_1 = KstateTraitT::norm_factor(*basis.vec_index()[bra_kstate_idx]) * KstateTraitT::norm_factor(*basis.vec_index()[ket_kstate_idx]);
                const size_t bra_n_least_replication_shift = KstateTraitT::n_least_replication_shift(*basis.vec_index()[bra_kstate_idx]);
                const size_t bra_n_replicas = n_sites / bra_n_least_replication_shift;
                const std::complex<double> value = std::conj(eigen_vector(bra_kstate_idx)) * eigen_vector(ket_kstate_idx);
                result_accumulator(ket_kernel_site_1_idx, bra_kernel_site_1_idx) += pre_norm_1 * bra_n_replicas * value;
            }
        } // end of bra_kernel_site_1_idx loop
    } // end of n_delta loop
}

template<typename KstateTraitT>
arma::cx_mat
calculate_reduced_density_operator_1(
        const kbasis::Basis<KstateTraitT>& basis,
        const arma::cx_vec& eigen_vector,
        unsigned n_threads) {
    // *********** asserts ****************************************************************
    static_assert(kstate_trait::IsTraitKstate<KstateTraitT>::value);
    static_assert(KstateTraitT::is_kstate_trait);
    // *********** helper types ***********************************************************
    using KstateT = typename KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    // *********** helper values **********************************************************
    const auto site_basis_dim = SiteStateTraitT::site_basis_dim();
    // *********** [[expect]] *************************************************************
    assert(eigen_vector.n_rows == basis.size());
    // *********** prepare ****************************************************************
    std::vector<arma::cx_mat> density_matrix_all(n_threads);
    for (unsigned i = 0; i < n_threads; i++) {
        density_matrix_all[i] = arma::cx_mat(site_basis_dim, site_basis_dim, arma::fill::zeros);
    }
    // *********** filling ****************************************************************
    //const auto tp_fill_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
#pragma omp parallel for num_threads(n_threads) schedule(guided)
    for (arma::uword ket_kstate_idx = 0; ket_kstate_idx < eigen_vector.n_rows; ket_kstate_idx++) {
        const auto tid = omp_get_thread_num();
        calculate_reduced_density_operator_1_impl(basis, eigen_vector, ket_kstate_idx, density_matrix_all[tid]);
    }
    // const auto tp_fill_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //std::cout << "fill took     : " << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_fill_2 - tp_fill_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
    // *********** reduction **************************************************************
    //const auto tp_reduce_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
    for (unsigned d = 1; d < n_threads; d *=2) {
#pragma omp parallel num_threads(n_threads)
        {
            const auto tid = omp_get_thread_num();
            if (tid % (2 * d) == 0) {
                const auto idx1 = tid;
                const auto idx2 = tid + d;
                if (idx2 < n_threads) {
                    density_matrix_all[idx1] += density_matrix_all[idx2];
                }
            }
        }
    }
    //const auto tp_reduce_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //std::cout << "reduce took   : " << std::chrono::duration_cast<std::chrono::nanoseconds>(tp_reduce_2 - tp_reduce_1).count() / 1e6 << "ms" << std::endl; // performance debug sake
    // *********** return *****************************************************************
    return density_matrix_all[0];
}

} // end of namespace bfpt_common;
