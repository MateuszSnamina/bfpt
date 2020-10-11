#pragma once

#include <koperator_trait/trait_koperator.hpp>

#include <kstate_trait/trait_kstate.hpp>

#include <kbasis/basis.hpp>

#include <kstate_impl/range_op_unique_shift.hpp>

#include <chainkernel/operator_kernel.hpp>

#include <extensions/adaptors.hpp>

#include <armadillo>

#include <type_traits>
#include <cassert>
#include <complex>

//#include<chrono> // performance debug sake

// #######################################################################
// ## KernelDrivenKstateOperatorMatrix                                  ##
// #######################################################################

namespace koperator_impl {

template<typename _KstateTraitT>
class KernelDrivenKstateOperatorMatrix {
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
public:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    using SiteStateT = typename KstateT::SiteStateT;
    using BasisT = kbasis::Basis<KstateTraitT>;
public:
    KernelDrivenKstateOperatorMatrix(
            const size_t n_sites,
            chainkernel::OperatorKernel1<SiteStateTraitT> operator_kernel_1,
            chainkernel::OperatorKernel12<SiteStateTraitT> operator_kernel_12);
    void fill_kn_operator_builder_matrix_coll(
            const BasisT& basis,
            size_t n_col,
            arma::sp_cx_mat& kn_operator_builder_matrix,
            const unsigned k_n) const;
private:
    const size_t _n_sites;
    const chainkernel::OperatorKernel1<SiteStateTraitT> _operator_kernel_1;
    const chainkernel::OperatorKernel12<SiteStateTraitT> _operator_kernel_12;
};

}  // namespace koperator_impl


// #######################################################################
// ## KernelDrivenKstateOperatorMatrix -- impl                          ##
// #######################################################################

namespace koperator_impl {

template<typename _KstateTraitT>
KernelDrivenKstateOperatorMatrix<_KstateTraitT>::KernelDrivenKstateOperatorMatrix(
        const size_t n_sites,
        chainkernel::OperatorKernel1<SiteStateTraitT> operator_kernel_1,
        chainkernel::OperatorKernel12<SiteStateTraitT> operator_kernel_12)
    : _n_sites(n_sites),
      _operator_kernel_1(operator_kernel_1),
      _operator_kernel_12(operator_kernel_12) {
}

/*
 * In this implementation we assume thar
 * in the basis there are only elements that are unique shifted!
 */
template<typename _KstateTraitT>
void
KernelDrivenKstateOperatorMatrix<_KstateTraitT>::fill_kn_operator_builder_matrix_coll(
        const BasisT& basis,
        const size_t ket_kstate_idx,
        arma::sp_cx_mat& kn_operator_builder_matrix,
        const unsigned k_n) const {
    using namespace std::complex_literals;
    assert(kn_operator_builder_matrix.n_cols == kn_operator_builder_matrix.n_rows);
    assert(kn_operator_builder_matrix.n_rows == basis.size());
    const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
    assert(ket_kstate_ptr);
    const auto& ket_kstate_range = KstateTraitT::to_range(*ket_kstate_ptr);
    // ********** OFF-DIAG, KERNEL1 *********************************************
    for (size_t n_delta = 0; n_delta < _n_sites; n_delta++) {
        const auto ket_kernel_site_1 = *std::next(std::begin(ket_kstate_range), n_delta);
        const chainkernel::StateKernel1<SiteStateTraitT> ket_kernel{ket_kernel_site_1};
        const auto equal_range = _operator_kernel_1._half_off_diag_info.equal_range(ket_kernel);
        for (auto off_diag_node_it = equal_range.first; off_diag_node_it != equal_range.second; ++off_diag_node_it) {
            const auto& ket_kernel_re = off_diag_node_it->first;
            [[maybe_unused]] const auto& ket_kernel_site_1_re = ket_kernel_re.state_1;
            assert(ket_kernel_site_1_re == ket_kernel_site_1);
            const auto& couple_info = off_diag_node_it->second;
            const auto& kernel_coupling_coef = couple_info.coef;
            const auto& bra_kernel = couple_info.kernel_state;
            const auto& bra_kernel_site_1 = bra_kernel.state_1;
            const auto refined_holder_1 = extension::boost::adaptors::refined(n_delta, bra_kernel_site_1); // Must outlive bra_kstate_range.
            const auto bra_kstate_range = ket_kstate_range | refined_holder_1;
            const size_t bra_kstate_n_unique_shift = kstate::n_unique_shift(bra_kstate_range);
            const auto bra_kstate_range_unique_shifted = bra_kstate_range | extension::boost::adaptors::rotated(bra_kstate_n_unique_shift); // equivalent to `kstate::make_unique_shift(bra_kstate)`
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(bra_kstate_range_unique_shifted)) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm_1 = _n_sites * KstateTraitT::norm_factor(*basis.vec_index()[bra_kstate_idx]) * KstateTraitT::norm_factor(*basis.vec_index()[ket_kstate_idx]);
                const size_t bra_n_least_replication_shift = KstateTraitT::n_least_replication_shift(*basis.vec_index()[bra_kstate_idx]);
                const size_t bra_n_replicas = _n_sites / bra_n_least_replication_shift;
                const int exponent_n = (int)bra_kstate_n_unique_shift;
                const double exponent_r = 2 * arma::datum::pi * k_n / _n_sites * exponent_n;
                assert(k_n * bra_n_least_replication_shift % _n_sites == 0);
                const std::complex<double> neo_sum_phase_factors = std::exp(1.0i * exponent_r) * (double)bra_n_replicas;
                const std::complex<double> pre_norm_2 = pre_norm_1 * neo_sum_phase_factors;
                kn_operator_builder_matrix(bra_kstate_idx, ket_kstate_idx) += pre_norm_2 * kernel_coupling_coef;
            }
        } // end of `_half_off_diag_info` equal_range loop
    } // end of `Delta` loop
    // ********** OFF-DIAG, KERNEL12 ********************************************
    //std::chrono::high_resolution_clock::time_point tp_u_1, tp_u_2; // performance debug sake
    //std::chrono::high_resolution_clock::time_point tp_nu_1, tp_nu_2; // performance debug sake
    //double unique_shift_time = 0.0, not_unique_shift_time = 0.0; // performance debug sake
    //const auto tp_offdiag_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //tp_nu_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const auto ket_kernel_site_1 = *std::next(std::begin(ket_kstate_range), n_delta);
        const auto ket_kernel_site_2 = *std::next(std::begin(ket_kstate_range), n_delta_p1);
        const chainkernel::StateKernel12<SiteStateTraitT> ket_kernel{ket_kernel_site_1, ket_kernel_site_2};
        const auto equal_range = _operator_kernel_12._half_off_diag_info.equal_range(ket_kernel);
        for (auto off_diag_node_it = equal_range.first; off_diag_node_it != equal_range.second; ++off_diag_node_it) {
            const auto& ket_kernel_re = off_diag_node_it->first;
            [[maybe_unused]] const auto& ket_kernel_site_1_re = ket_kernel_re.state_1;
            [[maybe_unused]] const auto& ket_kernel_site_2_re = ket_kernel_re.state_2;
            assert(ket_kernel_site_1_re == ket_kernel_site_1);
            assert(ket_kernel_site_2_re == ket_kernel_site_2);
            const auto& couple_info = off_diag_node_it->second;
            const auto& kernel_coupling_coef = couple_info.coef;
            const auto& bra_kernel = couple_info.kernel_state;
            const auto& bra_kernel_site_1 = bra_kernel.state_1;
            const auto& bra_kernel_site_2 = bra_kernel.state_2;
            const auto refined_holder_1 = extension::boost::adaptors::refined(n_delta, bra_kernel_site_1); // Must outlive bra_kstate_range.
            const auto refined_holder_2 = extension::boost::adaptors::refined(n_delta_p1, bra_kernel_site_2); // Must outlive bra_kstate_range.
            const auto bra_kstate_range = ket_kstate_range | refined_holder_1 | refined_holder_2;
            //tp_nu_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
            //not_unique_shift_time += std::chrono::duration_cast<std::chrono::nanoseconds>(tp_nu_2 - tp_nu_1).count(); // performance debug sake
            //tp_u_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
            const size_t bra_kstate_n_unique_shift = kstate::n_unique_shift(bra_kstate_range);
            const auto bra_kstate_range_unique_shifted = bra_kstate_range | extension::boost::adaptors::rotated(bra_kstate_n_unique_shift); // equivalent to `kstate::make_unique_shift(bra_kstate)`
            //tp_u_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
            //unique_shift_time += std::chrono::duration_cast<std::chrono::nanoseconds>(tp_u_2 - tp_u_1).count(); // performance debug sake
            //tp_nu_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(bra_kstate_range_unique_shifted)) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm_1 = _n_sites * KstateTraitT::norm_factor(*basis.vec_index()[bra_kstate_idx]) * KstateTraitT::norm_factor(*basis.vec_index()[ket_kstate_idx]);
                const size_t bra_n_least_replication_shift = KstateTraitT::n_least_replication_shift(*basis.vec_index()[bra_kstate_idx]);
                const size_t bra_n_replicas = _n_sites / bra_n_least_replication_shift;
                const int exponent_n = (int)bra_kstate_n_unique_shift;
                const double exponent_r = 2 * arma::datum::pi * k_n / _n_sites * exponent_n;
                assert(k_n * bra_n_least_replication_shift % _n_sites == 0);
                const std::complex<double> neo_sum_phase_factors = std::exp(1.0i * exponent_r) * (double)bra_n_replicas;
                const std::complex<double> pre_norm_2 = pre_norm_1 * neo_sum_phase_factors;
                kn_operator_builder_matrix(bra_kstate_idx, ket_kstate_idx) += pre_norm_2 * kernel_coupling_coef;
            }
        } // end of `_half_off_diag_info` equal_range loop
    } // end of `Delta` loop
    //const auto tp_offdiag_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //tp_nu_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
    //not_unique_shift_time += std::chrono::duration_cast<std::chrono::nanoseconds>(tp_nu_2 - tp_nu_1).count(); // performance debug sake
    //const double offdiag_time = std::chrono::duration_cast<std::chrono::nanoseconds>(tp_offdiag_2 - tp_offdiag_1).count(); // performance debug sake
    //std::cout << "[OFF-DIAG] [TIMING]: offdiag_time          : " << offdiag_time << std::endl; // performance debug sake
    //std::cout << "[OFF-DIAG] [TIMING]: unique_shift_time     : " << unique_shift_time << std::endl; // performance debug sake
    //std::cout << "[OFF-DIAG] [TIMING]: not_unique_shift_time : " << not_unique_shift_time << std::endl; // performance debug sake
    // ********** ON-DIAG, KERNEL1 **********************************************
    for (size_t n_delta = 0; n_delta < _n_sites; n_delta++) {
        const auto ket_kernel_site_1 = *std::next(std::begin(ket_kstate_range), n_delta);
        const chainkernel::StateKernel1<SiteStateTraitT> ket_kernel{ket_kernel_site_1};
        if (_operator_kernel_1._diag_info.count(ket_kernel)) {
            const auto kernel_diag_coef = _operator_kernel_1._diag_info.at(ket_kernel);
            kn_operator_builder_matrix(ket_kstate_idx, ket_kstate_idx) += kernel_diag_coef / 2; // factor '/2' is as we build matrix M such as H = M + M^T.
        }
    }  // end of `Delta` loop
    // ********** ON-DIAG, KERNEL12 *********************************************
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const auto ket_kernel_site_1 = *std::next(std::begin(ket_kstate_range), n_delta);
        const auto ket_kernel_site_2 = *std::next(std::begin(ket_kstate_range), n_delta_p1);
        const chainkernel::StateKernel12<SiteStateTraitT> ket_kernel{ket_kernel_site_1, ket_kernel_site_2};
        if (_operator_kernel_12._diag_info.count(ket_kernel)) {
            const auto kernel_diag_coef = _operator_kernel_12._diag_info.at(ket_kernel);
            kn_operator_builder_matrix(ket_kstate_idx, ket_kstate_idx) += kernel_diag_coef / 2; // factor '/2' is as we build matrix M such as H = M + M^T.
        }
    }  // end of `Delta` loop
}

}  // namespace koperator_impl

// #######################################################################
// ## TraitsFor koperator_impl::KernelDrivenKstateOperatorMatrix        ##
// #######################################################################

namespace koperator_trait {

template<typename _KstateTraitT>
struct TraitKoperator<koperator_impl::KernelDrivenKstateOperatorMatrix<_KstateTraitT>> {
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
    // the is_kstate_trait flag:
    static constexpr bool is_koperator_trait = true;
    // helper types:
    using KstateTraitT = _KstateTraitT;
    using KoperatorT = koperator_impl::KernelDrivenKstateOperatorMatrix<_KstateTraitT>;
    using BasisT = kbasis::Basis<KstateTraitT>;
    // function being the public API:
    static void fill_kn_operator_builder_matrix_coll(
            const KoperatorT& koperator,
            const BasisT& basis,
            size_t n_col,
            arma::sp_cx_mat& kn_operator_builder_matrix,
            const unsigned k_n) {
        koperator.fill_kn_operator_builder_matrix_coll(basis,
                                                       n_col,
                                                       kn_operator_builder_matrix,
                                                       k_n);
    }

};

}
