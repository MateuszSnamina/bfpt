#ifndef BFPT_COMMON_KERNEL_DRIVEN_KSTATE_OPERATOR_HPP
#define BFPT_COMMON_KERNEL_DRIVEN_KSTATE_OPERATOR_HPP

#include <bfpt_common/operator_kernel.hpp>
#include <bfpt_common/i_kstate_operator_matrix.hpp>

#include <kstate/unique_shift.hpp>
#include <kstate/kstate_abstract.hpp>
#include <kstate/remove_cvref.hpp>
#include <kstate/is_base_of_template.hpp>
#include <kstate/kstate_abstract.hpp>
#include <extensions/adaptors.hpp>

#include <armadillo>

#include <type_traits>
#include <cassert>
#include <complex>

//#include<chrono> // performance debug sake

// #######################################################################
// ## KernelDrivenKstateOperatorMatrix                                  ##
// #######################################################################

namespace bfpt_common {

template<typename _KstateT>
class KernelDrivenKstateOperatorMatrix : public bfpt_common::IKstateOperatorMatrix<_KstateT> {
    static_assert(!std::is_array_v<_KstateT>);
    static_assert(!std::is_function_v<_KstateT>);
    static_assert(!std::is_void_v<std::decay<_KstateT>>);
    static_assert(!std::is_null_pointer_v<std::decay<_KstateT>>);
    static_assert(!std::is_enum_v<std::decay<_KstateT>>);
    static_assert(!std::is_union_v<std::decay<_KstateT>>);
    static_assert(std::is_class_v<std::decay<_KstateT>>);
    static_assert(!std::is_pointer_v<std::decay<_KstateT>>);
    static_assert(!std::is_member_object_pointer_v<_KstateT>);
    static_assert(!std::is_member_function_pointer_v<_KstateT>);
    static_assert(!std::is_const_v<_KstateT>);
    static_assert(!std::is_volatile_v<_KstateT>);
    static_assert(!std::is_reference_v<_KstateT>);
    static_assert(kstate::is_base_of_template_v<_KstateT, kstate::Kstate>);
public:
    using KstateT = _KstateT;
    using SiteStateT = typename kstate::remove_cvref_t<KstateT>::SiteType;
    using BasisT = kstate::Basis<KstateT>;
public:
    KernelDrivenKstateOperatorMatrix(
            const size_t n_sites,
            OperatorKernel1<SiteStateT> operator_kernel_1,
            OperatorKernel12<SiteStateT> operator_kernel_12);
    void fill_kn_operator_builder_matrix_coll(
            const BasisT& basis,
            size_t n_col,
            arma::sp_cx_mat& kn_operator_builder_matrix,
            const unsigned k_n) const override;
private:
    const size_t _n_sites;
    const OperatorKernel1<SiteStateT> _operator_kernel_1;
    const OperatorKernel12<SiteStateT> _operator_kernel_12;
};

}  // namespace bfpt_common


// #######################################################################
// ## KernelDrivenKstateOperatorMatrix -- impl                          ##
// #######################################################################

namespace bfpt_common {

template<typename _SiteStateT>
KernelDrivenKstateOperatorMatrix<_SiteStateT>::KernelDrivenKstateOperatorMatrix(
        const size_t n_sites,
        OperatorKernel1<SiteStateT> operator_kernel_1,
        OperatorKernel12<SiteStateT> operator_kernel_12)
    : _n_sites(n_sites),
      _operator_kernel_1(operator_kernel_1),
      _operator_kernel_12(operator_kernel_12) {
}

/*
 * In this implementation we assume thar
 * in the basis there are only elements that are unique shifted!
 */
template<typename _SiteStateT>
void
KernelDrivenKstateOperatorMatrix<_SiteStateT>::fill_kn_operator_builder_matrix_coll(
        const BasisT& basis,
        const size_t ket_kstate_idx,
        arma::sp_cx_mat& kn_operator_builder_matrix,
        const unsigned k_n) const {
    using namespace std::complex_literals;
    assert(kn_operator_builder_matrix.n_cols == kn_operator_builder_matrix.n_rows);
    assert(kn_operator_builder_matrix.n_rows == basis.size());
    const auto ket_kstate_ptr = basis.vec_index()[ket_kstate_idx];
    assert(ket_kstate_ptr);
    const auto& ket_kstate = (*ket_kstate_ptr).to_range();
    // ********** OFF-DIAG, KERNEL1 *********************************************
    for (size_t n_delta = 0; n_delta < _n_sites; n_delta++) {
        const auto ket_kernel_site_1 = *std::next(std::begin(ket_kstate), n_delta);
        const StateKernel1<SiteStateT> ket_kernel{ket_kernel_site_1};
        const auto equal_range = _operator_kernel_1._half_off_diag_info.equal_range(ket_kernel);
        for (auto off_diag_node_it = equal_range.first; off_diag_node_it != equal_range.second; ++off_diag_node_it) {
            const auto& ket_kernel_re = off_diag_node_it->first;
            [[maybe_unused]] const auto& ket_kernel_site_1_re = ket_kernel_re.state_1;
            assert(ket_kernel_site_1_re == ket_kernel_site_1);
            const auto& couple_info = off_diag_node_it->second;
            const auto& kernel_coupling_coef = couple_info.coef;
            const auto& bra_kernel = couple_info.kernel_state;
            const auto& bra_kernel_site_1 = bra_kernel.state_1;
            const auto bra_kstate = ket_kstate
                    | extension::boost::adaptors::refined(n_delta, bra_kernel_site_1);
            const size_t bra_n_unique_shift = kstate::n_unique_shift(bra_kstate);
            const auto bra_kstate_unique_shifted = bra_kstate | extension::boost::adaptors::rotated(bra_n_unique_shift); // equivalent to `kstate::make_unique_shift(bra_kstate)`
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(bra_kstate_unique_shifted)) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm_1 = _n_sites * basis.vec_index()[bra_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
                const size_t bra_n_least_replication_shift = basis.vec_index()[bra_kstate_idx]->n_least_replication_shift();
                const size_t bra_n_replicas = _n_sites / bra_n_least_replication_shift;
                const int exponent_n = (int)bra_n_unique_shift;
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
        const auto ket_kernel_site_1 = *std::next(std::begin(ket_kstate), n_delta);
        const auto ket_kernel_site_2 = *std::next(std::begin(ket_kstate), n_delta_p1);
        const StateKernel12<SiteStateT> ket_kernel{ket_kernel_site_1, ket_kernel_site_2};
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
            const auto bra_kstate = ket_kstate
                    | extension::boost::adaptors::refined(n_delta, bra_kernel_site_1)
                    | extension::boost::adaptors::refined(n_delta_p1, bra_kernel_site_2);
            //tp_nu_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
            //not_unique_shift_time += std::chrono::duration_cast<std::chrono::nanoseconds>(tp_nu_2 - tp_nu_1).count(); // performance debug sake
            //tp_u_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
            const size_t bra_n_unique_shift = kstate::n_unique_shift(bra_kstate);
            const auto bra_kstate_unique_shifted = bra_kstate | extension::boost::adaptors::rotated(bra_n_unique_shift); // equivalent to `kstate::make_unique_shift(bra_kstate)`
            //tp_u_2 = std::chrono::high_resolution_clock::now(); // performance debug sake
            //unique_shift_time += std::chrono::duration_cast<std::chrono::nanoseconds>(tp_u_2 - tp_u_1).count(); // performance debug sake
            //tp_nu_1 = std::chrono::high_resolution_clock::now(); // performance debug sake
            if (const auto& bra_kstate_optional_idx = basis.find_element_and_get_its_ra_index(bra_kstate_unique_shifted)) {
                const auto bra_kstate_idx = *bra_kstate_optional_idx;
                double pre_norm_1 = _n_sites * basis.vec_index()[bra_kstate_idx]->norm_factor() * basis.vec_index()[ket_kstate_idx]->norm_factor();
                const size_t bra_n_least_replication_shift = basis.vec_index()[bra_kstate_idx]->n_least_replication_shift();
                const size_t bra_n_replicas = _n_sites / bra_n_least_replication_shift;
                const int exponent_n = (int)bra_n_unique_shift;
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
        const auto ket_kernel_site_1 = *std::next(std::begin(ket_kstate), n_delta);
        const StateKernel1<SiteStateT> ket_kernel{ket_kernel_site_1};
        if (_operator_kernel_1._diag_info.count(ket_kernel)) {
            const auto kernel_diag_coef = _operator_kernel_1._diag_info.at(ket_kernel);
            const double pre_norm_1 = _n_sites * ket_kstate_ptr->norm_factor() * ket_kstate_ptr->norm_factor();
            const double pre_norm_2 = pre_norm_1 * (_n_sites / ket_kstate_ptr->n_least_replication_shift());
            //TODO is not pre_norm_2 ALWAYS equal to 1.0?
            kn_operator_builder_matrix(ket_kstate_idx, ket_kstate_idx) += pre_norm_2 * kernel_diag_coef / 2; // factor '/2' is as we build matrix M such as H = M + M^T.
        }
    }  // end of `Delta` loop
    // ********** ON-DIAG, KERNEL12 *********************************************
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const auto ket_kernel_site_1 = *std::next(std::begin(ket_kstate), n_delta);
        const auto ket_kernel_site_2 = *std::next(std::begin(ket_kstate), n_delta_p1);
        const StateKernel12<SiteStateT> ket_kernel{ket_kernel_site_1, ket_kernel_site_2};
        if (_operator_kernel_12._diag_info.count(ket_kernel)) {
            const auto kernel_diag_coef = _operator_kernel_12._diag_info.at(ket_kernel);
            const double pre_norm_1 = _n_sites * ket_kstate_ptr->norm_factor() * ket_kstate_ptr->norm_factor();
            const double pre_norm_2 = pre_norm_1 * (_n_sites / ket_kstate_ptr->n_least_replication_shift());
            //TODO is not pre_norm_2 ALWAYS equal to 1.0?
            kn_operator_builder_matrix(ket_kstate_idx, ket_kstate_idx) += pre_norm_2 * kernel_diag_coef / 2; // factor '/2' is as we build matrix M such as H = M + M^T.
        }
    }  // end of `Delta` loop
}

}  // namespace bfpt_common

#endif
