#ifndef BFPT_COMMON_GENERIC_KSTATE_HAMILTONIAN_HPP
#define BFPT_COMMON_GENERIC_KSTATE_HAMILTONIAN_HPP

#include <bfpt_common/operator_kernel.hpp>
#include <bfpt_common/i_kstate_operator.hpp>
#include <bfpt_common/i_kstate_populator.hpp>
#include <bfpt_common/generate_pt_basis.hpp> //TODO check if needed!

#include <kstate/unique_shift.hpp>
#include <kstate/kstate_abstract.hpp>
#include <kstate/remove_cvref.hpp>
#include <kstate/is_base_of_template.hpp>
#include <kstate/kstate_abstract.hpp>
#include <extensions/adaptors.hpp>

#include <armadillo> //TODO check if needed!

#include <type_traits>
#include <cassert>
#include <complex> //TODO check if needed!

//#include<chrono> // performance debug sake

// #######################################################################
// ## GenericKstateHamiltonian                                          ##
// #######################################################################

namespace bfpt_common {

template<typename _KstateT>
class GenericKstatePopulator :
        public bfpt_common::IKstatePopulator<_KstateT> {
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
    GenericKstatePopulator(
            const size_t n_sites,
            OperatorKernel1<SiteStateT> hamiltonian_kernel_1,
            OperatorKernel12<SiteStateT> hamiltonian_kernel_12);
    kstate::KstateSet<KstateT> get_coupled_states(
            const KstateT& generator) const override;
private:
    const size_t _n_sites;
    const OperatorKernel1<SiteStateT> _hamiltonian_kernel_1;
    const OperatorKernel12<SiteStateT> _hamiltonian_kernel_12;
};

}  // namespace bfpt_common

// #######################################################################
// ## GenericKstateHamiltonian                                          ##
// #######################################################################

using namespace std::complex_literals;//TODO remove!

namespace bfpt_common {

template<typename _SiteStateT>
GenericKstatePopulator<_SiteStateT>::GenericKstatePopulator(
        const size_t n_sites,
        OperatorKernel1<SiteStateT> hamiltonian_kernel_1,
        OperatorKernel12<SiteStateT> hamiltonian_kernel_12)
    : _n_sites(n_sites),
      _hamiltonian_kernel_1(hamiltonian_kernel_1),
      _hamiltonian_kernel_12(hamiltonian_kernel_12) {
}

template<typename _SiteStateT>
kstate::KstateSet<typename GenericKstatePopulator<_SiteStateT>::KstateT>
GenericKstatePopulator<_SiteStateT>::get_coupled_states(
        const KstateT& generator) const {
    kstate::KstateSet<KstateT> result;
    assert(generator.n_sites() == _n_sites);
    const auto generator_range = generator.to_range();
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const auto ket_kernel_site_1 = *std::next(std::begin(generator_range), n_delta);
        const auto ket_kernel_site_2 = *std::next(std::begin(generator_range), n_delta_p1);
        const StateKernel12<SiteStateT> ket_kernel{ket_kernel_site_1, ket_kernel_site_2};
        const auto equal_range = _hamiltonian_kernel_12._full_off_diag_info.equal_range(ket_kernel);
        for (auto off_diag_node_it = equal_range.first; off_diag_node_it != equal_range.second; ++off_diag_node_it) {
            const auto& ket_12_re = off_diag_node_it->first;
            [[maybe_unused]] const auto& ket_kernel_site_1_re = ket_12_re.state_1;
            [[maybe_unused]] const auto& ket_kernel_site_2_re = ket_12_re.state_2;
            assert(ket_kernel_site_1_re == ket_kernel_site_1);
            assert(ket_kernel_site_2_re == ket_kernel_site_2);
            const auto& couple_info = off_diag_node_it->second;
            //[[maybe_unused]] const auto& kernel_coupling_coef = couple_info.coef;
            const auto& bra_kernel = couple_info.kernel_state;
            const auto& bra_kernel_site_1 = bra_kernel.state_1;
            const auto& bra_kernel_site_2 = bra_kernel.state_2;
            //TODO: if (kernel_coupling_coef !=0 ){ FILL }
            const auto conjugated_range =
                    generator_range |
                    extension::boost::adaptors::refined(n_delta, bra_kernel_site_1) |
                    extension::boost::adaptors::refined(n_delta_p1, bra_kernel_site_2);
            const auto conjugated_range_unique_shifted = kstate::make_unique_shift(conjugated_range);
            const auto conjugated_kstate_ptr = std::make_shared<KstateT>(conjugated_range_unique_shifted, kstate::ctr_from_range);
            result.insert(conjugated_kstate_ptr);
        } // end of `_full_off_diag_info` equal_range loop
    }  // end of `Delta` loop
    return result;
}

}  // namespace bfpt_common

#endif
