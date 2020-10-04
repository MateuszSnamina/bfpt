#ifndef BFPT_COMMON_KERNEL_DRIVEN_KSTATE_BASIS_POPULATOR_HPP
#define BFPT_COMMON_KERNEL_DRIVEN_KSTATE_BASIS_POPULATOR_HPP

#include <bfpt_common/operator_kernel.hpp>
#include <bfpt_common/i_kstate_basis_populator.hpp>

#include <kstate/unique_shift.hpp>
#include <kstate/kstate_abstract.hpp>
#include <kstate/kstate_concrete.hpp>
#include <extensions/adaptors.hpp>

#include <type_traits>
#include <cassert>

//#include<chrono> // performance debug sake

// #######################################################################
// ## KernelDrivenKstateBasisPopulator                                  ##
// #######################################################################

namespace bfpt_common {

template<typename _KstateTraitT>
class KernelDrivenKstateBasisPopulator : public bfpt_common::IKstateBasisPopulator<_KstateTraitT> {
    static_assert(kstate::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
public:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    using SiteStateT = typename KstateT::SiteStateT;
    using BasisT = kstate::Basis<KstateTraitT>;
public:
    KernelDrivenKstateBasisPopulator(
            const size_t n_sites,
            OperatorKernel1<SiteStateTraitT> operator_kernel_1,
            OperatorKernel12<SiteStateTraitT> operator_kernel_12);
    kstate::KstateSet<KstateT> get_coupled_states(
            const KstateT& generator) const override;
private:
    const size_t _n_sites;
    const OperatorKernel1<SiteStateTraitT> _operator_kernel_1;
    const OperatorKernel12<SiteStateTraitT> _operator_kernel_12;
};

}  // namespace bfpt_common

// #######################################################################
// ## KernelDrivenKstateBasisPopulator -- impl                          ##
// #######################################################################

namespace bfpt_common {

template<typename _KstateTraitT>
KernelDrivenKstateBasisPopulator<_KstateTraitT>::KernelDrivenKstateBasisPopulator(
        const size_t n_sites,
        OperatorKernel1<SiteStateTraitT> operator_kernel_1,
        OperatorKernel12<SiteStateTraitT> operator_kernel_12)
    : _n_sites(n_sites),
      _operator_kernel_1(operator_kernel_1),
      _operator_kernel_12(operator_kernel_12) {
}

template<typename _KstateTraitT>
kstate::KstateSet<typename KernelDrivenKstateBasisPopulator<_KstateTraitT>::KstateT>
KernelDrivenKstateBasisPopulator<_KstateTraitT>::get_coupled_states(
        const KstateT& generator) const {
    kstate::KstateSet<KstateT> result;
    assert(generator.n_sites() == _n_sites);
    // ********** OFF-DIAG, KERNEL12 ********************************************
    const auto generator_range = generator.to_range();
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const auto ket_kernel_site_1 = *std::next(std::begin(generator_range), n_delta);
        const auto ket_kernel_site_2 = *std::next(std::begin(generator_range), n_delta_p1);
        const StateKernel12<SiteStateTraitT> ket_kernel{ket_kernel_site_1, ket_kernel_site_2};
        const auto equal_range = _operator_kernel_12._full_off_diag_info.equal_range(ket_kernel);
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
    // ********** OFF-DIAG, KERNEL1  ********************************************
    //TODO kernel 1!
    // **************************************************************************
    return result;
}

}  // namespace bfpt_common

#endif
