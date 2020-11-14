#pragma once

#include <kpopulator_trait/trait_kpopulator.hpp>

#include <kstate_view_amend_spec/amend_spec.hpp>

#include <kstate_trait/trait_kstate.hpp>
#include <kstate_trait/kstate_stl.hpp>

#include <kbasis/basis.hpp>

#include <chainkernel/operator_kernel.hpp>

#include <type_traits>
#include <functional>
#include <cassert>

//#include<chrono> // performance debug sake

// #######################################################################
// ## KernelDrivenKstateBasisPopulator                                  ##
// #######################################################################

namespace kpopulator_impl {

template <typename _KstateTraitT>
class KernelDrivenKstateBasisPopulator {
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);

   public:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    using SiteStateT = typename KstateT::SiteStateT;
    using BasisT = kbasis::Basis<KstateTraitT>;

   public:
    KernelDrivenKstateBasisPopulator(
        const size_t n_sites,
        chainkernel::OperatorKernel1<SiteStateTraitT> operator_kernel_1,
        chainkernel::OperatorKernel12<SiteStateTraitT> operator_kernel_12);
    KernelDrivenKstateBasisPopulator(
        const size_t n_sites,
        chainkernel::OperatorKernel1<SiteStateTraitT> operator_kernel_1,
        chainkernel::OperatorKernel12<SiteStateTraitT> operator_kernel_12,
        std::function<bool(KstateT)> acceptance_predicate);

    kstate_trait::KstateSet<KstateTraitT> get_coupled_states(
        const KstateT& generator,
        const unsigned n_k) const;

   private:
    const size_t _n_sites;
    const chainkernel::OperatorKernel1<SiteStateTraitT> _operator_kernel_1;
    const chainkernel::OperatorKernel12<SiteStateTraitT> _operator_kernel_12;
    const std::function<bool(KstateT)> _acceptance_predicate;
};

}  // namespace kpopulator_impl

// #######################################################################
// ## KernelDrivenKstateBasisPopulator -- impl                          ##
// #######################################################################

namespace kpopulator_impl {

template <typename _KstateTraitT>
KernelDrivenKstateBasisPopulator<_KstateTraitT>::KernelDrivenKstateBasisPopulator(
    const size_t n_sites,
    chainkernel::OperatorKernel1<SiteStateTraitT> operator_kernel_1,
    chainkernel::OperatorKernel12<SiteStateTraitT> operator_kernel_12,
    std::function<bool(KstateT)> acceptance_predicate)
    : _n_sites(n_sites),
      _operator_kernel_1(operator_kernel_1),
      _operator_kernel_12(operator_kernel_12),
      _acceptance_predicate(acceptance_predicate) {
}

template <typename _KstateTraitT>
KernelDrivenKstateBasisPopulator<_KstateTraitT>::KernelDrivenKstateBasisPopulator(
    const size_t n_sites,
    chainkernel::OperatorKernel1<SiteStateTraitT> operator_kernel_1,
    chainkernel::OperatorKernel12<SiteStateTraitT> operator_kernel_12)
    : _n_sites(n_sites),
      _operator_kernel_1(operator_kernel_1),
      _operator_kernel_12(operator_kernel_12),
      _acceptance_predicate([](const KstateT&) -> bool { return true; }) {
}

template <typename _KstateTraitT>
kstate_trait::KstateSet<typename KernelDrivenKstateBasisPopulator<_KstateTraitT>::KstateTraitT>
KernelDrivenKstateBasisPopulator<_KstateTraitT>::get_coupled_states(
    const KstateT& generator,
    const unsigned n_k) const {
    kstate_trait::KstateSet<KstateTraitT> result;
    assert(KstateTraitT::n_sites(generator) == _n_sites);
    // ********** OFF-DIAG, KERNEL12 ********************************************
    const auto generator_view = KstateTraitT::to_view(generator);
    for (size_t n_delta = 0, n_delta_p1 = 1; n_delta < _n_sites; n_delta++, n_delta_p1 = (n_delta + 1) % _n_sites) {
        const auto generator_site_1 = KstateTraitT::view_n_th_site_state(generator_view, n_delta);
        const auto generator_site_2 = KstateTraitT::view_n_th_site_state(generator_view, n_delta_p1);
        const chainkernel::StateKernel12<SiteStateTraitT> ket_kernel{generator_site_1, generator_site_2};
        const auto equal_range = _operator_kernel_12._full_off_diag_info.equal_range(ket_kernel);
        for (auto off_diag_node_it = equal_range.first; off_diag_node_it != equal_range.second; ++off_diag_node_it) {
            const auto& ket_12_re = off_diag_node_it->first;
            [[maybe_unused]] const auto& ket_kernel_site_1_re = ket_12_re.state_1;
            [[maybe_unused]] const auto& ket_kernel_site_2_re = ket_12_re.state_2;
            assert(ket_kernel_site_1_re == generator_site_1);
            assert(ket_kernel_site_2_re == generator_site_2);
            const auto& couple_info = off_diag_node_it->second;
            //[[maybe_unused]] const auto& kernel_coupling_coef = couple_info.coef;
            const auto& bra_kernel = couple_info.kernel_state;
            const auto& bra_kernel_site_1 = bra_kernel.state_1;
            const auto& bra_kernel_site_2 = bra_kernel.state_2;
            //TODO: if (kernel_coupling_coef !=0 ){ FILL }
            const auto refined_holder_1 = kstate_view_amend_spec::refined(n_delta, bra_kernel_site_1);             // Must outlive conjugated_view.
            const auto refined_holder_2 = kstate_view_amend_spec::refined(n_delta_p1, bra_kernel_site_2);          // Must outlive conjugated_view.
            const auto conjugated_view_preproduct = KstateTraitT::refined_view(generator_view, refined_holder_1);  // Must outlive conjugated_view.
            const auto conjugated_view = KstateTraitT::refined_view(conjugated_view_preproduct, refined_holder_2);
            if (KstateTraitT::is_prolific(conjugated_view, n_k)) {
                const size_t conjugated_view_n_unique_shift = KstateTraitT::view_n_unique_shift(conjugated_view);
                const auto rotation_spec = kstate_view_amend_spec::rotated(conjugated_view_n_unique_shift);              // Must outlive conjugated_view_unique_shifted.
                const auto conjugated_view_unique_shifted = KstateTraitT::rotated_view(conjugated_view, rotation_spec);  // equivalent to `kstate::make_unique_shift(conjugated_view)`
                const auto conjugated_kstate_ptr = KstateTraitT::shared_from_view(conjugated_view_unique_shifted);
                if (_acceptance_predicate(*conjugated_kstate_ptr)) {
                    result.insert(conjugated_kstate_ptr);
                }
            }
        }  // end of `_full_off_diag_info` equal_range loop
    }      // end of `Delta` loop
    // ********** OFF-DIAG, KERNEL1  ********************************************
    //TODO kernel 1!
    // **************************************************************************
    return result;
}

}  // namespace kpopulator_impl

// #######################################################################
// ## TraitsFor kpopulator_impl::KernelDrivenKstateBasisPopulator       ##
// #######################################################################

namespace kpopulator_trait {

template <typename _KstateTraitT>
struct TraitKpopulator<kpopulator_impl::KernelDrivenKstateBasisPopulator<_KstateTraitT>> {
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
    // the is_kstate_trait flag:
    static constexpr bool is_kpopulator_trait = true;
    // helper types:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename KstateTraitT::KstateT;
    using KpopulatorT = kpopulator_impl::KernelDrivenKstateBasisPopulator<_KstateTraitT>;
    using BasisT = kbasis::Basis<KstateTraitT>;
    // function being the public API:
    static kstate_trait::KstateSet<KstateTraitT> get_coupled_states(
        const KpopulatorT& kpopulator,
        const KstateT& generator,
        const unsigned n_k) {
        return kpopulator.get_coupled_states(generator, n_k);
    }
};

}  // namespace kpopulator_trait
