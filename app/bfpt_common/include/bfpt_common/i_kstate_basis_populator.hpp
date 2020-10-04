#ifndef BFPT_COMMON_I_KSTATE_BASIS_POPULATOR_HPP
#define BFPT_COMMON_I_KSTATE_BASIS_POPULATOR_HPP

#include <kstate/basis.hpp>
#include <kstate/kstate_stl.hpp>

// #######################################################################
// ##  IKstateBasisPopulator                                            ##
// #######################################################################

namespace bfpt_common {

template<typename _KstateTraitT>
class IKstateBasisPopulator {
    static_assert(kstate::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
   public:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    using SiteStateT = typename KstateT::SiteStateT;
    using BasisT = kstate::Basis<KstateTraitT>;
   public:
    virtual kstate::KstateSet<KstateT> get_coupled_states(
        const KstateT& generator) const = 0;
};

}  // namespace bfpt_common

#endif
