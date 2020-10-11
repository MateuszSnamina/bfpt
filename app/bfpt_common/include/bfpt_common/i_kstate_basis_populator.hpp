#pragma once

#include <kbasis/basis.hpp>

#include <kstate_trait/kstate_stl.hpp>

// #######################################################################
// ##  IKstateBasisPopulator                                            ##
// #######################################################################

namespace bfpt_common {

template<typename _KstateTraitT>
class IKstateBasisPopulator {
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
   public:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;
    using SiteStateTraitT = typename KstateT::SiteStateTraitT;
    using SiteStateT = typename KstateT::SiteStateT;
    using BasisT = kbasis::Basis<KstateTraitT>;
   public:
    virtual kstate_trait::KstateSet<KstateTraitT> get_coupled_states(
        const KstateT& generator) const = 0;
};

}  // namespace bfpt_common
