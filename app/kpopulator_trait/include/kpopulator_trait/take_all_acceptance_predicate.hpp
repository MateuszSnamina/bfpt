#pragma once

#include<kstate_trait/trait_kstate.hpp>

// #######################################################################
// ## TakeAllAcceptancePredicate                                        ##
// #######################################################################

namespace kpopulator_trait {

template <typename _KstateTraitT>
class TakeAllAcceptancePredicate{
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);

   public:
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;

   public:
    bool operator()(const KstateT&) const {
        return true;
    };
};

}  // namespace kpopulator_trait
