#ifndef BFPT_COMMON_I_DYNAMIC_UNIQUE_KSTATE_POPULATOR_HPP
#define BFPT_COMMON_I_DYNAMIC_UNIQUE_KSTATE_POPULATOR_HPP

#include <kstate/basis.hpp>
#include <kstate/kstate_concrete.hpp>

// #######################################################################
// ##  IDynamicUniqueKstatePopulator                                    ##
// #######################################################################

namespace bfpt_common {

template <typename Element> // TODO: change to: template <typename _SiteType>
class IDynamicUniqueKstatePopulator {
   public:
    virtual void push_back_coupled_states_to_basis(
        const kstate::DynamicUniqueKstate<Element>& generator,
        kstate::Basis<kstate::DynamicUniqueKstate<Element>>& basis) const = 0;
};

}  // namespace bfpt_common

#endif
