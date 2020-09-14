#ifndef BFPT_COMMON_I_DYNAMIC_UNIQUE_KSTATE_POPULATOR_HPP
#define BFPT_COMMON_I_DYNAMIC_UNIQUE_KSTATE_POPULATOR_HPP

#include <kstate/remove_cvref.hpp>
#include <kstate/is_base_of_template.hpp>
#include <kstate/kstate_abstract.hpp>
#include <kstate/basis.hpp>

// #######################################################################
// ##  IDynamicUniqueKstatePopulator                                    ##
// #######################################################################

namespace bfpt_common {

template<typename _KstateT>
class IKstatePopulator {
    static_assert(!std::is_const<_KstateT>::value);
    static_assert(!std::is_volatile<_KstateT>::value);
    static_assert(!std::is_reference<_KstateT>::value);
    static_assert(kstate::is_base_of_template_v<_KstateT, kstate::Kstate>);
   public:
    using KstateT = _KstateT;
    using SiteStateT = typename kstate::remove_cvref_t<KstateT>::SiteType;
    using BasisT = kstate::Basis<KstateT>;
   public:
    virtual void push_back_coupled_states_to_basis(
        const KstateT& generator,
        BasisT& basis) const = 0;
};

}  // namespace bfpt_common

#endif
