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
    virtual void push_back_coupled_states_to_basis(
        const KstateT& generator,
        BasisT& basis) const = 0;
};

}  // namespace bfpt_common

#endif
