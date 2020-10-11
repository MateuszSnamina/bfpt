#ifndef KSTATE_TRAIT_KSTATE_HPP
#define KSTATE_TRAIT_KSTATE_HPP

#include <type_traits>

namespace kstate_trait {

template<typename _KstateT>
struct TraitKstate {
    static constexpr bool is_kstate_trait = false;
    using KstateT = _KstateT;
};

template<typename T>
struct IsTraitKstate : std::false_type{
};

template<typename T>
struct IsTraitKstate<TraitKstate<T>> : std::true_type{
};

} // end of namespace kstate

#endif
