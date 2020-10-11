#ifndef KSTATE_TRAIT_SITE_STATE_HPP
#define KSTATE_TRAIT_SITE_STATE_HPP

#include <type_traits>

namespace kstate_trait {

template<typename _SiteStateT>
struct TraitSiteState {
    static constexpr bool is_site_state_trait = false;
    using SiteStateT = _SiteStateT;
};

template<typename T>
struct IsTraitSiteState : std::false_type {
};

template<typename T>
struct IsTraitSiteState<TraitSiteState<T>> : std::true_type {
};

}

#endif
