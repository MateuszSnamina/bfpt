#pragma once

#include <type_traits>

namespace kstate_trait {

template <typename _SiteStateT>
struct TraitSiteState {
    static constexpr bool is_site_state_trait = false;
    using SiteStateT = _SiteStateT;
};

template <typename T>
struct IsTraitSiteState : std::false_type {
};

template <typename T>
struct IsTraitSiteState<TraitSiteState<T>> : std::true_type {
};

}  // namespace kstate_trait
