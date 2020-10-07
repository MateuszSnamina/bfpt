#pragma once

#include <kstate/trait_site_state.hpp>

namespace kstate {

template<>
struct TraitSiteState<int> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = int;
};

} // end of namespace kstate
