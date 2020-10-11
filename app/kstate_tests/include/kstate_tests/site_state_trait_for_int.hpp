#pragma once

#include <kstate_trait/trait_site_state.hpp>

namespace kstate_trait {

template<>
struct TraitSiteState<int> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = int;
};

} // end of namespace kstate
