#pragma once

#include <kstate/trait_site_state.hpp>

namespace kstate {

template<>
struct TraitSiteState<double> {
    static constexpr bool is_site_state_trait = true;
    using SiteStateT = double;
};

} // end of namespace kstate
