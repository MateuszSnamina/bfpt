#ifndef KSTATE_TRAIT_SITE_STATE_HPP
#define KSTATE_TRAIT_SITE_STATE_HPP

namespace kstate {

template<typename _SiteStateT>
struct TraitSiteState {
    static constexpr bool is_site_state_trait = false;
    using SiteStateT = _SiteStateT;
};

}

#endif
