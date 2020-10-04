#ifndef KSTATE_TRAIT_KSTATE_HPP
#define KSTATE_TRAIT_KSTATE_HPP

namespace kstate {

template<typename _KstateT>
struct TraitKstate {
    static constexpr bool is_kstate_trait = false;
    using KstateT = _KstateT;
};

}

#endif
