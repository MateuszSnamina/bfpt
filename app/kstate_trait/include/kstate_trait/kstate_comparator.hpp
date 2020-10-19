#pragma once

#include <kstate_trait/trait_kstate.hpp>

#include <boost/algorithm/string/predicate.hpp>

#include <memory>
#include <type_traits>

// #######################################################################
// ## RangeComparator                                                   ##
// #######################################################################

namespace kstate_trait {

template<typename _KstateTraitT>
struct ViewComparator {
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;

    template <typename View1T, typename View2T>
    bool operator()(const View1T& lhs,
                    const View2T& rhs) const {
        return KstateTraitT::view_compare_less(lhs, rhs);
    }
};

}  // namespace kstate

// #######################################################################
// ## KstateComparator                                                  ##
// #######################################################################

namespace kstate_trait {

template<typename _KstateTraitT>
struct KstateComparator {
    static_assert(kstate_trait::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;

    using is_transparent = std::true_type;

    bool operator()(const std::shared_ptr<KstateT>& lhs, const std::shared_ptr<KstateT>& rhs) const {
            return kstate_trait::ViewComparator<KstateTraitT>()(KstateTraitT::to_range(*lhs), KstateTraitT::to_range(*rhs));
    }

    template<typename ViewT>
    bool operator()(const std::shared_ptr<KstateT>& lhs, const ViewT& rhs) const {
            return kstate_trait::ViewComparator<KstateTraitT>()(KstateTraitT::to_range(*lhs), rhs);
    }

    template<typename ViewT>
    bool operator()(const ViewT& lhs, const std::shared_ptr<KstateT>& rhs) const {
            return kstate_trait::ViewComparator<KstateTraitT>()(lhs, KstateTraitT::to_range(*rhs));
    }

};

}  // namespace kstate
