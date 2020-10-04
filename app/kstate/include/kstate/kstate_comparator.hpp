#ifndef KSTATE_KSTATE_COMPARATOR_HPP
#define KSTATE_KSTATE_COMPARATOR_HPP

#include <kstate/trait_kstate.hpp>

#include <boost/algorithm/string/predicate.hpp>

#include <memory>
#include <type_traits>

// #######################################################################
// ## RangeComparator                                                   ##
// #######################################################################

namespace kstate {

struct RangeComparator {
    template <typename ConstRangeType1, typename ConstRangeType2>
    bool operator()(const ConstRangeType1& lhs,
                    const ConstRangeType2& rhs) const {
        return boost::lexicographical_compare(lhs, rhs);
    }
};

}  // namespace kstate

// #######################################################################
// ## KstateComparator                                                  ##
// #######################################################################

namespace kstate {

template<typename _KstateTraitT>
struct KstateComparator {
    static_assert(kstate::IsTraitKstate<_KstateTraitT>::value);
    static_assert(_KstateTraitT::is_kstate_trait);
    using KstateTraitT = _KstateTraitT;
    using KstateT = typename _KstateTraitT::KstateT;

    using is_transparent = std::true_type;

    bool operator()(std::shared_ptr<KstateT> lhs, std::shared_ptr<KstateT> rhs) const {
            return kstate::RangeComparator()(lhs->to_range(), rhs->to_range());
    }

    template<typename ForwardRangeT>
    bool operator()(std::shared_ptr<KstateT> lhs, ForwardRangeT rhs) const {
            return kstate::RangeComparator()(lhs->to_range(), rhs);
    }

    template<typename ForwardRangeT>
    bool operator()(ForwardRangeT lhs, std::shared_ptr<KstateT> rhs) const {
            return kstate::RangeComparator()(lhs, rhs->to_range());
    }

};

}  // namespace kstate

#endif // KSTATE_KSTATE_COMPARATOR_HPP
