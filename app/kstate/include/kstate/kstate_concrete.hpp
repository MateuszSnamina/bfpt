#ifndef KSTATE_KSTATE_CONCRETE_HPP
#define KSTATE_KSTATE_CONCRETE_HPP

#include <kstate/trait_site_state.hpp>
#include <kstate/trait_kstate.hpp>
#include <kstate/kstate_abstract.hpp>

#include <extensions/range_streamer.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm/search.hpp>
#include <boost/range/any_range.hpp>

#include <cassert>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// #######################################################################
// ## DynamicKstate                                                     ##
// #######################################################################

/*
 * DynamicKstate<SiteStateT> is a concrete subclass of
 * Kstate<SiteStateT, boost::random_access_traversal_tag>
 * that uses std::vector as an internal buffer.
 */

namespace kstate {

// Helper tag classes:
struct CtrFromRange {};
struct CtrFromBuffer {};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static CtrFromRange ctr_from_range{};
static CtrFromBuffer ctr_from_buffer{};
#pragma GCC diagnostic pop

// ***********************************************************************

// Helper functions:
template <typename Range>
std::vector<typename boost::range_value<Range>::type>
init_vector_from_range(
        const Range& rng) {
    using Element = typename boost::range_value<Range>::type;
    using Vector = std::vector<Element>;
    return Vector(std::begin(rng), std::end(rng));
}

// ***********************************************************************

template <typename _SiteStateTraitT>
struct DynamicKstateTypes {
    static_assert(IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    //    static_assert(!std::is_const<_SiteStateT>::value);
    //    static_assert(!std::is_volatile<_SiteStateT>::value);
    //    static_assert(!std::is_reference<_SiteStateT>::value);//TODO remove
    //    using SiteStateT = _SiteState;
    using BufferType = typename std::vector<SiteStateT>;
    using IteratorType = typename BufferType::iterator;
    using ConstIteratorType = typename BufferType::const_iterator;
    using RangeType = typename boost::iterator_range<IteratorType>;
    using ConstRangeType = typename boost::iterator_range<ConstIteratorType>;
    using AnyRangeType = typename boost::any_range<SiteStateT, boost::random_access_traversal_tag>;
    using ConstAnyRangeType = typename boost::any_range<const SiteStateT, boost::random_access_traversal_tag>;
};

template <typename _SiteStateTraitT>
class DynamicKstate : public SpeedyKstate<_SiteStateTraitT, typename DynamicKstateTypes<_SiteStateTraitT>::ConstRangeType> {
    static_assert(IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferType = typename DynamicKstateTypes<SiteStateTraitT>::BufferType;
    using IteratorType = typename DynamicKstateTypes<SiteStateTraitT>::IteratorType;
    using ConstIteratorType = typename DynamicKstateTypes<SiteStateTraitT>::ConstIteratorType;
    using RangeType = typename DynamicKstateTypes<SiteStateTraitT>::RangeType;
    using ConstRangeType = typename DynamicKstateTypes<SiteStateTraitT>::ConstRangeType;
    using AnyRangeType = typename DynamicKstateTypes<SiteStateTraitT>::AnyRangeType;
    using ConstAnyRangeType = typename DynamicKstateTypes<SiteStateTraitT>::ConstAnyRangeType;

public:
    DynamicKstate(BufferType&&, CtrFromBuffer);
    template <typename OtherRangeType>
    DynamicKstate(const OtherRangeType&, CtrFromRange);

public:
    ConstRangeType to_range() const override;
    size_t n_sites() const override;

protected:
    const BufferType _v;
};

// ***********************************************************************

template <typename _SiteStateTraitT>
DynamicKstate<_SiteStateTraitT>::DynamicKstate(
        DynamicKstate<_SiteStateTraitT>::BufferType&& v,
        CtrFromBuffer) :
    _v(std::move(v)) {
}

template <typename _SiteStateTraitT>
template <typename OtherRangeType>
DynamicKstate<_SiteStateTraitT>::DynamicKstate(const OtherRangeType& r, CtrFromRange) :
    _v(init_vector_from_range(r)) {
}

// ***********************************************************************

template <typename _SiteStateTraitT>
typename DynamicKstate<_SiteStateTraitT>::ConstRangeType
DynamicKstate<_SiteStateTraitT>::to_range() const {
    return _v;
}

template <typename _SiteStateTraitT>
size_t
DynamicKstate<_SiteStateTraitT>::n_sites() const {
    return _v.size();
}

}  // namespace kstate

// #######################################################################
// ## TraitsFor DynamicKstate                                           ##
// #######################################################################

namespace kstate {

template<typename _SiteStateTraitT>
struct TraitKstate<DynamicKstate<_SiteStateTraitT>> {
    static constexpr bool is_kstate_trait = true;
    using KstateT = DynamicKstate<_SiteStateTraitT>;
};

}

// #######################################################################
// ## StaticKstate                                                      ##
// #######################################################################

/*
 * DynamicKstate<SiteStateT, N> is a concrete subclass of
 * Kstate<SiteStateT, boost::random_access_traversal_tag>
 * that uses std::array as an internal buffer.
 */

namespace kstate {
//TODO
}

#endif
