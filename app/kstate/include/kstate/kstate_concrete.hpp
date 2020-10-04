#ifndef KSTATE_KSTATE_CONCRETE_HPP
#define KSTATE_KSTATE_CONCRETE_HPP

#include <kstate/kstate_abstract.hpp>
#include <kstate/trait_kstate.hpp>

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

template <typename _SiteStateTrait>
struct DynamicKstateTypes {
    static_assert(_SiteStateTrait::is_site_state_trait);
    using SiteStateTrait = _SiteStateTrait;
    using SiteStateT = typename _SiteStateTrait::SiteStateT;
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

template <typename _SiteStateTrait>
class DynamicKstate : public SpeedyKstate<_SiteStateTrait, typename DynamicKstateTypes<_SiteStateTrait>::ConstRangeType> {
    //    static_assert(!std::is_const<_SiteStateT>::value);
    //    static_assert(!std::is_volatile<_SiteStateT>::value);
    //    static_assert(!std::is_reference<_SiteStateT>::value);//TODO remove

public:
    using SiteStateTrait = _SiteStateTrait;
    using SiteStateT = typename _SiteStateTrait::SiteStateT;
    using BufferType = typename DynamicKstateTypes<SiteStateTrait>::BufferType;
    using IteratorType = typename DynamicKstateTypes<SiteStateTrait>::IteratorType;
    using ConstIteratorType = typename DynamicKstateTypes<SiteStateTrait>::ConstIteratorType;
    using RangeType = typename DynamicKstateTypes<SiteStateTrait>::RangeType;
    using ConstRangeType = typename DynamicKstateTypes<SiteStateTrait>::ConstRangeType;
    using AnyRangeType = typename DynamicKstateTypes<SiteStateTrait>::AnyRangeType;
    using ConstAnyRangeType = typename DynamicKstateTypes<SiteStateTrait>::ConstAnyRangeType;

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

template <typename _SiteStateTrait>
DynamicKstate<_SiteStateTrait>::DynamicKstate(
        DynamicKstate<_SiteStateTrait>::BufferType&& v,
        CtrFromBuffer) :
    _v(std::move(v)) {
}

template <typename _SiteStateTrait>
template <typename OtherRangeType>
DynamicKstate<_SiteStateTrait>::DynamicKstate(const OtherRangeType& r, CtrFromRange) :
    _v(init_vector_from_range(r)) {
}

// ***********************************************************************

template <typename _SiteStateTrait>
typename DynamicKstate<_SiteStateTrait>::ConstRangeType
DynamicKstate<_SiteStateTrait>::to_range() const {
    return _v;
}

template <typename _SiteStateTrait>
size_t
DynamicKstate<_SiteStateTrait>::n_sites() const {
    return _v.size();
}

}  // namespace kstate

// #######################################################################
// ## TraitsFor DynamicKstate                                           ##
// #######################################################################

namespace kstate {

template<typename _SiteStateTrait>
struct TraitKstate<DynamicKstate<_SiteStateTrait>> {
    static constexpr bool is_kstate_trait = true;
    using KstateT = DynamicKstate<_SiteStateTrait>;
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
