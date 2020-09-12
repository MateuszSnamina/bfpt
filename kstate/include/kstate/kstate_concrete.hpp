#ifndef KSTATE_KSTATE_CONCRETE_HPP
#define KSTATE_KSTATE_CONCRETE_HPP

#include <kstate/kstate_abstract.hpp>

#include <extensions/adaptors.hpp>
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
 * DynamicKstate<SiteType> is a concrete subclass of
 * Kstate<SiteType, boost::random_access_traversal_tag>
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

template <typename SiteType>
struct DynamicKstateTypes {
    static_assert(!std::is_const<SiteType>::value);
    static_assert(!std::is_volatile<SiteType>::value);
    static_assert(!std::is_reference<SiteType>::value);
    using BufferType = typename std::vector<SiteType>;
    using IteratorType = typename BufferType::iterator;
    using ConstIteratorType = typename BufferType::const_iterator;
    using RangeType = typename boost::iterator_range<IteratorType>;
    using ConstRangeType = typename boost::iterator_range<ConstIteratorType>;
    using AnyRangeType = typename boost::any_range<SiteType, boost::random_access_traversal_tag>;
    using ConstAnyRangeType = typename boost::any_range<const SiteType, boost::random_access_traversal_tag>;
};

template <typename SiteType>
class DynamicKstate : public SpeedyKstate<typename DynamicKstateTypes<SiteType>::ConstRangeType> {
    static_assert(!std::is_const<SiteType>::value);
    static_assert(!std::is_volatile<SiteType>::value);
    static_assert(!std::is_reference<SiteType>::value);

public:
    using BufferType = typename DynamicKstateTypes<SiteType>::BufferType;
    using IteratorType = typename DynamicKstateTypes<SiteType>::IteratorType;
    using ConstIteratorType = typename DynamicKstateTypes<SiteType>::ConstIteratorType;
    using RangeType = typename DynamicKstateTypes<SiteType>::RangeType;
    using ConstRangeType = typename DynamicKstateTypes<SiteType>::ConstRangeType;
    using AnyRangeType = typename DynamicKstateTypes<SiteType>::AnyRangeType;
    using ConstAnyRangeType = typename DynamicKstateTypes<SiteType>::ConstAnyRangeType;

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

template <typename SiteType>
DynamicKstate<SiteType>::DynamicKstate(DynamicKstate<SiteType>::BufferType&& v,
                                       CtrFromBuffer)
    : _v(std::move(v)) {
}

template <typename SiteType>
template <typename OtherRangeType>
DynamicKstate<SiteType>::DynamicKstate(const OtherRangeType& r, CtrFromRange)
    : _v(init_vector_from_range(r)) {
}

// ***********************************************************************

template <typename SiteType>
typename DynamicKstate<SiteType>::ConstRangeType
DynamicKstate<SiteType>::to_range() const {
    //boost::any_range<const SiteType, boost::random_access_traversal_tag> xxx(_v); // TODO ???
    return _v;
}

template <typename SiteType>
size_t
DynamicKstate<SiteType>::n_sites() const {
    return _v.size();
}

}  // namespace kstate

// #######################################################################
// ## UniqueDynamicKstate                                               ##
// #######################################################################

namespace kstate {

template <typename SiteType>
class DynamicUniqueKstate : public SpeedyKstate<typename DynamicKstateTypes<SiteType>::ConstRangeType> {

public:
    using BufferType = typename DynamicKstateTypes<SiteType>::BufferType;
    using IteratorType = typename DynamicKstateTypes<SiteType>::IteratorType;
    using ConstIteratorType = typename DynamicKstateTypes<SiteType>::ConstIteratorType;
    using RangeType = typename DynamicKstateTypes<SiteType>::RangeType;
    using ConstRangeType = typename DynamicKstateTypes<SiteType>::ConstRangeType;
    using AnyRangeType = typename DynamicKstateTypes<SiteType>::AnyRangeType;
    using ConstAnyRangeType = typename DynamicKstateTypes<SiteType>::ConstAnyRangeType;

public:
    template <typename SomeRangeType>
    DynamicUniqueKstate(const SomeRangeType& v, CtrFromRange);

public:
    ConstRangeType to_range() const override;
    size_t n_sites() const override;

protected:
    const BufferType _v;
};

// ***********************************************************************

template <typename SiteType>
template <typename SomeRangeType>
DynamicUniqueKstate<SiteType>::DynamicUniqueKstate(const SomeRangeType& r, CtrFromRange)
    : _v(init_vector_from_range(r | extension::boost::adaptors::rotated(n_unique_shift(r)))) {
}

// ***********************************************************************

template <typename SiteType>
typename DynamicUniqueKstate<SiteType>::ConstRangeType
DynamicUniqueKstate<SiteType>::to_range() const {
    return _v;
}

template <typename SiteType>
size_t
DynamicUniqueKstate<SiteType>::n_sites() const {
    return _v.size();
}

}  // namespace kstate

// #######################################################################
// ## StaticKstate                                                      ##
// #######################################################################

/*
 * DynamicKstate<SiteType, N> is a concrete subclass of
 * Kstate<SiteType, boost::random_access_traversal_tag>
 * that uses std::array as an internal buffer.
 */

namespace kstate {
//TODO
}

#endif
