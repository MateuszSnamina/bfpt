#ifndef KSTATE_KSTATE_HPP
#define KSTATE_KSTATE_HPP

#include <extensions/adaptors.hpp>
#include <extensions/range_streamer.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm/search.hpp>
#include <boost/range/any_range.hpp>

#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace {

template <typename ForwardRange>
size_t n_unique_shift(const ForwardRange& rng) {
    using Difference = typename boost::range_difference<ForwardRange>::type;
    const Difference d = std::distance(std::begin(rng), std::end(rng));
    size_t i = 0;
    for (size_t _ = 1; boost::numeric_cast<Difference>(_) < d; _++) {
        const bool result_cmp = boost::lexicographical_compare(
            rng | extension::boost::adaptors::rotated(i),
            rng | extension::boost::adaptors::rotated(_));
        if (result_cmp) {
            i = _;
        }
    }
    return i;
}

}  // namespace

namespace kstate {

// #######################################################################
// ## Kstate                                                            ##
// #######################################################################

template <typename SiteType, typename TraversalTag = boost::random_access_traversal_tag>
class Kstate {
   public:  // helper types:
    using AnyRangeType = typename boost::any_range<SiteType, TraversalTag>;
    using ConstAnyRangeType = typename boost::any_range<const SiteType, TraversalTag>;

   public:  // base on the two functions:
    virtual ConstAnyRangeType to_range() const = 0;
    virtual size_t n_sites() const = 0;

   public:  // Kstate implements:
    size_t n_least_replication_shift() const;
    bool is_prolific(int n_k) const;
    template <typename OtherConstRangeType>
    bool compare(const Kstate<OtherConstRangeType>& other) const;
    template <typename OtherConstRangeType>
    std::optional<size_t> tranlational_compare(const Kstate<OtherConstRangeType>& other) const;
    std::string to_str() const;

   public:
    virtual ~Kstate() = default;
};

// ***********************************************************************

template <typename SiteType, typename TraversalTag>
size_t Kstate<SiteType, TraversalTag>::n_least_replication_shift() const {
    assert(n_sites() > 0);
    const auto r = to_range();
    const auto rdr = r | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    const auto it = boost::range::search(rdr, r);
    const auto _ = std::distance(std::begin(rdr), it);
    assert(_ >= 0);
    return static_cast<size_t>(_ + 1);
}

template <typename SiteType, typename TraversalTag>
bool Kstate<SiteType, TraversalTag>::is_prolific(int n_k) const {
    return !((n_least_replication_shift() * n_k) % n_sites());
}

template <typename SiteType, typename TraversalTag>
template <typename OtherConstRangeType>
bool Kstate<SiteType, TraversalTag>::compare(
    const Kstate<OtherConstRangeType>& other) const {
    return boost::range::equal(to_range(), other.to_range());
}

template <typename SiteType, typename TraversalTag>
template <typename OtherConstRangeType>
std::optional<size_t>
Kstate<SiteType, TraversalTag>::tranlational_compare(
    const Kstate<OtherConstRangeType>& other) const {
    if (n_sites() != other.n_sites()) {
        return std::nullopt;
    }
    const auto r1 = to_range();
    const auto r2 = other.to_range();
    const auto r2d = r2 | extension::boost::adaptors::doubled;
    const auto it = boost::range::search(r2d, r1);
    return it == std::end(r2d)
               ? std::optional<size_t>()
               : static_cast<size_t>(std::distance(std::begin(r2d), it));
}

template <typename SiteType, typename TraversalTag>
std::string
Kstate<SiteType, TraversalTag>::to_str() const {
    return extension::boost::RangeStringStreamer()
        .set_stream_preparer([](std::ostream& s) { s << "⦃"; })
        .set_stream_sustainer([](std::ostream& s, size_t i) {})
        .set_stream_separer([](std::ostream& s) { s << "∙"; })
        .set_stream_finisher([](std::ostream& s) { s << "⦄"; })
        .stream(to_range())
        .str();
}

// #######################################################################
// ## KstateUniqueView                                                  ##
// #######################################################################

template <typename SiteType, typename TraversalTag = boost::random_access_traversal_tag>
class KstateUniqueView : public Kstate<SiteType, TraversalTag> {
   public:  // helper types:
    using AnyRangeType = typename boost::any_range<SiteType, TraversalTag>;
    using ConstAnyRangeType = typename boost::any_range<const SiteType, TraversalTag>;

   public:
    KstateUniqueView(const Kstate<SiteType, TraversalTag>&);

   public:
    ConstAnyRangeType to_range() const override;
    size_t n_sites() const override;

   private:
    const Kstate<SiteType, TraversalTag>& _r;
    const size_t _n_unique_shift;
};

// ***********************************************************************

template <typename SiteType, typename TraversalTag>
KstateUniqueView<SiteType, TraversalTag>::KstateUniqueView(const Kstate<SiteType, TraversalTag>& r)
    : _r(r),
      _n_unique_shift(n_unique_shift(r.to_range())) {
}

// ***********************************************************************

template <typename SiteType, typename TraversalTag>
typename KstateUniqueView<SiteType, TraversalTag>::ConstAnyRangeType
KstateUniqueView<SiteType, TraversalTag>::to_range() const {
    return _r.to_range() | extension::boost::adaptors::rotated(_n_unique_shift);
}

template <typename SiteType, typename TraversalTag>
size_t
KstateUniqueView<SiteType, TraversalTag>::n_sites() const {
    return _r.n_sites();
}

// #######################################################################
// ## DynamicKstate                                                     ##
// #######################################################################

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

/*
 * DynamicKstate uses std::vect as am internall buffer.
 */

template <typename SiteType>
struct DynamicKstateTypes {
    using BufferType = typename std::vector<SiteType>;
    using IteratorType = typename BufferType::iterator;
    using ConstIteratorType = typename BufferType::const_iterator;
    using RangeType = typename boost::iterator_range<IteratorType>;
    using ConstRangeType = typename boost::iterator_range<ConstIteratorType>;
    using AnyRangeType = typename boost::any_range<SiteType, boost::random_access_traversal_tag>;
    using ConstAnyRangeType = typename boost::any_range<const SiteType, boost::random_access_traversal_tag>;
};

template <typename SiteType>
class DynamicKstate : public Kstate<SiteType> {
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
    ConstAnyRangeType to_range() const override;
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
typename DynamicKstate<SiteType>::ConstAnyRangeType
DynamicKstate<SiteType>::to_range() const {
    boost::any_range<const SiteType, boost::random_access_traversal_tag> xxx(_v);
    return _v;
}

template <typename SiteType>
size_t
DynamicKstate<SiteType>::n_sites() const {
    return _v.size();
}

// #######################################################################
// ## UniqueDynamicKstate                                               ##
// #######################################################################

template <typename SiteType>
class DynamicUniqueKstate : public Kstate<SiteType> {
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
    ConstAnyRangeType to_range() const override;
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
typename DynamicUniqueKstate<SiteType>::ConstAnyRangeType
DynamicUniqueKstate<SiteType>::to_range() const {
    return _v;
}

template <typename SiteType>
size_t
DynamicUniqueKstate<SiteType>::n_sites() const {
    return _v.size();
}

// #######################################################################
// ## UniqueDynamicKstate                                               ##
// #######################################################################

// template <typename SiteType>
// class DynamicUniqueKstate
//     : public DynamicKstate<SiteType> {
//   using BufferType = typename DynamicKstateTypes<SiteType>::BufferType;
//   using IteratorType = typename DynamicKstateTypes<SiteType>::IteratorType;
//   using ConstIteratorType =
//       typename DynamicKstateTypes<SiteType>::ConstIteratorType;
//   using RangeType = typename DynamicKstateTypes<SiteType>::RangeType;
//   using ConstRangeType = typename
//   DynamicKstateTypes<SiteType>::ConstRangeType;

//  public:
//   template <typename SomeRangeType>
//   DynamicUniqueKstate(const SomeRangeType& v, CtrFromRange);

//  public:
//   using DynamicKstate<SiteType>::to_range;
//   using DynamicKstate<SiteType>::n_sites;
// };

// // ***********************************************************************

// template <typename SiteType>
// template <typename SomeRangeType>
// DynamicUniqueKstate<SiteType>::DynamicUniqueKstate(const SomeRangeType& r,
//                                                    CtrFromRange)
//     : DynamicKstate<SiteType>(
//           r | extension::boost::adaptors::rotated(n_unique_shift(r)),
//           ctr_from_range) {}

}  // namespace kstate

#endif