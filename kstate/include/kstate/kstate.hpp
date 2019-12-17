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

/*
 * Kstate<SiteType> class provides an abstraction for 1D cyclic quantum
 * state, such as the ground state of 1D Heisenberg spin chain (chain with
 * the periodic boundary condition imposed). The corresponding site state
 * is prescribed by SiteType class. (Following the example, SiteType
 * provides an abstraction for the spin Hilber space)
 * 
 * A Kstate instance, together with $k$ value, define a Bloch quantum
 * state. ($k$ value is not a part of the instance.)
 * 
 * The Kstate<SiteType> is an abstract base class.
 * Basing on the implementations the following two member functions:
 *   - to_any_range() const,
 *   - n_sites() const
 * the class provides an implementation of the
 * following member functions being the state descriptors:
 *   - n_least_replication_shift() const
 *   - is_prolific(int n_k) const
 *   - to_str() const
 * as well the following member functions allowing the states comparison:
 *   - compare_kstate(const Kstate<OtherSiteType>& other) const
 *   - translational_compare_kstate(const Kstate<OtherSiteType>& other) const.
 * 
 * The class may be fancy-formatted using KstateStreamer helper class.
 * 
 * Kstate<SiteType> is the high-level API class, that relies on
 * polymorphic boost::any_range class. For implementation within
 * polymorphic-less framework check SpeedyKstate class.
 */

template <typename SiteType, typename TraversalTag = boost::random_access_traversal_tag>
class Kstate {
   public:  // helper types:
    using AnyRangeType = typename boost::any_range<SiteType, TraversalTag>;
    using ConstAnyRangeType = typename boost::any_range<const SiteType, TraversalTag>;

   public:
    virtual ConstAnyRangeType to_any_range() const = 0;
    virtual size_t n_sites() const = 0;

   public:  // Kstate descriptors:
    virtual size_t n_least_replication_shift() const;
    virtual bool is_prolific(int n_k) const;
    virtual std::string to_str() const;

   public:  // Convenient binary functions:
    template <typename OtherSiteType>
    bool compare_kstate(const Kstate<OtherSiteType>& other) const;
    template <typename OtherSiteType>
    std::optional<size_t> translational_compare_kstate(const Kstate<OtherSiteType>& other) const;

   public:
    virtual ~Kstate() = default;
};

// ***********************************************************************

template <typename SiteType, typename TraversalTag>
size_t Kstate<SiteType, TraversalTag>::n_least_replication_shift() const {
    assert(n_sites() > 0);
    const auto r = to_any_range();
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
template <typename OtherSiteType>
bool Kstate<SiteType, TraversalTag>::compare_kstate(const Kstate<OtherSiteType>& other) const {
    assert(n_sites() == other.n_sites());
    return boost::range::equal(to_any_range(), other.to_any_range());
}

template <typename SiteType, typename TraversalTag>
template <typename OtherSiteType>
std::optional<size_t>
Kstate<SiteType, TraversalTag>::translational_compare_kstate(const Kstate<OtherSiteType>& other) const {
    assert(n_sites() == other.n_sites());
    const auto r1 = to_any_range();
    const auto r2 = other.to_any_range();
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
        .stream(to_any_range())
        .str();
}

// #######################################################################
// ## KstateUniqueView                                                  ##
// #######################################################################

// template <typename SiteType, typename TraversalTag = boost::random_access_traversal_tag>
// class KstateUniqueView : public Kstate<SiteType, TraversalTag> {
//    public:  // helper types:
//     using AnyRangeType = typename boost::any_range<SiteType, TraversalTag>;
//     using ConstAnyRangeType = typename boost::any_range<const SiteType, TraversalTag>;

//    public:
//     KstateUniqueView(const Kstate<SiteType, TraversalTag>&);

//    public:
//     ConstAnyRangeType to_range() const override;
//     size_t n_sites() const override;

//    private:
//     const Kstate<SiteType, TraversalTag>& _r;
//     const size_t _n_unique_shift;
// };

// // ***********************************************************************

// template <typename SiteType, typename TraversalTag>
// KstateUniqueView<SiteType, TraversalTag>::KstateUniqueView(const Kstate<SiteType, TraversalTag>& r)
//     : _r(r),
//       _n_unique_shift(n_unique_shift(r.to_range())) {
// }

// // ***********************************************************************

// template <typename SiteType, typename TraversalTag>
// typename KstateUniqueView<SiteType, TraversalTag>::ConstAnyRangeType
// KstateUniqueView<SiteType, TraversalTag>::to_range() const {
//     return _r.to_range() | extension::boost::adaptors::rotated(_n_unique_shift);
// }

// template <typename SiteType, typename TraversalTag>
// size_t
// KstateUniqueView<SiteType, TraversalTag>::n_sites() const {
//     return _r.n_sites();
// }

// #######################################################################
// ## SpeedyKstate                                                      ##
// #######################################################################

/*
 * SpeedyKstate<ConstRangeType> class is to model the same physical
 * abstraction as Kstate<SiteType> class does. The two classes provide
 * alternative APIs and alternative implementations for a similar
 * functionality.
 * 
 * The Kstate<SiteType> implementation relies a high level API, build
 * on top of the polymorphic abstraction layer (served by boost::any_range
 * class). The implementation may not be a desirable when execution speed
 * is the top priority.
 *
 * SpeedyKstate<ConstRangeType> is to take advantage of a polymorphic-less
 * API allowing high-speed implementations. The used approach relies
 * on ranges being treated as instances of template parameter classes,
 * rather than instances of boost::any_range class. Within the framework
 * more efficient implementations are possible.
 * 
 * As the template-only approach used in SpeedyKstate<ConstRangeType>
 * may be perceived as a tuned alternative for the polymorphic approach
 * used in Kstate<SiteType> then the former is formally treated
 * as subclass of the latter.
 * 
 * SpeedyKstate<ConstRangeType> overrides the following descriptor-type
 * member functions:
 *  - n_least_replication_shift() const.
 *  - is_prolific(int n_k) const.
 *  - to_str() const.
 * And provides the following two member functions:
 *  - compare_range(const OtherConstRangeType& other) const,
 *  - translational_compare_range(const OtherConstRangeType& other) const,
 * being alternatives to:
 *   - compare_kstate(const Kstate<OtherSiteType>& other) const,
 *   - translational_compare_kstate(const Kstate<OtherSiteType>& other) const.
 * 
 * Member function implementations defined in SpeedyKstate<ConstRangeType>
 * are based on to_range() member function. This is the origin of differences
 * between the implementations and their counterparts from Kstate<SiteType>,
 * as the latter may use only to_any_range() member function.
 * 
 * The SpeedyKstate<ConstRangeType> is conceived to being the layer between
 * Kstate<SiteType> abstract base class and its concrete sub-classes, like
 *  - DynamicKstate<SiteType>, and
 *  - StaticKstate<SiteType, N>
 * classes.
 * 
 */

template <typename ConstRangeType>
class SpeedyKstate : public Kstate<typename boost::range_value<ConstRangeType>::type, typename boost::range_traversal<ConstRangeType>::type> {
   public:  // helper types:
    using SiteType = typename boost::range_value<ConstRangeType>::type;
    using TraversalTag = typename boost::range_traversal<ConstRangeType>::type;
    using AnyRangeType = typename Kstate<SiteType, TraversalTag>::AnyRangeType;
    using ConstAnyRangeType = typename Kstate<SiteType, TraversalTag>::ConstAnyRangeType;

   public:  // base on the two functions:
    virtual ConstRangeType to_range() const = 0;
    virtual size_t n_sites() const = 0;

   public:  // Kstate implements:
    virtual ConstAnyRangeType to_any_range() const override;

   public:  // Fast implementations implements:
    size_t n_least_replication_shift() const override;
    bool is_prolific(int n_k) const override;
    template <typename OtherConstRangeType>
    bool compare_range(const OtherConstRangeType& other) const;
    template <typename OtherConstRangeType>
    std::optional<size_t> translational_compare_range(const OtherConstRangeType& other) const;
    std::string to_str() const override;

   public:
    virtual ~SpeedyKstate() = default;
};

// ***********************************************************************

template <typename ConstRangeType>
typename SpeedyKstate<ConstRangeType>::ConstAnyRangeType
SpeedyKstate<ConstRangeType>::to_any_range() const {
    return to_range();
}

template <typename ConstRangeType>
size_t SpeedyKstate<ConstRangeType>::n_least_replication_shift() const {
    assert(this->n_sites() > 0);
    const auto r = to_range();
    const auto rdr = r | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    const auto it = boost::range::search(rdr, r);
    const auto _ = std::distance(std::begin(rdr), it);
    assert(_ >= 0);
    return static_cast<size_t>(_ + 1);
}

template <typename ConstRangeType>
bool SpeedyKstate<ConstRangeType>::is_prolific(int n_k) const {
    return !((n_least_replication_shift() * n_k) % this->n_sites());
}

template <typename ConstRangeType>
template <typename OtherConstRangeType>
bool SpeedyKstate<ConstRangeType>::compare_range(const OtherConstRangeType& other) const {
    assert(this->n_sites() == boost::size(other));
    return boost::range::equal(to_range(), other);
}

template <typename ConstRangeType>
template <typename OtherConstRangeType>
std::optional<size_t>
SpeedyKstate<ConstRangeType>::translational_compare_range(const OtherConstRangeType& other) const {
    assert(this->n_sites() == boost::size(other));
    const auto r1 = to_range();
    const auto r2 = other;
    const auto r2d = r2 | extension::boost::adaptors::doubled;
    const auto it = boost::range::search(r2d, r1);
    return it == std::end(r2d)
               ? std::optional<size_t>()
               : static_cast<size_t>(std::distance(std::begin(r2d), it));
}

template <typename ConstRangeType>
std::string
SpeedyKstate<ConstRangeType>::to_str() const {
    return extension::boost::RangeStringStreamer()
        .set_stream_preparer([](std::ostream& s) { s << "⦃"; })
        .set_stream_sustainer([](std::ostream& s, size_t i) {})
        .set_stream_separer([](std::ostream& s) { s << "∙"; })
        .set_stream_finisher([](std::ostream& s) { s << "⦄"; })
        .stream(to_range())
        .str();
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
class DynamicKstate : public SpeedyKstate<typename DynamicKstateTypes<SiteType>::ConstRangeType> {
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
class DynamicUniqueKstate : public SpeedyKstate<typename DynamicKstateTypes<SiteType>::ConstRangeType> {
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