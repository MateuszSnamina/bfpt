#ifndef KSTATE_KSTATE_ABSTRACT_HPP
#define KSTATE_KSTATE_ABSTRACT_HPP

#include <kstate/trait_site_state.hpp>

#include <extensions/adaptors.hpp>
#include <extensions/range_streamer.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm/search.hpp>
#include <boost/range/any_range.hpp>

#include <cassert>
#include <cmath>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

// #######################################################################
// ## Kstate                                                            ##
// #######################################################################

/*
 * Kstate<SiteStateT> class provides an abstraction for 1D cyclic quantum
 * state, such as the ground state of 1D Heisenberg spin chain (chain with
 * the periodic boundary condition imposed). The corresponding site state
 * is prescribed by SiteStateT class. (Following the example, SiteState
 * provides an abstraction for the spin Hilber space)
 *
 * A Kstate instance, together with $k$ value, define a Bloch quantum
 * state. ($k$ value is not a part of the instance.)
 *
 * The Kstate<SiteStateT> is an abstract base class.
 * Basing on the implementations the following two member functions:
 *   - to_any_range() const,
 *   - n_sites() const
 * the class provides an implementation of the
 * following member functions being the state descriptors:
 *   - n_least_replication_shift() const
 *   - norm_factor() const
 *   - is_prolific(int n_k) const
 *   - to_str() const
 * as well the following member functions allowing the states comparison:
 *   - compare_kstate(const Kstate<SiteStateT>& other) const
 *   - translational_compare_kstate(const Kstate<SiteStateT>& other) const.
 *   - compare_any_range(const ConstAnyRangeType& other) const;
 *   - translational_compare_any_range(const ConstAnyRangeType& other) const;
 *
 * The class may be fancy-formatted using KstateStreamer helper class.
 *
 * Kstate<SiteStateT> is the high-level API class, that relies on
 * polymorphic boost::any_range class. For implementation within
 * polymorphic-less framework check SpeedyKstate class.
 */

namespace kstate {

template <typename _SiteStateTraitT, typename _TraversalTagT = boost::random_access_traversal_tag>
class Kstate {
    static_assert(IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(std::is_same_v<_TraversalTagT, boost::random_access_traversal_tag> || std::is_same_v<_TraversalTagT, boost::forward_traversal_tag>);
public: // Helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using TraversalTagT = _TraversalTagT;
    using AnyRangeType = typename boost::any_range<SiteStateT, TraversalTagT>;
    using ConstAnyRangeType = typename boost::any_range<const SiteStateT, TraversalTagT>;
public:
    virtual ConstAnyRangeType to_any_range() const = 0;
    virtual size_t n_sites() const = 0;
public:  // Kstate descriptors:
    virtual size_t n_least_replication_shift() const;
    virtual double norm_factor() const;
    virtual bool is_prolific(int n_k) const;
    virtual std::string to_str() const;
public:  // Convenient binary functions:
    bool compare_kstate(const Kstate<SiteStateTraitT, TraversalTagT>& other) const;
    std::optional<size_t> translational_compare_kstate(const Kstate<SiteStateTraitT, TraversalTagT>& other) const;
    bool compare_any_range(const ConstAnyRangeType& other) const;
    std::optional<size_t> translational_compare_any_range(const ConstAnyRangeType& other) const;
public:
    virtual ~Kstate() = default;
};

// ***********************************************************************

template <typename _SiteStateTraitT, typename _TraversalTagT>
size_t
Kstate<_SiteStateTraitT, _TraversalTagT>::n_least_replication_shift() const {
    assert(n_sites() > 0);
    const auto r = to_any_range();
    const auto rdr = r | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    const auto it = boost::range::search(rdr, r);
    const auto _ = std::distance(std::begin(rdr), it);
    assert(_ >= 0);
    return static_cast<size_t>(_ + 1);
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
double
Kstate<_SiteStateTraitT, _TraversalTagT>::norm_factor() const {
    return std::sqrt(n_least_replication_shift()) / n_sites();
    // The result is equal to 1 / std::sqrt(n_least_replication_shift) / n_replicas;
    // where n_replicas = n_sites / n_least_replication_shift
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
bool
Kstate<_SiteStateTraitT, _TraversalTagT>::is_prolific(int n_k) const {
    return !((n_least_replication_shift() * n_k) % n_sites());
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
bool
Kstate<_SiteStateTraitT, _TraversalTagT>::compare_kstate(const Kstate<_SiteStateTraitT, _TraversalTagT>& other) const {
    return compare_any_range(other.to_any_range());
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
std::optional<size_t>
Kstate<_SiteStateTraitT, _TraversalTagT>::translational_compare_kstate(const Kstate<_SiteStateTraitT, _TraversalTagT>& other) const {
    return translational_compare_any_range(other.to_any_range());
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
bool
Kstate<_SiteStateTraitT, _TraversalTagT>::compare_any_range(const ConstAnyRangeType& other) const {
    assert(n_sites() == boost::size(other));
    return boost::range::equal(to_any_range(), other);
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
std::optional<size_t>
Kstate<_SiteStateTraitT, _TraversalTagT>::translational_compare_any_range(const ConstAnyRangeType& other) const {
    assert(n_sites() == boost::size(other));
    const auto& r1 = to_any_range();
    const auto& r2 = other;
    const auto r2d = r2 | extension::boost::adaptors::doubled;
    const auto it = boost::range::search(r2d, r1);
    return it == std::end(r2d)
            ? std::optional<size_t>()
            : static_cast<size_t>(std::distance(std::begin(r2d), it));
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
std::string
Kstate<_SiteStateTraitT, _TraversalTagT>::to_str() const {
    using namespace extension::boost::stream_pragma;
    const auto range_stream_settings = RSS<SiteStateT>()
            .set_string_preparer("⦃")
            .set_null_sustainer()
            .set_string_separer("∙")
            .set_string_finisher("⦄");
    return (to_any_range() | range_stream_settings).str();
}

}  // namespace kstate

// #######################################################################
// ## SpeedyKstate                                                      ##
// #######################################################################

/*
 * SpeedyKstate<ConstRangeType> class is to model the same physical
 * abstraction as Kstate<SiteStateT> class does. The two classes provide
 * alternative APIs and alternative implementations for a similar
 * functionality.
 *
 * The Kstate<SiteStateT> implementation relies a high level API, build
 * on top of the polymorphic abstraction layer (served by boost::any_range
 * class). The implementation may not be a desirable when execution speed
 * is the top priority.
 *
 * SpeedyKstate<ConstRangeType> is to take advantage of a polymorphic-less
 * API allowing high-speed implementations. The used approach relies
 * on ranges being treated as instances of template parameter classes,
 * rather than instances of boost::any_range class. Within the framework
 * more efficient implementations are possible.
 *
 * As the template-only approach used in SpeedyKstate<ConstRangeType>
 * may be perceived as a tuned alternative for the polymorphic approach
 * used in Kstate<SiteStateT> then the former is formally treated
 * as subclass of the latter.
 *
 * SpeedyKstate<ConstRangeType> overrides the following descriptor-type
 * member functions:
 *  - n_least_replication_shift() const,
 *  - norm_factor() const,
 *  - is_prolific(int n_k) const,
 *  - to_str() const.
 * And provides the following two member functions:
 *  - compare_range(const OtherConstRangeType& other) const,
 *  - translational_compare_range(const OtherConstRangeType& other) const,
 * being alternatives to:
 *   - compare_kstate(const Kstate<OtherSiteStateT>& other) const,
 *   - translational_compare_kstate(const Kstate<OtherSiteStateT>& other) const.
 *
 * Member function implementations defined in SpeedyKstate<ConstRangeType>
 * are based on to_range() member function. This is the origin of differences
 * between the implementations and their counterparts from Kstate<SiteStateT>,
 * as the latter may use only to_any_range() member function.
 *
 * The SpeedyKstate<ConstRangeType> is conceived to being the layer between
 * Kstate<SiteStateT> abstract base class and its concrete sub-classes, like
 *  - DynamicKstate<SiteStateT>, and
 *  - StaticKstate<SiteStateT, N>
 * classes.
 *
 */

namespace kstate {

template <typename _SiteStateTraitT, typename _ConstRangeType>
class SpeedyKstate : public Kstate<
        _SiteStateTraitT,
        typename boost::range_traversal<_ConstRangeType>::type> {
    static_assert(IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(std::is_same_v<typename _SiteStateTraitT::SiteStateT, typename boost::range_value<_ConstRangeType>::type>);
public:  // Helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using ConstRangeType = _ConstRangeType;
    using SiteStateT = typename boost::range_value<ConstRangeType>::type;
    using TraversalTagT = typename boost::range_traversal<ConstRangeType>::type;
    static_assert(!std::is_const<SiteStateT>::value);
    static_assert(!std::is_volatile<SiteStateT>::value);
    static_assert(!std::is_reference<SiteStateT>::value);
    static_assert(std::is_same<TraversalTagT, boost::random_access_traversal_tag>::value ||
    std::is_same<TraversalTagT, boost::forward_traversal_tag>::value);
    using AnyRangeType = typename Kstate<SiteStateTraitT, TraversalTagT>::AnyRangeType;
    using ConstAnyRangeType = typename Kstate<SiteStateTraitT, TraversalTagT>::ConstAnyRangeType;
public:  // Base on the two functions:
    virtual ConstRangeType to_range() const = 0;
    virtual size_t n_sites() const = 0;
public:  // SpeedyKstate implements pure virtual base class function:
    virtual ConstAnyRangeType to_any_range() const override;
public:  // SpeedyKstate overrides for fast implementations:
    size_t n_least_replication_shift() const override;
    double norm_factor() const override;
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

template <typename _SiteStateTraitT, typename _ConstRangeType>
typename SpeedyKstate<_SiteStateTraitT, _ConstRangeType>::ConstAnyRangeType
SpeedyKstate<_SiteStateTraitT, _ConstRangeType>::to_any_range() const {
    return to_range();
}

template <typename _SiteStateTraitT, typename _ConstRangeType>
size_t
SpeedyKstate<_SiteStateTraitT, _ConstRangeType>::n_least_replication_shift() const {
    assert(this->n_sites() > 0);
    const auto r = to_range();
    const auto rdr = r | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    const auto it = boost::range::search(rdr, r);
    const auto _ = std::distance(std::begin(rdr), it);
    assert(_ >= 0);
    return static_cast<size_t>(_ + 1);
}

template <typename _SiteStateTraitT, typename _ConstRangeType>
double
SpeedyKstate<_SiteStateTraitT, _ConstRangeType>::norm_factor() const {
    return std::sqrt(n_least_replication_shift()) / n_sites();
}

template <typename _SiteStateTraitT, typename _ConstRangeType>
bool
SpeedyKstate<_SiteStateTraitT, _ConstRangeType>::is_prolific(int n_k) const {
    return !((n_least_replication_shift() * n_k) % this->n_sites());
}

template <typename _SiteStateTraitT, typename _ConstRangeType>
template <typename OtherConstRangeType>
bool
SpeedyKstate<_SiteStateTraitT, _ConstRangeType>::compare_range(const OtherConstRangeType& other) const {
    assert(this->n_sites() == boost::size(other));
    return boost::range::equal(to_range(), other);
}

template <typename _SiteStateTraitT, typename _ConstRangeType>
template <typename OtherConstRangeType>
std::optional<size_t>
SpeedyKstate<_SiteStateTraitT, _ConstRangeType>::translational_compare_range(const OtherConstRangeType& other) const {
    assert(this->n_sites() == boost::size(other));
    const auto r1 = to_range();
    const auto r2 = other;
    const auto r2d = r2 | extension::boost::adaptors::doubled;
    const auto it = boost::range::search(r2d, r1);
    return it == std::end(r2d)
            ? std::optional<size_t>()
            : static_cast<size_t>(std::distance(std::begin(r2d), it));
}

template <typename _SiteStateTraitT, typename _ConstRangeType>
std::string
SpeedyKstate<_SiteStateTraitT, _ConstRangeType>::to_str() const {
    using namespace extension::boost::stream_pragma;
    const auto range_stream_settings = RSS<SiteStateT>()
            .set_string_preparer("⦃")
            .set_null_sustainer()
            .set_string_separer("∙")
            .set_string_finisher("⦄");
    return (to_range() | range_stream_settings).str();
}

}  // namespace kstate

#endif
