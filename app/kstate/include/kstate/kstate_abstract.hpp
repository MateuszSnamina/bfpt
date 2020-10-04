#ifndef KSTATE_KSTATE_ABSTRACT_HPP
#define KSTATE_KSTATE_ABSTRACT_HPP

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
 * Kstate<SiteState> class provides an abstraction for 1D cyclic quantum
 * state, such as the ground state of 1D Heisenberg spin chain (chain with
 * the periodic boundary condition imposed). The corresponding site state
 * is prescribed by SiteState class. (Following the example, SiteState
 * provides an abstraction for the spin Hilber space)
 *
 * A Kstate instance, together with $k$ value, define a Bloch quantum
 * state. ($k$ value is not a part of the instance.)
 *
 * The Kstate<SiteState> is an abstract base class.
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
 *   - compare_kstate(const Kstate<SiteState>& other) const
 *   - translational_compare_kstate(const Kstate<SiteState>& other) const.
 *   - compare_any_range(const ConstAnyRangeType& other) const;
 *   - translational_compare_any_range(const ConstAnyRangeType& other) const;
 *
 * The class may be fancy-formatted using KstateStreamer helper class.
 *
 * Kstate<SiteState> is the high-level API class, that relies on
 * polymorphic boost::any_range class. For implementation within
 * polymorphic-less framework check SpeedyKstate class.
 */

namespace kstate {

template <typename _SiteStateTrait, typename _TraversalTag = boost::random_access_traversal_tag>
class Kstate {
    static_assert(_SiteStateTrait::is_site_state_trait);
    static_assert(std::is_same_v<_TraversalTag, boost::random_access_traversal_tag> || std::is_same_v<_TraversalTag, boost::forward_traversal_tag>);
    //    static_assert(!std::is_array_v<_SiteState>);
    //    static_assert(!std::is_function_v<_SiteState>);
    //    static_assert(!std::is_void_v<std::decay<_SiteState>>);
    //    static_assert(!std::is_null_pointer_v<std::decay<_SiteState>>);
    //    static_assert(std::is_enum_v<std::decay<_SiteState>> || std::is_union_v<std::decay<_SiteState>> || std::is_class_v<std::decay<_SiteState>>);
    //    static_assert(!std::is_pointer_v<std::decay<_SiteState>>);
    //    static_assert(!std::is_member_object_pointer_v<_SiteState>);
    //    static_assert(!std::is_member_function_pointer_v<_SiteState>);
    //    static_assert(!std::is_const_v<_SiteState>);
    //    static_assert(!std::is_volatile_v<_SiteState>);
    //    static_assert(!std::is_reference_v<_SiteState>);//TODO remove
public: // Helper types:
    using SiteStateTrait = _SiteStateTrait;
    using SiteState = typename _SiteStateTrait::SiteStateT;
    using TraversalTag = _TraversalTag;
    using AnyRangeType = typename boost::any_range<SiteState, TraversalTag>;
    using ConstAnyRangeType = typename boost::any_range<const SiteState, TraversalTag>;

public:
    virtual ConstAnyRangeType to_any_range() const = 0;
    virtual size_t n_sites() const = 0;

public:  // Kstate descriptors:
    virtual size_t n_least_replication_shift() const;
    virtual double norm_factor() const;
    virtual bool is_prolific(int n_k) const;
    virtual std::string to_str() const;

public:  // Convenient binary functions:
    bool compare_kstate(const Kstate<SiteStateTrait, TraversalTag>& other) const;
    std::optional<size_t> translational_compare_kstate(const Kstate<SiteStateTrait, TraversalTag>& other) const;
    bool compare_any_range(const ConstAnyRangeType& other) const;
    std::optional<size_t> translational_compare_any_range(const ConstAnyRangeType& other) const;

public:
    virtual ~Kstate() = default;
};

// ***********************************************************************

template <typename _SiteStateTrait, typename _TraversalTag>
size_t
Kstate<_SiteStateTrait, _TraversalTag>::n_least_replication_shift() const {
    assert(n_sites() > 0);
    const auto r = to_any_range();
    const auto rdr = r | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    const auto it = boost::range::search(rdr, r);
    const auto _ = std::distance(std::begin(rdr), it);
    assert(_ >= 0);
    return static_cast<size_t>(_ + 1);
}

template <typename _SiteStateTrait, typename _TraversalTag>
double
Kstate<_SiteStateTrait, _TraversalTag>::norm_factor() const {
    return std::sqrt(n_least_replication_shift()) / n_sites();
    // The result is equal to 1 / std::sqrt(n_least_replication_shift) / n_replicas;
    // where n_replicas = n_sites / n_least_replication_shift
}

template <typename _SiteStateTrait, typename _TraversalTag>
bool
Kstate<_SiteStateTrait, _TraversalTag>::is_prolific(int n_k) const {
    return !((n_least_replication_shift() * n_k) % n_sites());
}

template <typename _SiteStateTrait, typename _TraversalTag>
bool
Kstate<_SiteStateTrait, _TraversalTag>::compare_kstate(const Kstate<_SiteStateTrait, _TraversalTag>& other) const {
    return compare_any_range(other.to_any_range());
}

template <typename _SiteStateTrait, typename _TraversalTag>
std::optional<size_t>
Kstate<_SiteStateTrait, _TraversalTag>::translational_compare_kstate(const Kstate<_SiteStateTrait, _TraversalTag>& other) const {
    return translational_compare_any_range(other.to_any_range());
}

template <typename _SiteStateTrait, typename _TraversalTag>
bool
Kstate<_SiteStateTrait, _TraversalTag>::compare_any_range(const ConstAnyRangeType& other) const {
    assert(n_sites() == boost::size(other));
    return boost::range::equal(to_any_range(), other);
}

template <typename _SiteStateTrait, typename _TraversalTag>
std::optional<size_t>
Kstate<_SiteStateTrait, _TraversalTag>::translational_compare_any_range(const ConstAnyRangeType& other) const {
    assert(n_sites() == boost::size(other));
    const auto& r1 = to_any_range();
    const auto& r2 = other;
    const auto r2d = r2 | extension::boost::adaptors::doubled;
    const auto it = boost::range::search(r2d, r1);
    return it == std::end(r2d)
            ? std::optional<size_t>()
            : static_cast<size_t>(std::distance(std::begin(r2d), it));
}

template <typename _SiteStateTrait, typename _TraversalTag>
std::string
Kstate<_SiteStateTrait, _TraversalTag>::to_str() const {
    using namespace extension::boost::stream_pragma;
    const auto range_stream_settings = RSS<SiteState>()
            .set_string_preparer("⦃")
            .set_null_sustainer()
            .set_string_separer("∙")
            .set_string_finisher("⦄");
    return (to_any_range() | range_stream_settings).str();
}

}  // namespace kstate

// #######################################################################
// ## KstateUniqueView                                                  ##
// #######################################################################

// template <typename _SiteState, typename _TraversalTag = boost::random_access_traversal_tag>
// class KstateUniqueView : public Kstate<_SiteState, _TraversalTag> {
//    public:  // Helper types:
//    using SiteState = _SiteState;
//    using TraversalTag = _TraversalTag;
//    using AnyRangeType = typename boost::any_range<SiteState, TraversalTag>;
//    using ConstAnyRangeType = typename boost::any_range<const SiteState, TraversalTag>;

//    public:
//     KstateUniqueView(const Kstate<SiteState, TraversalTag>&);

//    public:
//     ConstAnyRangeType to_range() const override;
//     size_t n_sites() const override;

//    private:
//     const Kstate<SiteState, TraversalTag>& _r;
//     const size_t _n_unique_shift;
// };

// // ***********************************************************************

// template <typename _SiteState, typename _TraversalTag>
// KstateUniqueView<_SiteState, _TraversalTag>::KstateUniqueView(const Kstate<SiteState, TraversalTag>& r)
//     : _r(r),
//       _n_unique_shift(n_unique_shift(r.to_range())) {
// }

// // ***********************************************************************

// template <typename _SiteState, typename _TraversalTag>
// typename KstateUniqueView<_SiteState, _TraversalTag>::ConstAnyRangeType
// KstateUniqueView<_SiteState, _TraversalTag>::to_range() const {
//     return _r.to_range() | extension::boost::adaptors::rotated(_n_unique_shift);
// }

// template <typename _SiteState, typename _TraversalTag>
// size_t
// KstateUniqueView<_SiteState, _TraversalTag>::n_sites() const {
//     return _r.n_sites();
// }

// #######################################################################
// ## SpeedyKstate                                                      ##
// #######################################################################

/*
 * SpeedyKstate<ConstRangeType> class is to model the same physical
 * abstraction as Kstate<SiteState> class does. The two classes provide
 * alternative APIs and alternative implementations for a similar
 * functionality.
 *
 * The Kstate<SiteState> implementation relies a high level API, build
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
 * used in Kstate<SiteState> then the former is formally treated
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
 *   - compare_kstate(const Kstate<OtherSiteState>& other) const,
 *   - translational_compare_kstate(const Kstate<OtherSiteState>& other) const.
 *
 * Member function implementations defined in SpeedyKstate<ConstRangeType>
 * are based on to_range() member function. This is the origin of differences
 * between the implementations and their counterparts from Kstate<SiteState>,
 * as the latter may use only to_any_range() member function.
 *
 * The SpeedyKstate<ConstRangeType> is conceived to being the layer between
 * Kstate<SiteState> abstract base class and its concrete sub-classes, like
 *  - DynamicKstate<SiteState>, and
 *  - StaticKstate<SiteState, N>
 * classes.
 *
 */

namespace kstate {

template <typename _SiteStateTrait, typename _ConstRangeType>
class SpeedyKstate : public Kstate<
        _SiteStateTrait,
        typename boost::range_traversal<_ConstRangeType>::type> {
    static_assert(_SiteStateTrait::is_site_state_trait);
    static_assert(std::is_same_v<typename _SiteStateTrait::SiteStateT, typename boost::range_value<_ConstRangeType>::type>);

public:  // Helper types:
    using SiteStateTrait = _SiteStateTrait;
    using ConstRangeType = _ConstRangeType;
    using SiteState = typename boost::range_value<ConstRangeType>::type;
    using TraversalTag = typename boost::range_traversal<ConstRangeType>::type;
    static_assert(!std::is_const<SiteState>::value);
    static_assert(!std::is_volatile<SiteState>::value);
    static_assert(!std::is_reference<SiteState>::value);
    static_assert(std::is_same<TraversalTag, boost::random_access_traversal_tag>::value ||
    std::is_same<TraversalTag, boost::forward_traversal_tag>::value);
    using AnyRangeType = typename Kstate<SiteStateTrait, TraversalTag>::AnyRangeType;
    using ConstAnyRangeType = typename Kstate<SiteStateTrait, TraversalTag>::ConstAnyRangeType;

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

template <typename _SiteStateTrait, typename _ConstRangeType>
typename SpeedyKstate<_SiteStateTrait, _ConstRangeType>::ConstAnyRangeType
SpeedyKstate<_SiteStateTrait, _ConstRangeType>::to_any_range() const {
    return to_range();
}

template <typename _SiteStateTrait, typename _ConstRangeType>
size_t
SpeedyKstate<_SiteStateTrait, _ConstRangeType>::n_least_replication_shift() const {
    assert(this->n_sites() > 0);
    const auto r = to_range();
    const auto rdr = r | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    const auto it = boost::range::search(rdr, r);
    const auto _ = std::distance(std::begin(rdr), it);
    assert(_ >= 0);
    return static_cast<size_t>(_ + 1);
}

template <typename _SiteStateTrait, typename _ConstRangeType>
double
SpeedyKstate<_SiteStateTrait, _ConstRangeType>::norm_factor() const {
    return std::sqrt(n_least_replication_shift()) / n_sites();
}

template <typename _SiteStateTrait, typename _ConstRangeType>
bool
SpeedyKstate<_SiteStateTrait, _ConstRangeType>::is_prolific(int n_k) const {
    return !((n_least_replication_shift() * n_k) % this->n_sites());
}

template <typename _SiteStateTrait, typename _ConstRangeType>
template <typename OtherConstRangeType>
bool
SpeedyKstate<_SiteStateTrait, _ConstRangeType>::compare_range(const OtherConstRangeType& other) const {
    assert(this->n_sites() == boost::size(other));
    return boost::range::equal(to_range(), other);
}

template <typename _SiteStateTrait, typename _ConstRangeType>
template <typename OtherConstRangeType>
std::optional<size_t>
SpeedyKstate<_SiteStateTrait, _ConstRangeType>::translational_compare_range(const OtherConstRangeType& other) const {
    assert(this->n_sites() == boost::size(other));
    const auto r1 = to_range();
    const auto r2 = other;
    const auto r2d = r2 | extension::boost::adaptors::doubled;
    const auto it = boost::range::search(r2d, r1);
    return it == std::end(r2d)
            ? std::optional<size_t>()
            : static_cast<size_t>(std::distance(std::begin(r2d), it));
}

template <typename _SiteStateTrait, typename _ConstRangeType>
std::string
SpeedyKstate<_SiteStateTrait, _ConstRangeType>::to_str() const {
    using namespace extension::boost::stream_pragma;
    const auto range_stream_settings = RSS<SiteState>()
            .set_string_preparer("⦃")
            .set_null_sustainer()
            .set_string_separer("∙")
            .set_string_finisher("⦄");
    return (to_range() | range_stream_settings).str();
}

}  // namespace kstate

#endif
