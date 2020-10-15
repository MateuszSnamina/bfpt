#pragma once

#include <kstate_op_range/op_range_compare.hpp>
#include <kstate_op_range/op_range_least_replication_shift.hpp>

#include <kstate_trait/trait_site_state.hpp>

#include <extensions/adaptors.hpp>
#include <extensions/range_streamer.hpp>

#include <boost/algorithm/string/predicate.hpp>//TODO remove
#include <boost/numeric/conversion/cast.hpp>
#include <boost/range.hpp>
#include <boost/range/algorithm/search.hpp>//TODO remove
#include <boost/range/any_range.hpp>

#include <cassert>
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
 *   - compare_equality_kstate(const Kstate<SiteStateT>& other) const
 *   - compare_translational_equality_kstate(const Kstate<SiteStateT>& other) const.
 *   - compare_equality_any_range(const ConstAnyRangeT& other) const;
 *   - compare_translational_equality_any_range(const ConstAnyRangeT& other) const;
 *
 * The class may be fancy-formatted using KstateStreamer helper class.
 *
 * Kstate<SiteStateT> is the high-level API class, that relies on
 * polymorphic boost::any_range class. For implementation within
 * polymorphic-less framework check SpeedyKstate class.
 */

namespace kstate_impl {

template <typename _SiteStateTraitT, typename _TraversalTagT = boost::random_access_traversal_tag>
class Kstate {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(std::is_same_v<_TraversalTagT, boost::random_access_traversal_tag> || std::is_same_v<_TraversalTagT, boost::forward_traversal_tag>);
public: // Helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using TraversalTagT = _TraversalTagT;
    using AnyRangeT = typename boost::any_range<SiteStateT, TraversalTagT>;
    using ConstAnyRangeT = typename boost::any_range<const SiteStateT, TraversalTagT>;
public:
    virtual ConstAnyRangeT to_any_range() const noexcept = 0;
    virtual size_t n_sites() const noexcept = 0;
public:  // Kstate descriptors:
    virtual size_t n_least_replication_shift() const noexcept ;
    virtual double norm_factor() const noexcept ;
    virtual bool is_prolific(int n_k) const noexcept ;
    virtual std::string to_str() const noexcept ;
public:  // Convenient binary functions:
    bool compare_equality_kstate(const Kstate<SiteStateTraitT, TraversalTagT>& other) const noexcept ;
    std::optional<size_t> compare_translational_equality_kstate(const Kstate<SiteStateTraitT, TraversalTagT>& other) const noexcept ;
    bool compare_equality_any_range(const ConstAnyRangeT& other) const noexcept ;
    std::optional<size_t> compare_translational_equality_any_range(const ConstAnyRangeT& other) const noexcept ;
public:
    virtual ~Kstate() = default;
};

// ***********************************************************************

template <typename _SiteStateTraitT, typename _TraversalTagT>
size_t
Kstate<_SiteStateTraitT, _TraversalTagT>::n_least_replication_shift() const noexcept {
    return kstate_op_range::n_least_replication_shift(to_any_range());
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
double
Kstate<_SiteStateTraitT, _TraversalTagT>::norm_factor() const noexcept {
    return kstate_op_range::norm_factor(to_any_range());
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
bool
Kstate<_SiteStateTraitT, _TraversalTagT>::is_prolific(int n_k) const noexcept {
    return kstate_op_range::is_prolific(to_any_range(), n_k);
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
bool
Kstate<_SiteStateTraitT, _TraversalTagT>::compare_equality_kstate(const Kstate<_SiteStateTraitT, _TraversalTagT>& other) const noexcept {
    return compare_equality_any_range(other.to_any_range());
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
std::optional<size_t>
Kstate<_SiteStateTraitT, _TraversalTagT>::compare_translational_equality_kstate(const Kstate<_SiteStateTraitT, _TraversalTagT>& other) const noexcept {
    return compare_translational_equality_any_range(other.to_any_range());
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
bool
Kstate<_SiteStateTraitT, _TraversalTagT>::compare_equality_any_range(const ConstAnyRangeT& other) const noexcept {
    assert(n_sites() == boost::size(other));
    return kstate_op_range::compare_equality(to_any_range(), other);
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
std::optional<size_t>
Kstate<_SiteStateTraitT, _TraversalTagT>::compare_translational_equality_any_range(const ConstAnyRangeT& other) const noexcept {
    assert(this->n_sites() == boost::size(other));
    return kstate_op_range::compare_translational_equality(to_any_range(), other);
}

template <typename _SiteStateTraitT, typename _TraversalTagT>
std::string
Kstate<_SiteStateTraitT, _TraversalTagT>::to_str() const noexcept {
    using namespace extension::boost::stream_pragma;
    const auto range_stream_settings = RSS<SiteStateT>()
            .set_string_preparer("⦃")
            .set_null_sustainer()
            .set_string_separer("∙")
            .set_string_finisher("⦄");
    return (to_any_range() | range_stream_settings).str();
}

}  // namespace kstate_impl

// #######################################################################
// ## SpeedyKstate                                                      ##
// #######################################################################

/*
 * SpeedyKstate<ConstRangeT> class is to model the same physical
 * abstraction as Kstate<SiteStateT> class does. The two classes provide
 * alternative APIs and alternative implementations for a similar
 * functionality.
 *
 * The Kstate<SiteStateT> implementation relies a high level API, build
 * on top of the polymorphic abstraction layer (served by boost::any_range
 * class). The implementation may not be a desirable when execution speed
 * is the top priority.
 *
 * SpeedyKstate<ConstRangeT> is to take advantage of a polymorphic-less
 * API allowing high-speed implementations. The used approach relies
 * on ranges being treated as instances of template parameter classes,
 * rather than instances of boost::any_range class. Within the framework
 * more efficient implementations are possible.
 *
 * As the template-only approach used in SpeedyKstate<ConstRangeT>
 * may be perceived as a tuned alternative for the polymorphic approach
 * used in Kstate<SiteStateT> then the former is formally treated
 * as subclass of the latter.
 *
 * SpeedyKstate<ConstRangeT> overrides the following descriptor-type
 * member functions:
 *  - n_least_replication_shift() const,
 *  - norm_factor() const,
 *  - is_prolific(int n_k) const,
 *  - to_str() const.
 * And provides the following two member functions:
 *  - compare_equality_range(const OtherConstRangeT& other) const,
 *  - compare_translational_equality_range(const OtherConstRangeT& other) const,
 * being alternatives to:
 *   - compare_equality_kstate(const Kstate<OtherSiteStateT>& other) const,
 *   - compare_translational_equality_kstate(const Kstate<OtherSiteStateT>& other) const.
 *
 * Member function implementations defined in SpeedyKstate<ConstRangeT>
 * are based on to_range() member function. This is the origin of differences
 * between the implementations and their counterparts from Kstate<SiteStateT>,
 * as the latter may use only to_any_range() member function.
 *
 * The SpeedyKstate<ConstRangeT> is conceived to being the layer between
 * Kstate<SiteStateT> abstract base class and its concrete sub-classes, like
 *  - DynamicKstate<SiteStateT>, and
 *  - StaticKstate<SiteStateT, N>
 * classes.
 *
 */

namespace kstate_impl {

template <typename _SiteStateTraitT, typename _ConstRangeT>
class SpeedyKstate : public Kstate<
        _SiteStateTraitT,
        typename boost::range_traversal<_ConstRangeT>::type> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(std::is_same_v<typename _SiteStateTraitT::SiteStateT, typename boost::range_value<_ConstRangeT>::type>);
public:  // Helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using ConstRangeT = _ConstRangeT;
    using SiteStateT = typename boost::range_value<ConstRangeT>::type;
    using TraversalTagT = typename boost::range_traversal<ConstRangeT>::type;
    static_assert(!std::is_const<SiteStateT>::value);
    static_assert(!std::is_volatile<SiteStateT>::value);
    static_assert(!std::is_reference<SiteStateT>::value);
    static_assert(std::is_same<TraversalTagT, boost::random_access_traversal_tag>::value ||
    std::is_same<TraversalTagT, boost::forward_traversal_tag>::value);
    using AnyRangeT = typename Kstate<SiteStateTraitT, TraversalTagT>::AnyRangeT;
    using ConstAnyRangeT = typename Kstate<SiteStateTraitT, TraversalTagT>::ConstAnyRangeT;
public:  // Base on the two functions:
    virtual ConstRangeT to_range() const noexcept = 0;
    virtual size_t n_sites() const noexcept = 0;
public:  // SpeedyKstate implements pure virtual base class function:
    virtual ConstAnyRangeT to_any_range() const noexcept override;
public:  // SpeedyKstate overrides for fast implementations:
    size_t n_least_replication_shift() const noexcept override;
    double norm_factor() const noexcept override;
    bool is_prolific(int n_k) const noexcept override;
    template <typename OtherConstRangeT>
    bool compare_equality_range(const OtherConstRangeT& other) const noexcept;
    template <typename OtherConstRangeT>
    std::optional<size_t> compare_translational_equality_range(const OtherConstRangeT& other) const noexcept ;
    std::string to_str() const noexcept override;
public:
    virtual ~SpeedyKstate() = default;
};

// ***********************************************************************

template <typename _SiteStateTraitT, typename _ConstRangeT>
typename SpeedyKstate<_SiteStateTraitT, _ConstRangeT>::ConstAnyRangeT
SpeedyKstate<_SiteStateTraitT, _ConstRangeT>::to_any_range() const noexcept {
    return to_range();
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
size_t
SpeedyKstate<_SiteStateTraitT, _ConstRangeT>::n_least_replication_shift() const noexcept {
    return kstate_op_range::n_least_replication_shift(to_range());
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
double
SpeedyKstate<_SiteStateTraitT, _ConstRangeT>::norm_factor() const noexcept {
    return kstate_op_range::norm_factor(to_range());
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
bool
SpeedyKstate<_SiteStateTraitT, _ConstRangeT>::is_prolific(int n_k) const noexcept {
    return kstate_op_range::is_prolific(to_range(), n_k);
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
template <typename OtherConstRangeT>
bool
SpeedyKstate<_SiteStateTraitT, _ConstRangeT>::compare_equality_range(const OtherConstRangeT& other) const noexcept {
    assert(this->n_sites() == boost::size(other));
    return kstate_op_range::compare_equality(to_range(), other);
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
template <typename OtherConstRangeT>
std::optional<size_t>
SpeedyKstate<_SiteStateTraitT, _ConstRangeT>::compare_translational_equality_range(const OtherConstRangeT& other) const noexcept {
    //assert(this->n_sites() == boost::size(other));
    return kstate_op_range::compare_translational_equality(to_range(), other);
}

template <typename _SiteStateTraitT, typename _ConstRangeT>
std::string
SpeedyKstate<_SiteStateTraitT, _ConstRangeT>::to_str() const noexcept {
    using namespace extension::boost::stream_pragma;
    const auto range_stream_settings = RSS<SiteStateT>()
            .set_string_preparer("⦃")
            .set_null_sustainer()
            .set_string_separer("∙")
            .set_string_finisher("⦄");
    return (to_range() | range_stream_settings).str();
}

}  // namespace kstate_impl
