#pragma once

#include <kstate_impl/kstate_abstract.hpp>

#include <kstate_trait/trait_site_state.hpp>
#include <kstate_trait/trait_kstate.hpp>

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
#include <memory>

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
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    //    static_assert(!std::is_const<_SiteStateT>::value);
    //    static_assert(!std::is_volatile<_SiteStateT>::value);
    //    static_assert(!std::is_reference<_SiteStateT>::value);//TODO remove
    //    using SiteStateT = _SiteState;
    using BufferT = typename std::vector<SiteStateT>;
    using IteratorT = typename BufferT::iterator;
    using ConstIteratorT = typename BufferT::const_iterator;
    using RangeT = typename boost::iterator_range<IteratorT>;
    using ConstRangeT = typename boost::iterator_range<ConstIteratorT>;
    using AnyRangeT = typename boost::any_range<SiteStateT, boost::random_access_traversal_tag>;
    using ConstAnyRangeT = typename boost::any_range<const SiteStateT, boost::random_access_traversal_tag>;
};

template <typename _SiteStateTraitT>
class DynamicKstate : public SpeedyKstate<_SiteStateTraitT, typename DynamicKstateTypes<_SiteStateTraitT>::ConstRangeT> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = typename DynamicKstateTypes<SiteStateTraitT>::BufferT;
    using IteratorT = typename DynamicKstateTypes<SiteStateTraitT>::IteratorT;
    using ConstIteratorT = typename DynamicKstateTypes<SiteStateTraitT>::ConstIteratorT;
    using RangeT = typename DynamicKstateTypes<SiteStateTraitT>::RangeT;
    using ConstRangeT = typename DynamicKstateTypes<SiteStateTraitT>::ConstRangeT;
    using AnyRangeT = typename DynamicKstateTypes<SiteStateTraitT>::AnyRangeT;
    using ConstAnyRangeT = typename DynamicKstateTypes<SiteStateTraitT>::ConstAnyRangeT;

public:
    DynamicKstate(BufferT&&, CtrFromBuffer);
    template <typename OtherRangeT>
    DynamicKstate(const OtherRangeT&, CtrFromRange);

public:
    ConstRangeT to_range() const override;
    size_t n_sites() const override;

protected:
    const BufferT _v;
};

// ***********************************************************************

template <typename _SiteStateTraitT>
DynamicKstate<_SiteStateTraitT>::DynamicKstate(
        DynamicKstate<_SiteStateTraitT>::BufferT&& v,
        CtrFromBuffer) :
    _v(std::move(v)) {
}

template <typename _SiteStateTraitT>
template <typename OtherRangeT>
DynamicKstate<_SiteStateTraitT>::DynamicKstate(const OtherRangeT& r, CtrFromRange) :
    _v(init_vector_from_range(r)) {
}

// ***********************************************************************

template <typename _SiteStateTraitT>
typename DynamicKstate<_SiteStateTraitT>::ConstRangeT
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
// ## TraitsFor kstate_impl::DynamicKstate                              ##
// #######################################################################

namespace kstate_trait {

template<typename _SiteStateTraitT>
struct TraitKstate<kstate::DynamicKstate<_SiteStateTraitT>> {
    // the is_kstate_trait flag:
    static constexpr bool is_kstate_trait = true;
    // helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using KstateT = kstate::DynamicKstate<_SiteStateTraitT>;
    using ConstRangeT = typename kstate::DynamicKstateTypes<SiteStateTraitT>::ConstRangeT;
    using ConstAnyRangeT = typename kstate::DynamicKstateTypes<SiteStateTraitT>::ConstAnyRangeT;
    // function being the public API:
    template <typename OtherRangeT>
    static KstateT from_range(const OtherRangeT& range) {
        return KstateT(range, kstate::CtrFromRange{});
    }
    template <typename OtherRangeT>
    static std::shared_ptr<KstateT> shared_from_range(const OtherRangeT& range) {
        return std::make_shared<KstateT>(range, kstate::CtrFromRange{});
    }
    static ConstRangeT to_range(const KstateT& kstate) {
        return kstate.to_range();
    }
    static ConstAnyRangeT to_any_range(const KstateT& kstate) {
        return kstate.to_any_range();
    }
    static size_t n_sites(const KstateT& kstate){
        return kstate.is_prolific();
    }
    static size_t n_least_replication_shift(const KstateT& kstate) {
        return kstate.n_least_replication_shift();
    }
    static double norm_factor(const KstateT& kstate) {
        return kstate.norm_factor();
    }
    static bool is_prolific(const KstateT& kstate, int n_k) {
        return kstate.is_prolific(n_k);
    }
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
