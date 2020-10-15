#pragma once

#include <kstate_impl/kstate_abstract.hpp>

#include <kstate_trait/trait_site_state.hpp>
#include <kstate_trait/trait_kstate.hpp>

#include <boost/range.hpp>
#include <boost/range/any_range.hpp>

#include <iterator>
#include <type_traits>
#include <vector>
#include <array>
#include <memory>
#include <cassert>

namespace kstate_impl {

// Helper tag classes:
struct CtrFromRange {};
struct CtrFromBuffer {};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static CtrFromRange ctr_from_range{};
static CtrFromBuffer ctr_from_buffer{};
#pragma GCC diagnostic pop

} // end of namespace kstate_impl

// #######################################################################
// ## init_vector_from_range (a helper function)                        ##
// #######################################################################

namespace kstate_impl {

template <typename Range>
std::vector<typename boost::range_value<Range>::type>
init_vector_from_range(
        const Range& rng) {
    using Element = typename boost::range_value<Range>::type;
    using Vector = std::vector<Element>;
    return Vector(std::begin(rng), std::end(rng));
}

} // end of namespace kstate_impl

// #######################################################################
// ## init_array_from_range (a helper function)                         ##
// #######################################################################

// Solution coppied from:
// https://stackoverflow.com/questions/10929202/initialize-stdarray-with-a-range-pair-of-iterators

namespace kstate_impl {

template <std::size_t... Indices>
struct indices {
    using next = indices<Indices..., sizeof...(Indices)>;
};

template <std::size_t N>
struct build_indices {
    using type = typename build_indices<N-1>::type::next;
};

template <>
struct build_indices<0> {
    using type = indices<>;
};

template <std::size_t N>
using BuildIndices = typename build_indices<N>::type;

template <typename Iterator>
using ValueType = typename std::iterator_traits<Iterator>::value_type;

// internal overload with indices tag
template <std::size_t... I, typename RandomAccessIterator,
          typename Array = std::array<ValueType<RandomAccessIterator>, sizeof...(I)>>
Array make_array_impl(RandomAccessIterator first, indices<I...>) {
    return Array { { first[I]... } };
}

// externally visible interface
template <std::size_t N, typename RandomAccessIterator>
std::array<ValueType<RandomAccessIterator>, N>
make_array(RandomAccessIterator first, RandomAccessIterator last) {
    // last is not relevant if we're assuming the size is N
    // I'll assert it is correct anyway
    assert(last - first == N);
    return make_array_impl(first, BuildIndices<N> {});
}

template <std::size_t N, typename Range>
std::array<typename boost::range_value<Range>::type, N>
init_array_from_range(
        const Range& rng) {
    return make_array<N>(std::begin(rng), std::end(rng));
}

} // end of namespace kstate_impl

// #######################################################################
// ## DynamicKstate                                                     ##
// #######################################################################

/*
 * DynamicKstate<SiteStateT> is a concrete subclass of
 * Kstate<SiteStateT, boost::random_access_traversal_tag>
 * that uses std::vector as an internal buffer.
 */

namespace kstate_impl {

template <typename _SiteStateTraitT>
struct DynamicKstateTypes {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = typename std::vector<SiteStateT>;
    using IteratorT = typename BufferT::iterator;
    using ConstIteratorT = typename BufferT::const_iterator;
    using RangeT = typename boost::iterator_range<IteratorT>;
    using ConstRangeT = typename boost::iterator_range<ConstIteratorT>;
    using AnyRangeT = typename boost::any_range<SiteStateT, boost::random_access_traversal_tag>;
    using ConstAnyRangeT = typename boost::any_range<const SiteStateT, boost::random_access_traversal_tag>;
};

template <typename _SiteStateTraitT>
class DynamicKstate final : public SpeedyKstate<_SiteStateTraitT, typename DynamicKstateTypes<_SiteStateTraitT>::ConstRangeT> {
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
    ConstRangeT to_range() const noexcept override;
    size_t n_sites() const noexcept override;

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
DynamicKstate<_SiteStateTraitT>::to_range() const noexcept {
    return _v;
}

template <typename _SiteStateTraitT>
size_t
DynamicKstate<_SiteStateTraitT>::n_sites() const noexcept {
    return _v.size();
}

}  // namespace kstate_impl

// #######################################################################
// ## TraitsFor kstate_impl::DynamicKstate                              ##
// #######################################################################

namespace kstate_trait {

template<typename _SiteStateTraitT>
struct TraitKstate<kstate_impl::DynamicKstate<_SiteStateTraitT>> {
    // the is_kstate_trait flag:
    static constexpr bool is_kstate_trait = true;
    // helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using KstateT = kstate_impl::DynamicKstate<_SiteStateTraitT>;
    using ConstRangeT = typename kstate_impl::DynamicKstateTypes<SiteStateTraitT>::ConstRangeT;
    using ConstAnyRangeT = typename kstate_impl::DynamicKstateTypes<SiteStateTraitT>::ConstAnyRangeT;
    // function being the public API:
    template <typename OtherRangeT>
    static KstateT from_range(const OtherRangeT& range) {
        return KstateT(range, kstate_impl::CtrFromRange{});
    }
    template <typename OtherRangeT>
    static std::shared_ptr<KstateT> shared_from_range(const OtherRangeT& range) {
        return std::make_shared<KstateT>(range, kstate_impl::CtrFromRange{});
    }
    static ConstRangeT to_range(const KstateT& kstate) noexcept {
        return kstate.to_range();
    }
    static ConstAnyRangeT to_any_range(const KstateT& kstate) noexcept {
        return kstate.to_any_range();
    }
    static size_t n_sites(const KstateT& kstate) noexcept {
        return kstate.is_prolific();
    }
    static size_t n_least_replication_shift(const KstateT& kstate) noexcept {
        return kstate.n_least_replication_shift();
    }
    static double norm_factor(const KstateT& kstate) noexcept {
        return kstate.norm_factor();
    }
    static bool is_prolific(const KstateT& kstate, int n_k) noexcept {
        return kstate.is_prolific(n_k);
    }
};

} // end of namespace kstate_trait

// #######################################################################
// ## StaticKstate                                                      ##
// #######################################################################

/*
 * StaticKstate<SiteStateT, N> is a concrete subclass of
 * Kstate<SiteStateT, boost::random_access_traversal_tag>
 * that uses std::array as an internal buffer.
 */

namespace kstate_impl {

template <typename _SiteStateTraitT, std::size_t _N>
struct StaticKstateTypes {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    constexpr static size_t N = _N;
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = typename std::array<SiteStateT, N>;
    using IteratorT = typename BufferT::iterator;
    using ConstIteratorT = typename BufferT::const_iterator;
    using RangeT = typename boost::iterator_range<IteratorT>;
    using ConstRangeT = typename boost::iterator_range<ConstIteratorT>;
    using AnyRangeT = typename boost::any_range<SiteStateT, boost::random_access_traversal_tag>;
    using ConstAnyRangeT = typename boost::any_range<const SiteStateT, boost::random_access_traversal_tag>;
};

template <typename _SiteStateTraitT, std::size_t _N>
class StaticKstate final : public SpeedyKstate<_SiteStateTraitT, typename StaticKstateTypes<_SiteStateTraitT, _N>::ConstRangeT> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
public:
    using SiteStateTraitT = _SiteStateTraitT;
    constexpr static std::size_t N = _N;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = typename StaticKstateTypes<SiteStateTraitT, N>::BufferT;
    using IteratorT = typename StaticKstateTypes<SiteStateTraitT, N>::IteratorT;
    using ConstIteratorT = typename StaticKstateTypes<SiteStateTraitT, N>::ConstIteratorT;
    using RangeT = typename StaticKstateTypes<SiteStateTraitT, N>::RangeT;
    using ConstRangeT = typename StaticKstateTypes<SiteStateTraitT, N>::ConstRangeT;
    using AnyRangeT = typename StaticKstateTypes<SiteStateTraitT, N>::AnyRangeT;
    using ConstAnyRangeT = typename StaticKstateTypes<SiteStateTraitT, N>::ConstAnyRangeT;

public:
    StaticKstate(BufferT&&, CtrFromBuffer);
    template <typename OtherRangeT>
    StaticKstate(const OtherRangeT&, CtrFromRange);

public:
    ConstRangeT to_range() const noexcept override;
    size_t n_sites() const noexcept override;

protected:
    const BufferT _a;
};

// ***********************************************************************

template <typename _SiteStateTraitT, std::size_t _N>
StaticKstate<_SiteStateTraitT, _N>::StaticKstate(
        StaticKstate<_SiteStateTraitT, _N>::BufferT&& a,
        CtrFromBuffer) :
   _a(std::move(a)) {
}

template <typename _SiteStateTraitT, std::size_t _N>
template <typename OtherRangeT>
StaticKstate<_SiteStateTraitT, _N>::StaticKstate(const OtherRangeT& r, CtrFromRange) :
    _a(init_array_from_range<N>(r)) {
}

// ***********************************************************************

template <typename _SiteStateTraitT, std::size_t _N>
typename StaticKstate<_SiteStateTraitT, _N>::ConstRangeT
StaticKstate<_SiteStateTraitT, _N>::to_range() const noexcept {
    return _a;
}

template <typename _SiteStateTraitT, std::size_t _N>
size_t
StaticKstate<_SiteStateTraitT, _N>::n_sites() const noexcept {
    return N;
}

}  // namespace kstate_impl

// #######################################################################
// ## TraitsFor kstate_impl::StaticKstate                               ##
// #######################################################################

namespace kstate_trait {

template<typename _SiteStateTraitT, std::size_t _N>
struct TraitKstate<kstate_impl::StaticKstate<_SiteStateTraitT, _N>> {
    // the is_kstate_trait flag:
    static constexpr bool is_kstate_trait = true;
    // helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using KstateT = kstate_impl::StaticKstate<_SiteStateTraitT, _N>;
    using ConstRangeT = typename kstate_impl::StaticKstateTypes<SiteStateTraitT, _N>::ConstRangeT;
    using ConstAnyRangeT = typename kstate_impl::StaticKstateTypes<SiteStateTraitT, _N>::ConstAnyRangeT;
    // function being the public API:
    template <typename OtherRangeT>
    static KstateT from_range(const OtherRangeT& range) {
        return KstateT(range, kstate_impl::CtrFromRange{});
    }
    template <typename OtherRangeT>
    static std::shared_ptr<KstateT> shared_from_range(const OtherRangeT& range) {
        return std::make_shared<KstateT>(range, kstate_impl::CtrFromRange{});
    }
    static ConstRangeT to_range(const KstateT& kstate) noexcept {
        return kstate.to_range();
    }
    static ConstAnyRangeT to_any_range(const KstateT& kstate) noexcept {
        return kstate.to_any_range();
    }
    static size_t n_sites(const KstateT& kstate) noexcept {
        return kstate.is_prolific();
    }
    static size_t n_least_replication_shift(const KstateT& kstate) noexcept {
        return kstate.n_least_replication_shift();
    }
    static double norm_factor(const KstateT& kstate) noexcept {
        return kstate.norm_factor();
    }
    static bool is_prolific(const KstateT& kstate, int n_k) noexcept {
        return kstate.is_prolific(n_k);
    }
};

} // end of namespace kstate_trait

