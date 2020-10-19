#pragma once

#include <kstate_impl/kstate_abstract.hpp>
#include <kstate_impl/kstate_concrete_stl_helpers.hpp>
#include <kstate_impl/kstate_constructor_flavor_tag.hpp>

#include <kstate_op_range/op_range_raw_adaptors.hpp>
#include <kstate_op_range/op_range_unique_shift.hpp>
#include <kstate_op_range/op_range_least_replication_shift.hpp>

#include <kstate_view_amend_spec/amend_spec.hpp>

#include <kstate_trait/trait_site_state.hpp>
#include <kstate_trait/trait_kstate.hpp>

#include <boost/range.hpp>
#include <iterator>
#include <type_traits>
#include <vector>
#include <array>
#include <memory>
#include <cassert>

// #######################################################################
// ## DynamicStlKstate                                                  ##
// #######################################################################

/*
 * DynamicStlKstate<SiteStateT> is a concrete subclass of
 * Kstate<SiteStateT, boost::random_access_traversal_tag>
 * that uses std::vector as an internal buffer.
 */

namespace kstate_impl {

template <typename _SiteStateTraitT>
struct DynamicStlKstateTypes {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = typename std::vector<SiteStateT>;
    using IteratorT = typename BufferT::iterator;
    using ConstIteratorT = typename BufferT::const_iterator;
    using RangeT = typename boost::iterator_range<IteratorT>;
    using ConstRangeT = typename boost::iterator_range<ConstIteratorT>;
};

template <typename _SiteStateTraitT>
class DynamicStlKstate final : public Kstate<_SiteStateTraitT, typename DynamicStlKstateTypes<_SiteStateTraitT>::ConstRangeT> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = typename DynamicStlKstateTypes<SiteStateTraitT>::BufferT;
    using IteratorT = typename DynamicStlKstateTypes<SiteStateTraitT>::IteratorT;
    using ConstIteratorT = typename DynamicStlKstateTypes<SiteStateTraitT>::ConstIteratorT;
    using RangeT = typename DynamicStlKstateTypes<SiteStateTraitT>::RangeT;
    using ConstRangeT = typename DynamicStlKstateTypes<SiteStateTraitT>::ConstRangeT;
public:
    DynamicStlKstate(BufferT&&, CtrFromBuffer);
    template <typename OtherRangeT>
    DynamicStlKstate(const OtherRangeT&, CtrFromRange);

public:
    ConstRangeT to_range() const noexcept override;
    size_t n_sites() const noexcept override;

protected:
    const BufferT _v;
};

// ***********************************************************************

template <typename _SiteStateTraitT>
DynamicStlKstate<_SiteStateTraitT>::DynamicStlKstate(
        DynamicStlKstate<_SiteStateTraitT>::BufferT&& v,
        CtrFromBuffer) :
    _v(std::move(v)) {
}

template <typename _SiteStateTraitT>
template <typename OtherRangeT>
DynamicStlKstate<_SiteStateTraitT>::DynamicStlKstate(const OtherRangeT& r, CtrFromRange) :
    _v(helpers::init_vector_from_range(r)) {
}

// ***********************************************************************

template <typename _SiteStateTraitT>
typename DynamicStlKstate<_SiteStateTraitT>::ConstRangeT
DynamicStlKstate<_SiteStateTraitT>::to_range() const noexcept {
    return _v;
}

template <typename _SiteStateTraitT>
size_t
DynamicStlKstate<_SiteStateTraitT>::n_sites() const noexcept {
    return _v.size();
}

}  // namespace kstate_impl

// #######################################################################
// ## TraitsFor kstate_impl::DynamicStlKstate                           ##
// #######################################################################

namespace kstate_trait {

template<typename _SiteStateTraitT>
struct TraitKstate<kstate_impl::DynamicStlKstate<_SiteStateTraitT>> {
    // the is_kstate_trait flag:
    static constexpr bool is_kstate_trait = true;
    // helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using KstateT = kstate_impl::DynamicStlKstate<_SiteStateTraitT>;
    using ConstRangeT = typename kstate_impl::DynamicStlKstateTypes<SiteStateTraitT>::ConstRangeT;
    // function being the public API:
    // -----
    static size_t n_sites(const KstateT& kstate) noexcept {
        return kstate.n_sites();
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
    // -----
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
    // -----
    template<typename ViewT>
    static KstateT from_view(const ViewT& v) {
        return KstateT(v, kstate_impl::CtrFromRange{});
    }
    template<typename ViewT>
    static std::shared_ptr<KstateT> shared_from_view(const ViewT& v) {
        return std::make_shared<KstateT>(v, kstate_impl::CtrFromRange{});
    }
    static ConstRangeT to_view(const KstateT& kstate) noexcept {
        return kstate.to_range();
    }
    template<typename View1T, typename View2T>
    static bool view_compare_less(const View1T& v1, const View2T& v2) noexcept {
        return kstate_op_range::compare_less(v1, v2);
    }
    template<typename View1T, typename View2T>
    static bool view_compare_equality(const View1T& v1, const View2T& v2) noexcept {
        return kstate_op_range::compare_equality(v1, v2);
    }
    template<typename ViewT>
    static auto refined_view(const ViewT& v, const kstate_view_amend_spec::RefinedHolder<typename SiteStateTraitT::SiteStateT>& h) noexcept {
        return kstate_op_range::raw::refined(v, h);
    }
    template<typename ViewT>
    static auto rotated_view(const ViewT& v, const kstate_view_amend_spec::RotateHolder& h) noexcept {
        return kstate_op_range::raw::rotated(v, h);
    }
    template<typename ViewT>
    static auto view_n_least_replication_shift(const ViewT& v) noexcept {
        return kstate_op_range::n_least_replication_shift(v);
    }
    template<typename ViewT>
    static auto view_n_unique_shift(const ViewT& v) noexcept {
        return kstate_op_range::n_unique_shift(v);
    }
};

} // end of namespace kstate_trait

// #######################################################################
// ## StaticStlKstate                                                   ##
// #######################################################################

/*
 * StaticStlKstate<SiteStateT, N> is a concrete subclass of
 * Kstate<SiteStateT, boost::random_access_traversal_tag>
 * that uses std::array as an internal buffer.
 */

namespace kstate_impl {

template <typename _SiteStateTraitT, std::size_t _N>
struct StaticStlKstateTypes {
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
};

template <typename _SiteStateTraitT, std::size_t _N>
class StaticStlKstate final : public Kstate<_SiteStateTraitT, typename StaticStlKstateTypes<_SiteStateTraitT, _N>::ConstRangeT> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
public:
    using SiteStateTraitT = _SiteStateTraitT;
    constexpr static std::size_t N = _N;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = typename StaticStlKstateTypes<SiteStateTraitT, N>::BufferT;
    using IteratorT = typename StaticStlKstateTypes<SiteStateTraitT, N>::IteratorT;
    using ConstIteratorT = typename StaticStlKstateTypes<SiteStateTraitT, N>::ConstIteratorT;
    using RangeT = typename StaticStlKstateTypes<SiteStateTraitT, N>::RangeT;
    using ConstRangeT = typename StaticStlKstateTypes<SiteStateTraitT, N>::ConstRangeT;
public:
    StaticStlKstate(BufferT&&, CtrFromBuffer);
    template <typename OtherRangeT>
    StaticStlKstate(const OtherRangeT&, CtrFromRange);

public:
    ConstRangeT to_range() const noexcept override;
    size_t n_sites() const noexcept override;

protected:
    const BufferT _a;
};

// ***********************************************************************

template <typename _SiteStateTraitT, std::size_t _N>
StaticStlKstate<_SiteStateTraitT, _N>::StaticStlKstate(
        StaticStlKstate<_SiteStateTraitT, _N>::BufferT&& a,
        CtrFromBuffer) :
   _a(std::move(a)) {
}

template <typename _SiteStateTraitT, std::size_t _N>
template <typename OtherRangeT>
StaticStlKstate<_SiteStateTraitT, _N>::StaticStlKstate(const OtherRangeT& r, CtrFromRange) :
    _a(helpers::init_array_from_range<N>(r)) {
}

// ***********************************************************************

template <typename _SiteStateTraitT, std::size_t _N>
typename StaticStlKstate<_SiteStateTraitT, _N>::ConstRangeT
StaticStlKstate<_SiteStateTraitT, _N>::to_range() const noexcept {
    return _a;
}

template <typename _SiteStateTraitT, std::size_t _N>
size_t
StaticStlKstate<_SiteStateTraitT, _N>::n_sites() const noexcept {
    return N;
}

}  // namespace kstate_impl

// #######################################################################
// ## TraitsFor kstate_impl::StaticStlKstate                            ##
// #######################################################################


//TODO copy solution from DynamicState!

namespace kstate_trait {

template<typename _SiteStateTraitT, std::size_t _N>
struct TraitKstate<kstate_impl::StaticStlKstate<_SiteStateTraitT, _N>> {
    // the is_kstate_trait flag:
    static constexpr bool is_kstate_trait = true;
    // helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using KstateT = kstate_impl::StaticStlKstate<_SiteStateTraitT, _N>;
    using ConstRangeT = typename kstate_impl::StaticStlKstateTypes<SiteStateTraitT, _N>::ConstRangeT;
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
