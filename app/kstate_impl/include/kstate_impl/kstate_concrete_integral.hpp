#pragma once

#include <kstate_impl/kstate_abstract.hpp>
#include <kstate_impl/kstate_constructor_flavor_tag.hpp>

#include <kstate_trait/trait_site_state.hpp>
#include <kstate_trait/trait_kstate.hpp>

#include <kstate_op_integral/integral_bits_view.hpp>
#include <kstate_op_integral/integral_bits_buffer.hpp>
#include <kstate_op_integral/op_integral_bits_compare.hpp>
#include <kstate_op_integral/op_integral_bits_least_replication_shift.hpp>
#include <kstate_op_integral/op_integral_bits_unique_shift.hpp>
#include <kstate_op_integral/op_integral_bits_raw.hpp>

#include <cstdint>  // for types like: uint64_t
#include <type_traits>
#include <memory>
#include <cassert>

// #######################################################################
// ## init_vector_from_range (a helper function)                        ##
// #######################################################################

namespace kstate_impl::helpers {

template <typename SiteStateTraitT, typename IntegralT>
auto integral_bits_to_site_state_range(
        const kstate_op_integral::IntegralBitsDynamicBuffer<IntegralT>& integral_bits,
        unsigned char n_bits_per_site) noexcept {
    static_assert(kstate_trait::IsTraitSiteState<SiteStateTraitT>::value);
    static_assert(SiteStateTraitT::is_site_state_trait);
    assert(n_bits_per_site >= 1u);
    assert(SiteStateTraitT::site_basis_dim() <= (1u << n_bits_per_site));
    const IntegralT integral_number = integral_bits.get_number();
    const unsigned char n_all_bits = integral_bits.get_n_all_bits();
    return kstate_op_integral::raw::integral_to_chunk_numbers_range<IntegralT, unsigned>(integral_number, n_all_bits, n_bits_per_site) | boost::adaptors::transformed(SiteStateTraitT::from_index);
}

template <typename SiteStateTraitT, typename IntegralT>
using IntegralNumberToSiteStateRangeResult = decltype(
integral_bits_to_site_state_range<SiteStateTraitT, IntegralT>(
std::declval<kstate_op_integral::IntegralBitsDynamicBuffer<IntegralT>>(),
std::declval<unsigned char>()));

template <typename SiteStateTraitT, typename IntegralT, typename RangeT>
kstate_op_integral::IntegralBitsDynamicBuffer<IntegralT>
integral_bits_from_site_state_range(
        RangeT r,
        unsigned char n_bits_per_site) noexcept {
    static_assert(kstate_trait::IsTraitSiteState<SiteStateTraitT>::value);
    static_assert(SiteStateTraitT::is_site_state_trait);
    assert(n_bits_per_site >= 1u);
    assert(SiteStateTraitT::site_basis_dim() <= (1u << n_bits_per_site));
    const auto chunk_numbers_range = r | boost::adaptors::transformed(SiteStateTraitT::get_index);
    const IntegralT integral_number = kstate_op_integral::raw::integral_from_chunk_numbers_range<IntegralT, unsigned>(chunk_numbers_range, n_bits_per_site);
    const unsigned char n_all_bits = static_cast<unsigned char>(boost::size(r)) * n_bits_per_site;
    return kstate_op_integral::IntegralBitsDynamicBuffer<IntegralT>{integral_number, n_all_bits};
}

}  // namespace kstate_impl::helpers

// #######################################################################
// ## DynamicIntegral64Kstate                                           ##
// #######################################################################

namespace kstate_impl {

template <typename SiteStateTraitT>
using DynamicIntegral64KstateRange = helpers::IntegralNumberToSiteStateRangeResult<SiteStateTraitT, uint64_t>;

template <typename _SiteStateTraitT, unsigned char _n_bits_per_site>
class DynamicIntegral64Kstate final : public Kstate<_SiteStateTraitT, DynamicIntegral64KstateRange<_SiteStateTraitT>> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(_n_bits_per_site >= 1u);
    static_assert(_SiteStateTraitT::site_basis_dim() <= (1u << _n_bits_per_site));

public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t>;
    using RangeT = DynamicIntegral64KstateRange<SiteStateTraitT>;
    using ConstRangeT = DynamicIntegral64KstateRange<SiteStateTraitT>;
    constexpr static unsigned char n_bits_per_site = _n_bits_per_site;

public:
    DynamicIntegral64Kstate(BufferT, CtrFromBuffer);
    template <typename OtherRangeT>
    DynamicIntegral64Kstate(OtherRangeT, CtrFromRange);

public:
    ConstRangeT to_range() const noexcept override;
    size_t n_sites() const noexcept override;
    const BufferT _integral_bits;
};

// ***********************************************************************

template <typename _SiteStateTraitT, unsigned char _n_bits_per_site>
DynamicIntegral64Kstate<_SiteStateTraitT, _n_bits_per_site>::DynamicIntegral64Kstate(
        DynamicIntegral64Kstate<_SiteStateTraitT, _n_bits_per_site>::BufferT integral_bits,
        CtrFromBuffer) : _integral_bits(integral_bits) {
    assert(_integral_bits.get_n_all_bits() % n_bits_per_site == 0);
}

template <typename _SiteStateTraitT, unsigned char _n_bits_per_site>
template <typename OtherRangeT>
DynamicIntegral64Kstate<_SiteStateTraitT, _n_bits_per_site>::DynamicIntegral64Kstate(OtherRangeT r, CtrFromRange)
    : _integral_bits(helpers::integral_bits_from_site_state_range<SiteStateTraitT, uint64_t, OtherRangeT>(r, n_bits_per_site)) {
}

// ***********************************************************************

template <typename _SiteStateTraitT, unsigned char _n_bits_per_site>
typename DynamicIntegral64Kstate<_SiteStateTraitT, _n_bits_per_site>::ConstRangeT
DynamicIntegral64Kstate<_SiteStateTraitT, _n_bits_per_site>::to_range() const noexcept {
    return helpers::integral_bits_to_site_state_range<SiteStateTraitT, uint64_t>(_integral_bits, n_bits_per_site);
}

template <typename _SiteStateTraitT, unsigned char _n_bits_per_site>
size_t
DynamicIntegral64Kstate<_SiteStateTraitT, _n_bits_per_site>::n_sites() const noexcept {
    return _integral_bits.get_n_all_bits() / n_bits_per_site;
}

}  // end of namespace kstate_impl

// #######################################################################
// ## TraitsFor kstate_impl::DynamicIntegral64Kstate                    ##
// #######################################################################

namespace kstate_trait {

template <typename _SiteStateTraitT, unsigned char _n_bits_per_site>
struct TraitKstate<kstate_impl::DynamicIntegral64Kstate<_SiteStateTraitT, _n_bits_per_site>> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(_n_bits_per_site >= 1u);
    static_assert(_SiteStateTraitT::site_basis_dim() <= (1u << _n_bits_per_site));
    // the is_kstate_trait flag:
    static constexpr bool is_kstate_trait = true;
    // helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using KstateT = kstate_impl::DynamicIntegral64Kstate<_SiteStateTraitT, _n_bits_per_site>;
    using ConstRangeT = typename KstateT::ConstRangeT;
    // function being the public API:
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
        [[maybe_unused]] unsigned n_sites = boost::size(range);
        assert(n_sites * KstateT::n_bits_per_site <= 8 * sizeof(typename KstateT::BufferT::IntegralT));
        return KstateT(range, kstate_impl::CtrFromRange{});
    }
    template <typename OtherRangeT>
    static std::shared_ptr<KstateT> shared_from_range(const OtherRangeT& range) {
        [[maybe_unused]] unsigned n_sites = boost::size(range);
        assert(n_sites * KstateT::n_bits_per_site <= 8 * sizeof(typename KstateT::BufferT::IntegralT));
        return std::make_shared<KstateT>(range, kstate_impl::CtrFromRange{});
    }
    static ConstRangeT to_range(const KstateT& kstate) noexcept {
        return kstate.to_range();
    }
    // -----
    template <typename ViewT>
    static KstateT from_view(const ViewT& v) {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        static_assert(std::is_same_v<typename ViewT::IntegralT, typename KstateT::BufferT::IntegralT>);
        const typename KstateT::BufferT buffer{v.get_number(), v.get_n_all_bits()};
        return KstateT(buffer, kstate_impl::CtrFromBuffer{});
    }
    template <typename ViewT>
    static std::shared_ptr<KstateT> shared_from_view(const ViewT& v) {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        static_assert(std::is_same_v<typename ViewT::IntegralT, typename KstateT::BufferT::IntegralT>);
        const typename KstateT::BufferT buffer{v.get_number(), v.get_n_all_bits()};
        return std::make_shared<KstateT>(buffer, kstate_impl::CtrFromBuffer{});
    }
    static kstate_op_integral::IntegralBitsView<typename KstateT::BufferT>
    to_view(const KstateT& kstate) noexcept {
        return kstate_op_integral::IntegralBitsView<typename KstateT::BufferT>{kstate._integral_bits};
    }
    template <typename ViewT>
    static typename SiteStateTraitT::SiteStateT view_n_th_site_state(const ViewT& v, unsigned idx) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        assert(v.get_n_all_bits() % KstateT::n_bits_per_site == 0);
        assert(idx < v.get_n_all_bits() / KstateT::n_bits_per_site);
        const unsigned site_unsigned = kstate_op_integral::raw::extract_chunk_number<uint64_t, unsigned>(
                    v.get_number(),
                    idx * KstateT::n_bits_per_site,
                    KstateT::n_bits_per_site);
        const typename SiteStateTraitT::SiteStateT site_state = SiteStateTraitT::from_index(site_unsigned);
        return site_state;
    }
    template <typename View1T, typename View2T>
    static bool view_compare_less(const View1T& v1, const View2T& v2) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<View1T>::value);
        static_assert(kstate_op_integral::IsIntegralBits<View2T>::value);
        static_assert(std::is_same_v<typename View1T::IntegralT, typename KstateT::BufferT::IntegralT>);
        static_assert(std::is_same_v<typename View2T::IntegralT, typename KstateT::BufferT::IntegralT>);
        return kstate_op_integral::compare_less(v1, v2);
    }
    template <typename View1T, typename View2T>
    static bool view_compare_equality(const View1T& v1, const View2T& v2) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<View1T>::value);
        static_assert(kstate_op_integral::IsIntegralBits<View2T>::value);
        static_assert(std::is_same_v<typename View1T::IntegralT, typename KstateT::BufferT::IntegralT>);
        static_assert(std::is_same_v<typename View2T::IntegralT, typename KstateT::BufferT::IntegralT>);
        return kstate_op_integral::compare_equality(v1, v2);
    }
    template <typename ViewT>
    static kstate_op_integral::IntegralBitsRefinedView<ViewT, SiteStateTraitT, KstateT::n_bits_per_site>
    refined_view(const ViewT& v, const kstate_view_amend_spec::RefinedHolder<typename SiteStateTraitT::SiteStateT>& h) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        static_assert(std::is_same_v<typename ViewT::IntegralT, typename KstateT::BufferT::IntegralT>);
        return kstate_op_integral::IntegralBitsRefinedView<ViewT, SiteStateTraitT, KstateT::n_bits_per_site>{v, h};
    }
    template <typename ViewT>
    static kstate_op_integral::IntegralBitsRotatedView<ViewT, KstateT::n_bits_per_site>
    rotated_view(const ViewT& v, const kstate_view_amend_spec::RotateHolder& h) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        static_assert(std::is_same_v<typename ViewT::IntegralT, typename KstateT::BufferT::IntegralT>);
        return kstate_op_integral::IntegralBitsRotatedView<ViewT, KstateT::n_bits_per_site>{v, h};
    }
    template <typename ViewT>
    static size_t view_n_least_replication_shift(const ViewT& v) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        static_assert(std::is_same_v<typename ViewT::IntegralT, typename KstateT::BufferT::IntegralT>);
        return kstate_op_integral::n_least_replication_shift(v, KstateT::n_bits_per_site);
    }
    template <typename ViewT>
    static size_t view_n_unique_shift(const ViewT& v) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        static_assert(std::is_same_v<typename ViewT::IntegralT, typename KstateT::BufferT::IntegralT>);
        return kstate_op_integral::n_unique_shift(v, KstateT::n_bits_per_site);
    }
    template <typename ViewT>
    static size_t is_prolific(const ViewT& v, unsigned n_k) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        static_assert(std::is_same_v<typename ViewT::IntegralT, typename KstateT::BufferT::IntegralT>);
        return kstate_op_integral::is_prolific(v, KstateT::n_bits_per_site, n_k);
    }
};

}  // end of namespace kstate_trait
