#pragma once

#include <kstate_impl/kstate_abstract.hpp>
#include <kstate_impl/kstate_constructor_flavor_tag.hpp>

#include <kstate_trait/trait_site_state.hpp>
#include <kstate_trait/trait_kstate.hpp>

#include <kstate_op_integral/integral_bits_buffer.hpp>
#include <kstate_op_integral/integral_bits_view.hpp>
#include <kstate_op_integral/op_integral_bits_compare.hpp>
#include <kstate_op_integral/op_integral_bits_least_replication_shift.hpp>
#include <kstate_op_integral/op_integral_bits_unique_shift.hpp>

#include <kstate_op_integral/op_integral_bits_raw.hpp>
#include <kstate_op_range/op_range_raw_adaptors.hpp> //TODO REMOVE
#include <kstate_op_range/op_range_unique_shift.hpp> //TODO REMOVE
#include <kstate_op_range/op_range_least_replication_shift.hpp> //TODO REMOVE

#include <boost/range.hpp>
#include <boost/range/any_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <cstdint> // for types like: uint64_t
#include <iterator>
#include <type_traits>
#include <memory>
#include <cassert>

// #######################################################################
// ## init_vector_from_range (a helper function)                        ##
// #######################################################################

namespace kstate_impl::helpers {

constexpr bool bool_from_unsiged(unsigned u) {
    return (u ? true : false);
}//TODO remove

constexpr unsigned bool_to_unsigned(bool b) {
    return (b ? 1u : 0u);
}//TODO remove

template<typename SiteStateTraitT, typename IntegralT>
auto integral_bits_to_two_level_site_state_range(const kstate_op_integral::IntegralBitsDynamicBuffer<IntegralT>& integral_bits) noexcept {
    static_assert(kstate_trait::IsTraitSiteState<SiteStateTraitT>::value);
    static_assert(SiteStateTraitT::is_site_state_trait);
    static_assert(SiteStateTraitT::site_basis_dim() == 2);
    const IntegralT integral_number = integral_bits.get_number();
    const unsigned char n_all_bits = integral_bits.get_n_all_bits();
    return kstate_op_integral::raw::integral_to_chunk_numbers_range<IntegralT, unsigned>(integral_number, n_all_bits, 1u)
            | boost::adaptors::transformed(SiteStateTraitT::from_index);
}

template<typename SiteStateTraitT, typename IntegralT>
using IntegralNumberToTwoLevelSiteStateRangeResult = decltype(integral_bits_to_two_level_site_state_range<SiteStateTraitT, IntegralT>(std::declval<kstate_op_integral::IntegralBitsDynamicBuffer<IntegralT>>()));

template<typename SiteStateTraitT, typename IntegralT, typename RangeT>
kstate_op_integral::IntegralBitsDynamicBuffer<IntegralT>
integral_bits_from_two_level_site_state_range(RangeT r) noexcept {
    static_assert(kstate_trait::IsTraitSiteState<SiteStateTraitT>::value);
    static_assert(SiteStateTraitT::is_site_state_trait);
    static_assert(SiteStateTraitT::site_basis_dim() == 2);
    const auto chunk_numbers_range = r
            | boost::adaptors::transformed(SiteStateTraitT::get_index);
    const IntegralT integral_number = kstate_op_integral::raw::integral_from_chunk_numbers_range<IntegralT, unsigned>(chunk_numbers_range, 1u);
    const unsigned char n_all_bits = static_cast<unsigned char>(boost::size(r));
    return kstate_op_integral::IntegralBitsDynamicBuffer<IntegralT>{integral_number, n_all_bits};
}

} // end of namespace kstate_impl

// #######################################################################
// ## DynamicTwoLevelIntegralKstate                                     ##
// #######################################################################

namespace kstate_impl {

template<typename SiteStateTraitT>
using DynamicTwoLevelIntegral64KstateRange = helpers::IntegralNumberToTwoLevelSiteStateRangeResult<SiteStateTraitT, uint64_t>;

template <typename _SiteStateTraitT>
class DynamicTwoLevelIntegral64Kstate final : public Kstate<_SiteStateTraitT, DynamicTwoLevelIntegral64KstateRange<_SiteStateTraitT>> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(_SiteStateTraitT::site_basis_dim() == 2);
public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t>;
    using RangeT = DynamicTwoLevelIntegral64KstateRange<SiteStateTraitT>;
    using ConstRangeT = DynamicTwoLevelIntegral64KstateRange<SiteStateTraitT>;

public:
    DynamicTwoLevelIntegral64Kstate(BufferT, CtrFromBuffer);
    template <typename OtherRangeT>
    DynamicTwoLevelIntegral64Kstate(OtherRangeT, CtrFromRange);

public:
    ConstRangeT to_range() const noexcept override;
    size_t n_sites() const noexcept override;

    //protected://TODO think
    const BufferT _integral_bits;
};

// ***********************************************************************

template <typename _SiteStateTraitT>
DynamicTwoLevelIntegral64Kstate<_SiteStateTraitT>::DynamicTwoLevelIntegral64Kstate(
        DynamicTwoLevelIntegral64Kstate<_SiteStateTraitT>::BufferT integral_bits,
        CtrFromBuffer) :
    _integral_bits(integral_bits) {
}


template <typename _SiteStateTraitT>
template <typename OtherRangeT>
DynamicTwoLevelIntegral64Kstate<_SiteStateTraitT>::DynamicTwoLevelIntegral64Kstate(OtherRangeT r, CtrFromRange) :
    _integral_bits(helpers::integral_bits_from_two_level_site_state_range<SiteStateTraitT, uint64_t, OtherRangeT>(r)) {
}

// ***********************************************************************

template <typename _SiteStateTraitT>
typename DynamicTwoLevelIntegral64Kstate<_SiteStateTraitT>::ConstRangeT
DynamicTwoLevelIntegral64Kstate<_SiteStateTraitT>::to_range() const noexcept {
    return helpers::integral_bits_to_two_level_site_state_range<SiteStateTraitT, uint64_t>(_integral_bits);
}


template <typename _SiteStateTraitT>
size_t
DynamicTwoLevelIntegral64Kstate<_SiteStateTraitT>::n_sites() const noexcept {
    return _integral_bits.get_n_all_bits();
}

} // end of namespace kstate_impl

// #######################################################################
// ## TraitsFor kstate_impl::DynamicTwoLevelIntegral64Kstate            ##
// #######################################################################

namespace kstate_trait {

template<typename _SiteStateTraitT>
struct TraitKstate<kstate_impl::DynamicTwoLevelIntegral64Kstate<_SiteStateTraitT>> {
    // the is_kstate_trait flag:
    static constexpr bool is_kstate_trait = true;
    // helper types:
    using SiteStateTraitT = _SiteStateTraitT;
    using KstateT = kstate_impl::DynamicTwoLevelIntegral64Kstate<_SiteStateTraitT>;
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
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        const typename KstateT::BufferT buffer{v.get_number(), v.get_n_all_bits()};
        return KstateT(buffer, kstate_impl::CtrFromBuffer{});
    }
    template<typename ViewT>
    static std::shared_ptr<KstateT> shared_from_view(const ViewT& v) {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        const typename KstateT::BufferT buffer{v.get_number(), v.get_n_all_bits()};
        return std::make_shared<KstateT>(buffer, kstate_impl::CtrFromBuffer{});
    }
    static kstate_op_integral::IntegralBitsView<typename KstateT::BufferT>
    to_view(const KstateT& kstate) noexcept {
        return kstate_op_integral::IntegralBitsView<typename KstateT::BufferT>{kstate._integral_bits};
    }
    template<typename ViewT>
    static typename SiteStateTraitT::SiteStateT view_n_th_site_state(const ViewT& v, unsigned idx) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        const bool site_bit = kstate_op_integral::raw::extract_bit(v.get_number() , idx);
        const unsigned site_unsigned = kstate_impl::helpers::bool_to_unsigned(site_bit);
        const typename SiteStateTraitT::SiteStateT site_state = SiteStateTraitT::from_index(site_unsigned);
        return site_state;
    }
    template<typename View1T, typename View2T>
    static bool view_compare_less(const View1T& v1, const View2T& v2) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<View1T>::value);
        static_assert(kstate_op_integral::IsIntegralBits<View2T>::value);
        return kstate_op_integral::compare_less(v1, v2);
    }
    template<typename View1T, typename View2T>
    static bool view_compare_equality(const View1T& v1, const View2T& v2) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<View1T>::value);
        static_assert(kstate_op_integral::IsIntegralBits<View2T>::value);
        return kstate_op_integral::compare_equality(v1, v2);
    }
    template<typename ViewT>
    static kstate_op_integral::IntegralBitsRefinedView<ViewT, SiteStateTraitT>
    refined_view(const ViewT& v, const kstate_view_amend_spec::RefinedHolder<typename SiteStateTraitT::SiteStateT>& h) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        return kstate_op_integral::IntegralBitsRefinedView<ViewT, SiteStateTraitT>{v, h};
    }
    template<typename ViewT>
    static kstate_op_integral::IntegralBitsRotatedView<ViewT>
    rotated_view(const ViewT& v, const kstate_view_amend_spec::RotateHolder& h) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        return kstate_op_integral::IntegralBitsRotatedView<ViewT>{v, h};
    }
    template<typename ViewT>
    static size_t view_n_least_replication_shift(const ViewT& v) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        return kstate_op_integral::n_least_replication_shift(v, 1);
    }
    template<typename ViewT>
    static size_t view_n_unique_shift(const ViewT& v) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        return kstate_op_integral::n_unique_shift(v, 1);
    }
    template<typename ViewT>
    static size_t is_prolific(const ViewT& v, unsigned n_k) noexcept {
        static_assert(kstate_op_integral::IsIntegralBits<ViewT>::value);
        return kstate_op_integral::is_prolific(v, 1, n_k);
    }
};

} // end of namespace kstate_trait
