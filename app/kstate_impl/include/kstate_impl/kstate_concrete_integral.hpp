#pragma once

#include <kstate_impl/kstate_abstract.hpp>

#include <kstate_trait/trait_site_state.hpp>
#include <kstate_trait/trait_kstate.hpp>

#include <kstate_op_integral/integral_bits.hpp>
#include <kstate_op_integral/op_integral_bits.hpp>

#include <boost/range.hpp>
#include <boost/range/any_range.hpp>
#include <boost/range/adaptor/transformed.hpp>

//#include <cstdint> for types like: uint64_t
#include <iterator>
#include <type_traits>
#include <memory>
#include <cassert>

// #######################################################################
// ## init_vector_from_range (a helper function)                        ##
// #######################################################################

namespace kstate_impl {

bool unsigned_to_bool(unsigned u) {
    return (u ? true : false);
}

unsigned bool_to_unsigned(bool b) {
    return (b ? 1u : 0u);
}

template<typename SiteStateTraitT, typename IntegralT>
auto integral_number_to_two_level_site_state_range(const kstate_op_integral::IntegralBitsDynamic<IntegralT>& integral_bits) noexcept {
    const IntegralT integral_number = integral_bits.get_number();
    const unsigned char n_all_bits = integral_bits.get_n_all_bits();
    static_assert(kstate_trait::IsTraitSiteState<SiteStateTraitT>::value);
    static_assert(SiteStateTraitT::is_site_state_trait);
    static_assert(SiteStateTraitT::site_basis_dim() == 2);
    return kstate_op_integral::integral_to_bits_range(integral_number, n_all_bits)
            | boost::adaptors::transformed(bool_to_unsigned)
            | boost::adaptors::transformed(SiteStateTraitT::from_index);
}

template<typename SiteStateTraitT, typename IntegralT>
using IntegralNumberToTwoLevelSiteStateRangeResult = decltype(integral_number_to_two_level_site_state_range<SiteStateTraitT, IntegralT>(std::declval<kstate_op_integral::IntegralBitsDynamic<IntegralT>>()));

template<typename SiteStateTraitT, typename IntegralT, typename RangeT>
kstate_op_integral::IntegralBitsDynamic<IntegralT>
two_level_site_state_range_to_integral_number(RangeT r) noexcept {
    static_assert(kstate_trait::IsTraitSiteState<SiteStateTraitT>::value);
    static_assert(SiteStateTraitT::is_site_state_trait);
    const auto bits_range = r
            | boost::adaptors::transformed(SiteStateTraitT::get_index)
            | boost::adaptors::transformed(unsigned_to_bool);
    const IntegralT integral_number = kstate_op_integral::integral_from_bits_range<IntegralT>(bits_range);
    const unsigned char n_all_bits = static_cast<unsigned char>(boost::size(r));
    return kstate_op_integral::IntegralBitsDynamic<IntegralT>{integral_number, n_all_bits};
}

} // end of namespace kstate_impl

// #######################################################################
// ## DynamicTwoLevelIntegralKstate                                     ##
// #######################################################################

namespace kstate_impl {

template <typename _SiteStateTraitT>
class DynamicTwoLevelIntegralKstate final : public SpeedyKstate<_SiteStateTraitT, IntegralNumberToTwoLevelSiteStateRangeResult<_SiteStateTraitT, uint64_t>> {
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
public:
    using SiteStateTraitT = _SiteStateTraitT;
    using SiteStateT = typename _SiteStateTraitT::SiteStateT;
    using BufferT = kstate_op_integral::IntegralBitsDynamic<uint64_t>;
//    using IteratorT = typename StaticKstateTypes<SiteStateTraitT, N>::IteratorT;
//    using ConstIteratorT = typename StaticKstateTypes<SiteStateTraitT, N>::ConstIteratorT;
//    using RangeT = typename StaticKstateTypes<SiteStateTraitT, N>::RangeT;
//    using ConstRangeT = typename StaticKstateTypes<SiteStateTraitT, N>::ConstRangeT;
//    using AnyRangeT = typename StaticKstateTypes<SiteStateTraitT, N>::AnyRangeT;
//    using ConstAnyRangeT = typename StaticKstateTypes<SiteStateTraitT, N>::ConstAnyRangeT;

public:
//    DynamicTwoLevelIntegralKstate(BufferT, CtrFromBuffer);
//    template <typename OtherRangeT>
//    DynamicTwoLevelIntegralKstate(const OtherRangeT&, CtrFromRange);

public:
//    ConstRangeT to_range() const noexcept override;
    size_t n_sites() const noexcept override;

protected:
    const BufferT _integral_bits;
};

// ***********************************************************************

//template <typename _SiteStateTraitT>
//DynamicTwoLevelIntegralKstate<_SiteStateTraitT>::DynamicTwoLevelIntegralKstate(
//        DynamicTwoLevelIntegralKstate<_SiteStateTraitT>::BufferT integral_bits,
//        CtrFromBuffer) :
//   _integral_bits(integral_bits) {
//}

//template <typename _SiteStateTraitT>
//template <typename OtherRangeT>
//DynamicTwoLevelIntegralKstate<_SiteStateTraitT>::DynamicTwoLevelIntegralKstate(const OtherRangeT& r, CtrFromRange) :
//    _integral_bits(TODO) {
//}

// ***********************************************************************

//template <typename _SiteStateTraitT>
//typename StaticKstate<_SiteStateTraitT>::ConstRangeT
//DynamicTwoLevelIntegralKstate<_SiteStateTraitT>::to_range() const noexcept {
//    return _a;
//}

template <typename _SiteStateTraitT>
size_t
DynamicTwoLevelIntegralKstate<_SiteStateTraitT>::n_sites() const noexcept {
    return _integral_bits.get_n_all_bits();
}

} // end of namespace kstate_impl

