
#pragma once

#include<kstate_op_integral/integral_bits.hpp>
#include<kstate_op_integral/op_integral_bits.hpp>

#include<cassert>
/*
// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
size_t n_unique_shift(IntegralBitsT integral_bits) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
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

}  // namespace kstate_op_integral

// #######################################################################
// ## make_unique_shift                                                 ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
extension::boost::adaptors::RotatedRangeType<ForwardRange> make_unique_shift(IntegralBitsT integral_bits) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    return rng | extension::boost::adaptors::rotated(n_unique_shift(rng));
}

}  // namespace kstate_op_integral
*/
