
#pragma once

#include<kstate_op_integral/integral_bits.hpp>
#include<kstate_op_integral/op_integral_bits.hpp>

#include<cassert>

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
size_t n_unique_shift(IntegralBitsT integral_bits) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::BufferT));
    auto n2_buffer = integral_bits.get_buffer();
    size_t i = 0;
    typename IntegralBitsT::BufferT n_max_buffer = integral_bits.get_buffer();
    for (size_t _ = 1; _ < n_all_bits; _++) {
        n2_buffer = ::kstate_op_integral::rotate(n2_buffer, n_all_bits, 1);
        if (n_max_buffer < n2_buffer) {
            i = _;
            n_max_buffer = n2_buffer;
        }
    }
    return i;
}

}  // namespace kstate_op_integral
/*
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
