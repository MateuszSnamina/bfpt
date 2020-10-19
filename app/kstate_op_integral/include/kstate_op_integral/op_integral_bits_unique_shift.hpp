
#pragma once

#include<kstate_op_integral/integral_bits.hpp>
#include<kstate_op_integral/op_integral_bits_raw.hpp>

#include<cassert>

// #######################################################################
// ## _n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_integral {

template<typename _IntegralT>
struct UniqueShiftReceipt {
    using IntegralT = _IntegralT;
    static_assert(std::is_arithmetic_v<_IntegralT>);
    static_assert(std::is_integral_v<_IntegralT>);
    static_assert(std::is_unsigned_v<_IntegralT>);
    unsigned char n_roration;
    IntegralT buffer;
};

template <typename IntegralBitsT>
UniqueShiftReceipt<typename IntegralBitsT::IntegralT> _n_unique_shift(IntegralBitsT integral_bits) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    auto n2_buffer = integral_bits.get_number();
    unsigned char i = 0;
    typename IntegralBitsT::IntegralT n_max_buffer = integral_bits.get_number();
    for (unsigned char _ = 1; _ < n_all_bits; _++) {
        n2_buffer = ::kstate_op_integral::raw::rotate(n2_buffer, n_all_bits, 1);
        if (n_max_buffer < n2_buffer) {
            i = _;
            n_max_buffer = n2_buffer;
        }
    }
    return {i, n_max_buffer};
}

}  // namespace kstate_op_integral

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
unsigned char n_unique_shift(IntegralBitsT integral_bits) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    [[maybe_unused]] const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    return _n_unique_shift(integral_bits).n_roration;
}

}  // namespace kstate_op_integral

// #######################################################################
// ## make_unique_shift                                                 ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
typename IntegralBitsT::IntegralT make_unique_shift(IntegralBitsT integral_bits) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    [[maybe_unused]] const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    return _n_unique_shift(integral_bits).buffer;
}

}  // namespace kstate_op_integral
