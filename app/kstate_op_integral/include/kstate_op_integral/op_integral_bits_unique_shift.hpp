
#pragma once

#include<kstate_op_integral/is_integral_bits.hpp>
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
    unsigned char n_rotation;
    IntegralT buffer;
};

template <typename IntegralBitsT>
UniqueShiftReceipt<typename IntegralBitsT::IntegralT> _n_unique_shift(IntegralBitsT integral_bits, unsigned n_bits_per_site) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    assert(n_all_bits % n_bits_per_site == 0);
    const auto n_sites = n_all_bits / n_bits_per_site;
    UniqueShiftReceipt<typename IntegralBitsT::IntegralT> current_max{0, integral_bits.get_number()};
    {
        UniqueShiftReceipt<typename IntegralBitsT::IntegralT> proposed_max{0, integral_bits.get_number()};
        while (proposed_max.n_rotation + 1u < n_sites)
        {
            proposed_max.buffer = ::kstate_op_integral::raw::rotate(proposed_max.buffer, n_all_bits, n_bits_per_site);
            proposed_max.n_rotation++;
            if (current_max.buffer < proposed_max.buffer) {
                current_max = proposed_max;
            }
        }
    }
    return current_max;
}

}  // namespace kstate_op_integral

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
unsigned char n_unique_shift(IntegralBitsT integral_bits, unsigned n_bits_per_site) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    [[maybe_unused]] const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    assert(n_all_bits % n_bits_per_site == 0);
    return _n_unique_shift(integral_bits, n_bits_per_site).n_rotation;
}

}  // namespace kstate_op_integral

// #######################################################################
// ## make_unique_shift                                                 ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
typename IntegralBitsT::IntegralT make_unique_shift(IntegralBitsT integral_bits, unsigned n_bits_per_site) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    [[maybe_unused]] const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    assert(n_all_bits % n_bits_per_site == 0);
    return _n_unique_shift(integral_bits, n_bits_per_site).buffer;
}

}  // namespace kstate_op_integral
