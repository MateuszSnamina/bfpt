#pragma once

#include<kstate_op_integral/integral_bits.hpp>
#include<kstate_op_integral/op_integral_bits.hpp>

#include <optional>
#include <cassert>

// #######################################################################
// ## compare_XXX                                                       ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
bool compare_less(IntegralBitsT integral_bits_1, IntegralBitsT integral_bits_2) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    assert(integral_bits_1.get_n_all_bits() == integral_bits_2.get_n_all_bits());
    assert(integral_bits_1.get_n_all_bits() < 8 * sizeof(typename IntegralBitsT::BufferT));
    return integral_bits_1.get_buffer() < integral_bits_2.get_buffer();
}

template <typename IntegralBitsT>
bool compare_equality(IntegralBitsT n1, IntegralBitsT n2) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    assert(n1.get_n_all_bits() == n2.get_n_all_bits());
    assert(n1.get_n_all_bits() < 8 * sizeof(typename IntegralBitsT::BufferT));
    return n1.get_buffer() == n2.get_buffer();
}

template <typename IntegralBitsT>
std::optional<size_t> compare_translational_equality(IntegralBitsT integral_bits_1, IntegralBitsT integral_bits_2) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    assert(integral_bits_1.get_n_all_bits() == integral_bits_2.get_n_all_bits());
    const auto n_all_bits = integral_bits_1.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::BufferT));
    const auto n1_buffer = integral_bits_1.get_buffer();
    auto n2_buffer = integral_bits_2.get_buffer();
    for (size_t i = 0; i < integral_bits_1.get_n_all_bits(); i++) {
        if (n1_buffer == n2_buffer) {
            return i;
        }
        n2_buffer = ::kstate_op_integral::rotate(n2_buffer, n_all_bits, 1);
    }
    return std::nullopt;
}

}  // namespace kstate_op_integral

