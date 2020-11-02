#pragma once

#include <kstate_op_integral/is_integral_bits.hpp>
#include <kstate_op_integral/op_integral_bits_raw.hpp>

#include <type_traits>
#include <optional>
#include <cassert>

// #######################################################################
// ## compare_XXX                                                       ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBits1T, typename IntegralBits2T>
bool compare_less(IntegralBits1T integral_bits_1, IntegralBits2T integral_bits_2) noexcept {
    static_assert(IsIntegralBits<IntegralBits1T>::value);
    static_assert(IsIntegralBits<IntegralBits2T>::value);
    static_assert(std::is_same_v<typename IntegralBits1T::IntegralT, typename IntegralBits2T::IntegralT>);
    assert(integral_bits_1.get_n_all_bits() == integral_bits_2.get_n_all_bits());
    assert(integral_bits_1.get_n_all_bits() < 8 * sizeof(typename IntegralBits1T::IntegralT));
    return integral_bits_1.get_number() < integral_bits_2.get_number();
}

template <typename IntegralBits1T, typename IntegralBits2T>
bool compare_equality(IntegralBits1T n1, IntegralBits2T n2) noexcept {
    static_assert(IsIntegralBits<IntegralBits1T>::value);
    static_assert(IsIntegralBits<IntegralBits2T>::value);
    static_assert(std::is_same_v<typename IntegralBits1T::IntegralT, typename IntegralBits2T::IntegralT>);
    assert(n1.get_n_all_bits() == n2.get_n_all_bits());
    assert(n1.get_n_all_bits() < 8 * sizeof(typename IntegralBits1T::IntegralT));
    return n1.get_number() == n2.get_number();
}

template <typename IntegralBits1T, typename IntegralBits2T>
std::optional<size_t> compare_translational_equality(IntegralBits1T integral_bits_1, IntegralBits2T integral_bits_2) noexcept {
    static_assert(IsIntegralBits<IntegralBits1T>::value);
    static_assert(IsIntegralBits<IntegralBits2T>::value);
    static_assert(std::is_same_v<typename IntegralBits1T::IntegralT, typename IntegralBits2T::IntegralT>);
    assert(integral_bits_1.get_n_all_bits() == integral_bits_2.get_n_all_bits());
    const auto n_all_bits = integral_bits_1.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBits1T::IntegralT));
    const auto n1_buffer = integral_bits_1.get_number();
    auto n2_buffer = integral_bits_2.get_number();
    for (size_t i = 0; i < integral_bits_1.get_n_all_bits(); i++) {
        if (n1_buffer == n2_buffer) {
            return i;
        }
        n2_buffer = ::kstate_op_integral::raw::rotate_bit(n2_buffer, n_all_bits, 1);
    }
    return std::nullopt;
}

}  // namespace kstate_op_integral
