#pragma once

#include <kstate_op_integral/is_integral_bits.hpp>

#include <type_traits>
#include <cassert>

namespace kstate_op_integral {

template <typename _IntegralT>
struct IntegralBitsDynamicBuffer {
    static_assert(std::is_arithmetic_v<_IntegralT>);
    static_assert(std::is_integral_v<_IntegralT>);
    static_assert(std::is_unsigned_v<_IntegralT>);
    using IntegralT = _IntegralT;
    const IntegralT number;
    const unsigned char n_all_bits;
    const IntegralT& get_number_ref() const noexcept {
        return number;
    }
    IntegralT get_number() const noexcept {
        return number;
    }
    unsigned char get_n_all_bits() const noexcept {
        return n_all_bits;
    };
};

template <typename _IntegralT, std::size_t _N>
struct IntegralBitsStaticBuffer {
    static_assert(std::is_arithmetic_v<_IntegralT>);
    static_assert(std::is_integral_v<_IntegralT>);
    static_assert(std::is_unsigned_v<_IntegralT>);
    static_assert(_N < 8 * sizeof(_IntegralT));
    using IntegralT = _IntegralT;
    constexpr static std::size_t N = _N;
    const IntegralT number;
    const IntegralT& get_number_ref() const noexcept {
        return number;
    }
    IntegralT get_number() const noexcept {
        return number;
    }
    unsigned char get_n_all_bits() const noexcept {
        return N;
    };
};

template <typename _IntegralT>
struct IsIntegralBits<IntegralBitsDynamicBuffer<_IntegralT>> : std::true_type {};

template <typename _IntegralT, std::size_t _N>
struct IsIntegralBits<IntegralBitsStaticBuffer<_IntegralT, _N>> : std::true_type {};

}  // namespace kstate_op_integral
