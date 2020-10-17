#pragma once

#include <type_traits>
#include <cassert>

namespace kstate_op_integral {

template<typename _IntegralT>
struct IntegralBitsDynamic {
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

template<typename _IntegralT, std::size_t _N>
struct IntegralBitsStatic {
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

template<typename T>
struct IsIntegralBits : std::false_type {};

template<typename T>
struct IsIntegralBits<IntegralBitsDynamic<T>> : std::true_type {};

template<typename T, std::size_t N>
struct IsIntegralBits<IntegralBitsStatic<T, N>> : std::true_type {};

}
