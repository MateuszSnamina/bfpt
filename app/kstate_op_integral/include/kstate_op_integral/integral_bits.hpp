#pragma once

#include <type_traits>
#include <cassert>

namespace kstate_op_integral {

template<typename _BufferT>
struct IntegralBitsDynamic {
    static_assert(std::is_arithmetic_v<_BufferT>);
    static_assert(std::is_integral_v<_BufferT>);
    static_assert(std::is_unsigned_v<_BufferT>);
    using BufferT = _BufferT;
    const BufferT buffer;
    const unsigned char n_all_bits;
    const BufferT& get_buffer_ref() const noexcept {
        return buffer;
    }
    BufferT get_buffer() const noexcept {
        return buffer;
    }
    unsigned char get_n_all_bits() const noexcept {
        return n_all_bits;
    };
};

template<typename _BufferT, std::size_t _N>
struct IntegralBitsStatic {
    static_assert(std::is_arithmetic_v<_BufferT>);
    static_assert(std::is_integral_v<_BufferT>);
    static_assert(std::is_unsigned_v<_BufferT>);
    static_assert(_N < 8 * sizeof(_BufferT));
    using BufferT = _BufferT;
    constexpr static std::size_t N = _N;
    const BufferT buffer;
    const BufferT& get_buffer_ref() const noexcept {
        return buffer;
    }
    BufferT get_buffer() const noexcept {
        return buffer;
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
