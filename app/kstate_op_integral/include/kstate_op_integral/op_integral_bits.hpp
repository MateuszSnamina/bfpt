#pragma once

#include <boost/range.hpp>
#include <boost/range/any_range.hpp>

#include <iterator>
#include <type_traits>
#include <memory>
#include <cassert>

// #######################################################################
// ## init_vector_from_range (a helper function)                        ##
// #######################################################################

namespace kstate_op_integral {

template <typename T>
bool extract_bit(T n, unsigned idx) {
    static_assert(std::is_arithmetic_v<T>);
    static_assert(std::is_integral_v<T>);
    static_assert(std::is_unsigned_v<T>);
    assert(idx < 8 * sizeof(T));
    return static_cast<bool>(n & (static_cast<T>(1u) << idx));
}

template <typename T>
T rotate(T n, unsigned n_all_bits, unsigned idx_pivot_bit) {
    static_assert(std::is_arithmetic_v<T>);
    static_assert(std::is_integral_v<T>);
    static_assert(std::is_unsigned_v<T>);
    assert(n_all_bits < 8 * sizeof(T));
    assert(idx_pivot_bit < n_all_bits);
    // low bits idx:  [0, idx_pivot_bit)
    // high bits idx: [idx_pivot_bit, n_all_bits)
    const unsigned n_low_bits = idx_pivot_bit;
    const unsigned n_high_bits = n_all_bits - idx_pivot_bit;
    const T low_bits_mask = ((static_cast<T>(1u) << n_low_bits) - 1);
    const T high_bits_mask = ((static_cast<T>(1u) << n_high_bits) - 1) << n_low_bits;
    const T low_bits = n & low_bits_mask;
    const T high_bits = n & high_bits_mask;
    const T result = (high_bits >> n_low_bits) | (low_bits << n_high_bits);
    return result;
}

} // end of namespace kstate_op_integral

