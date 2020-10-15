#pragma once

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <iterator>
#include <type_traits>
#include <cassert>

// #######################################################################
// ## init_vector_from_range (a helper function)                        ##
// #######################################################################

namespace kstate_op_integral {

template <typename T>
bool extract_bit(T n, unsigned idx) noexcept {
    static_assert(std::is_arithmetic_v<T>);
    static_assert(std::is_integral_v<T>);
    static_assert(std::is_unsigned_v<T>);
    assert(idx < 8 * sizeof(T));
    return static_cast<bool>(n & (static_cast<T>(1u) << idx));
}

template <typename T>
T rotate(T n, unsigned n_all_bits, unsigned idx_pivot_bit) noexcept {
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

template <typename T>
T refine(T n, bool new_value, unsigned idx_bit) noexcept {
    static_assert(std::is_arithmetic_v<T>);
    static_assert(std::is_integral_v<T>);
    static_assert(std::is_unsigned_v<T>);
    const T mask = (static_cast<T>(1u) << idx_bit);
    const T result = (new_value ? n | mask : n & ~mask);
    return result;
}

template <typename T>
auto integral_to_bits_range(T n, unsigned n_all_bits) noexcept {
    static_assert(std::is_arithmetic_v<T>);
    static_assert(std::is_integral_v<T>);
    static_assert(std::is_unsigned_v<T>);
    assert(n_all_bits < 8 * sizeof(T));
    const auto r = boost::irange(0u, n_all_bits) | boost::adaptors::transformed(
                [n](unsigned idx_bit) -> bool {return extract_bit(n, idx_bit);}
    );
    return r;
}

template <typename T, typename RangeT>
T integral_from_bits_range(RangeT r) noexcept {
    static_assert(std::is_arithmetic_v<T>);
    static_assert(std::is_integral_v<T>);
    static_assert(std::is_unsigned_v<T>);
    using IteratorT = typename boost::range_iterator<RangeT>::type;
    using DerefIteratorT = decltype(*std::declval<IteratorT>());
    using RemoveCvRefDerefIteratorT = std::remove_cv_t<std::remove_reference_t<DerefIteratorT>>; //  std::remove_cvref_t<DerefIteratorT>
    static_assert(std::is_same_v<RemoveCvRefDerefIteratorT, bool>);
    [[maybe_unused]] unsigned n_all_bits = boost::size(r);
    assert(n_all_bits < 8 * sizeof(T));
    T mask = static_cast<T>(1u);
    T result = static_cast<T>(0u);
    for (bool b : r) {
        if (b) {
            result |= mask;
        };
        mask <<= 1;
    }
    return result;
}

} // end of namespace kstate_op_integral

