#pragma once

#include <boost/range.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include <iterator>
#include <type_traits>
#include <cassert>

// #######################################################################
// ## kstate_op_integral::raw (helper functions)                        ##
// #######################################################################

namespace kstate_op_integral::raw {

// #######################################################################
// ## extract_{bit,chunk_number}                                        ##
// #######################################################################

template <typename IntegralT>
bool extract_bit(IntegralT n, unsigned char idx) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    assert(idx < 8 * sizeof(IntegralT));
    const IntegralT mask = static_cast<IntegralT>(1u) << idx;
    return static_cast<bool>(n & mask);
}

template <typename IntegralT, typename ChunkIntegralT>
unsigned extract_chunk_number(IntegralT n, unsigned char idx_first_bit, unsigned char n_bits_in_number) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    static_assert(std::is_arithmetic_v<ChunkIntegralT>);
    static_assert(std::is_integral_v<ChunkIntegralT>);
    static_assert(std::is_unsigned_v<ChunkIntegralT>);
    static_assert(sizeof(ChunkIntegralT) <= sizeof(IntegralT));
    assert(n_bits_in_number + 1u < 8 * sizeof(IntegralT));
    assert(idx_first_bit + n_bits_in_number < 8 * sizeof(IntegralT));
    const IntegralT mask = (static_cast<IntegralT>(1u) << n_bits_in_number) - static_cast<IntegralT>(1u);
    const IntegralT result = (n & (mask << idx_first_bit)) >> idx_first_bit;
    return static_cast<ChunkIntegralT>(result);
}

// #######################################################################
// ## refine_{bit,chunk_number}                                         ##
// #######################################################################

template <typename IntegralT>
IntegralT refine_bit(IntegralT n, bool new_value, unsigned char idx_bit) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    const IntegralT mask = (static_cast<IntegralT>(1u) << idx_bit);
    const IntegralT result = (new_value ? n | mask : n & ~mask);
    return result;
}

template <typename IntegralT, typename ChunkIntegralT>
IntegralT refine_chunk_number(IntegralT n, unsigned char idx_first_bit, unsigned char n_bits_in_chunk_number, ChunkIntegralT n_chunk) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    static_assert(std::is_arithmetic_v<ChunkIntegralT>);
    static_assert(std::is_integral_v<ChunkIntegralT>);
    static_assert(std::is_unsigned_v<ChunkIntegralT>);
    static_assert(sizeof(ChunkIntegralT) <= sizeof(IntegralT));
    assert(n_bits_in_chunk_number + 1u < 8 * sizeof(IntegralT));
    assert(idx_first_bit + n_bits_in_chunk_number < 8 * sizeof(IntegralT));
    const IntegralT mask = (static_cast<IntegralT>(1u) << n_bits_in_chunk_number) - static_cast<IntegralT>(1u);
    n &= ~(mask << idx_first_bit);
    n |= (static_cast<IntegralT>(n_chunk) << idx_first_bit);
    return n;
}

// #######################################################################
// ## rotate_{bit,chunk_number}                                         ##
// #######################################################################

template <typename IntegralT>
IntegralT rotate_bit(IntegralT n, unsigned char n_all_bits, unsigned char idx_pivot_bit) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    assert(n_all_bits < 8 * sizeof(IntegralT));
    assert(idx_pivot_bit < n_all_bits);
    // low bits idx:  [0, idx_pivot_bit)
    // high bits idx: [idx_pivot_bit, n_all_bits)
    const unsigned char n_low_bits = idx_pivot_bit;
    const unsigned char n_high_bits = n_all_bits - idx_pivot_bit;
    const IntegralT low_bits_mask = ((static_cast<IntegralT>(1u) << n_low_bits) - 1);
    const IntegralT high_bits_mask = ((static_cast<IntegralT>(1u) << n_high_bits) - 1) << n_low_bits;
    const IntegralT low_bits = n & low_bits_mask;
    const IntegralT high_bits = n & high_bits_mask;
    const IntegralT result = (high_bits >> n_low_bits) | (low_bits << n_high_bits);
    return result;
}

template <typename IntegralT>
IntegralT rotate_chunk_number(IntegralT n, unsigned char n_all_bits, unsigned char idx_pivot_chunk_number, unsigned char n_bits_in_chunk_number) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    assert(n_bits_in_chunk_number != 0u);
    assert(n_all_bits < 8 * sizeof(IntegralT));
    assert(n_all_bits % n_bits_in_chunk_number == 0);
    const unsigned char idx_pivot_bit = idx_pivot_chunk_number * n_bits_in_chunk_number;
    assert(idx_pivot_bit < n_all_bits);
    return rotate_bit<IntegralT>(n, n_all_bits, idx_pivot_bit);
}

// #######################################################################
// ## integral_to_{bits,chunk_numbers}_range                            ##
// #######################################################################

template <typename IntegralT>
auto integral_to_bits_range(IntegralT n, unsigned char n_all_bits) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    assert(n_all_bits < 8 * sizeof(IntegralT));
    const auto transformation = [n](unsigned char idx_bit) -> bool {
        return extract_bit(n, idx_bit);
    };
    const auto result_range =
            boost::irange(static_cast<unsigned char>(0u), n_all_bits) |
            boost::adaptors::transformed(transformation);
    return result_range;
}

template <typename IntegralT, typename ChunkIntegralT>
auto integral_to_chunk_numbers_range(IntegralT n, unsigned char n_all_bits, unsigned char n_bits_in_chunk_number) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    static_assert(std::is_arithmetic_v<ChunkIntegralT>);
    static_assert(std::is_integral_v<ChunkIntegralT>);
    static_assert(std::is_unsigned_v<ChunkIntegralT>);
    static_assert(sizeof(ChunkIntegralT) <= sizeof(IntegralT));
    assert(n_all_bits < 8 * sizeof(IntegralT));
    assert(n_bits_in_chunk_number != 0);
    assert(n_all_bits % n_bits_in_chunk_number == 0);
    const auto transformation = [n, n_bits_in_chunk_number](unsigned char idx_bit) -> ChunkIntegralT {
        return extract_chunk_number<IntegralT, ChunkIntegralT>(n, idx_bit, n_bits_in_chunk_number);
    };
    const auto result_range =
            boost::irange(static_cast<unsigned char>(0u), n_all_bits, n_bits_in_chunk_number) |
            boost::adaptors::transformed(transformation);
    return result_range;
}

// #######################################################################
// ## integral_from_{bits,chunk_numbers}_range                          ##
// #######################################################################

template <typename IntegralT, typename RangeT>
IntegralT integral_from_bits_range(RangeT r) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    using IteratorT = typename boost::range_iterator<RangeT>::type;
    using DerefIteratorT = decltype(*std::declval<IteratorT>()); //TODO: make it better. boost::iterator_value<> ??
    using RemoveCvRefDerefIteratorT = std::remove_cv_t<std::remove_reference_t<DerefIteratorT>>; //TODO: std::remove_cvref_t<DerefIteratorIntegralT>
    static_assert(std::is_same_v<RemoveCvRefDerefIteratorT, bool>);
    [[maybe_unused]] const unsigned char n_all_bits = boost::size(r);
    assert(n_all_bits < 8 * sizeof(IntegralT));
    IntegralT mask = static_cast<IntegralT>(1u);
    IntegralT result = static_cast<IntegralT>(0u);
    for (bool b : r) {
        if (b) {
            result |= mask;
        };
        mask <<= 1;
    }
    return result;
}

template <typename IntegralT, typename ChunkIntegralT, typename RangeT>
IntegralT integral_from_chunk_numbers_range(RangeT chunk_numbers_range, unsigned char n_bits_in_chunk_number) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    static_assert(std::is_arithmetic_v<ChunkIntegralT>);
    static_assert(std::is_integral_v<ChunkIntegralT>);
    static_assert(std::is_unsigned_v<ChunkIntegralT>);
    static_assert(sizeof(ChunkIntegralT) <= sizeof(IntegralT));
    using IteratorT = typename boost::range_iterator<RangeT>::type;
    using DerefIteratorT = decltype(*std::declval<IteratorT>()); //TODO: make it better. boost::iterator_value<> ??
    using RemoveCvRefDerefIteratorT = std::remove_cv_t<std::remove_reference_t<DerefIteratorT>>; //TODO: std::remove_cvref_t<DerefIteratorIntegralT>
    static_assert(std::is_same_v<RemoveCvRefDerefIteratorT, ChunkIntegralT>);
    [[maybe_unused]] const unsigned char n_chunk_numbers = boost::size(chunk_numbers_range);
    assert(n_bits_in_chunk_number != 0);
    assert(n_chunk_numbers * n_bits_in_chunk_number < 8 * sizeof(IntegralT));
    IntegralT result = static_cast<IntegralT>(0u);
    unsigned char idx_first_bit = 0;
    for (ChunkIntegralT n_chunk : chunk_numbers_range) {
        result = refine_chunk_number<IntegralT, ChunkIntegralT>(result, idx_first_bit, n_bits_in_chunk_number, n_chunk);
        idx_first_bit += n_bits_in_chunk_number;
    }
    return result;
}

} // end of namespace kstate_op_integral::raw
