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

template <typename IntegralT>
bool extract_bit(IntegralT n, unsigned char idx) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    assert(idx < 8 * sizeof(IntegralT));
    const IntegralT mask = static_cast<IntegralT>(1u) << idx;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
    return static_cast<bool>(n & mask);
#pragma GCC diagnostic pop
    // TODO:
    // Find why Wmaybe-uninitialized is issued
    //
    // const auto rngdr = rng | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    // const auto it = boost::range::search(rngdr, rng);
    //
    // being part of
    // template <typename ForwardRange>
    // size_t kstate_op_range::n_least_replication_shift(const ForwardRange& rng)
    //
    // when instantiated for:
    // template <typename IntegralT> auto integral_to_bits_range(IntegralT n, unsigned char n_all_bits) noexcept {
}

template <typename IntegralT>
unsigned extract_chunk_number(IntegralT n, unsigned char idx_first_bit, unsigned char n_bits_in_number) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    assert(n_bits_in_number + 1u < 8 * sizeof(unsigned));
    assert(idx_first_bit + n_bits_in_number < 8 * sizeof(IntegralT));
    const IntegralT mask = (static_cast<IntegralT>(1u) << n_bits_in_number) - static_cast<IntegralT>(1u);
    const unsigned result = (n & (mask << idx_first_bit)) >> idx_first_bit;
    return result;
}

template <typename IntegralT>
void set_chunk_number(IntegralT& n, unsigned char idx_first_bit, unsigned char n_bits_in_number, unsigned n_chunk) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    assert(n_bits_in_number + 1u < 8 * sizeof(unsigned));
    assert(idx_first_bit + n_bits_in_number < 8 * sizeof(IntegralT));
    const IntegralT mask = (static_cast<IntegralT>(1u) << n_bits_in_number) - static_cast<IntegralT>(1u);
    n &= ~(mask << idx_first_bit);
    n |= (n_chunk << idx_first_bit);
}

template <typename IntegralT>
IntegralT rotate(IntegralT n, unsigned char n_all_bits, unsigned char idx_pivot_bit) noexcept {
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
IntegralT refine(IntegralT n, bool new_value, unsigned char idx_bit) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    const IntegralT mask = (static_cast<IntegralT>(1u) << idx_bit);
    const IntegralT result = (new_value ? n | mask : n & ~mask);
    return result;
}

template <typename IntegralT>
auto integral_to_bits_range(IntegralT n, unsigned char n_all_bits) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    assert(n_all_bits < 8 * sizeof(IntegralT));
    const auto r = boost::irange(static_cast<unsigned char>(0u), n_all_bits) | boost::adaptors::transformed(
                [n](unsigned char idx_bit) -> bool {return extract_bit(n, idx_bit);}
    );
    return r;
}

template <typename IntegralT, typename RangeT>
IntegralT integral_from_bits_range(RangeT r) noexcept {
    static_assert(std::is_arithmetic_v<IntegralT>);
    static_assert(std::is_integral_v<IntegralT>);
    static_assert(std::is_unsigned_v<IntegralT>);
    using IteratorT = typename boost::range_iterator<RangeT>::type;
    using DerefIteratorT = decltype(*std::declval<IteratorT>());
    using RemoveCvRefDerefIteratorT = std::remove_cv_t<std::remove_reference_t<DerefIteratorT>>; //  std::remove_cvref_t<DerefIteratorIntegralT>
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

} // end of namespace kstate_op_integral::raw
