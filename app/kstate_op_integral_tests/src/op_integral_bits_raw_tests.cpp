#include <kstate_op_integral/op_integral_bits_raw.hpp>

#include <boost/range.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cstdint> //for types like: uint64_t

TEST(KstateOpIntegralRaw, BitsExtractBit) {
    EXPECT_FALSE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 0));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 1));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 2));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 3));
    EXPECT_FALSE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 4));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 5));
    EXPECT_FALSE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 6));
    EXPECT_FALSE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 7));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 8));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110, 9));
}

TEST(KstateOpIntegralRaw, BitsRotated) {
    ASSERT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110, 10, 0), 0b1100101110);
    ASSERT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110, 10, 1), 0b0110010111);
    ASSERT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110, 10, 2), 0b1011001011);
    ASSERT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110, 10, 9), 0b1001011101);
}

TEST(KstateOpIntegralRaw, BitsRefined) {
    ASSERT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110, false, 3), 0b1100100110);
    ASSERT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110, true, 3), 0b1100101110);
    ASSERT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110, false, 4), 0b1100101110);
    ASSERT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110, true, 4), 0b1100111110);
}

TEST(KstateOpIntegralRaw, ExtractChunkNumber) {
    ASSERT_EQ(kstate_op_integral::raw::extract_chunk_number<uint64_t>(0b1100101110, 0, 4), 0b1110);
    ASSERT_EQ(kstate_op_integral::raw::extract_chunk_number<uint64_t>(0b1100101110, 1, 4), 0b0111);
    ASSERT_EQ(kstate_op_integral::raw::extract_chunk_number<uint64_t>(0b1100101110, 2, 4), 0b1011);
    ASSERT_EQ(kstate_op_integral::raw::extract_chunk_number<uint64_t>(0b1100101110, 3, 4), 0b0101);
    ASSERT_EQ(kstate_op_integral::raw::extract_chunk_number<uint64_t>(0b1100101110, 4, 4), 0b0010);
    ASSERT_EQ(kstate_op_integral::raw::extract_chunk_number<uint64_t>(0b1100101110, 5, 4), 0b1001);
    ASSERT_EQ(kstate_op_integral::raw::extract_chunk_number<uint64_t>(0b1100101110, 6, 4), 0b1100);
}

TEST(KstateOpIntegralRaw, SetChunkNumber) {
    {
      uint64_t n = 0b0;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 0, 5, 0b10101);
      ASSERT_EQ(n, 0b10101);
    }
    {
      uint64_t n = 0b0;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 0, 5, 0b01010);
      ASSERT_EQ(n, 0b01010);
    }
    {
      uint64_t n = 0b111111111;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 0, 5, 0b10101);
      ASSERT_EQ(n, 0b111110101);
    }
    {
      uint64_t n = 0b111111111;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 0, 5, 0b01010);
      ASSERT_EQ(n, 0b111101010);
    }
    {
      uint64_t n = 0b0;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 1, 5, 0b10101);
      ASSERT_EQ(n, 0b101010);
    }
    {
      uint64_t n = 0b0;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 1, 5, 0b01010);
      ASSERT_EQ(n, 0b010100);
    }
    {
      uint64_t n = 0b111111111;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 1, 5, 0b10101);
      ASSERT_EQ(n, 0b111101011);
    }
    {
      uint64_t n = 0b111111111;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 1, 5, 0b01010);
      ASSERT_EQ(n, 0b111010101);
    }
    {
      //           0b9876543210
      //           0b-11010----
      uint64_t n = 0b1100101110;
      kstate_op_integral::raw::set_chunk_number<uint64_t>(n, 4, 5, 0b11010);
      ASSERT_EQ(n, 0b1110101110);
    }
}


//template <typename IntegralT>
//void set_chunk_number(IntegralT& n, unsigned char idx_first_bit, unsigned char n_bits_in_number, unsigned n_chunk) noexcept {
//    static_assert(std::is_arithmetic_v<IntegralT>);
//    static_assert(std::is_integral_v<IntegralT>);
//    static_assert(std::is_unsigned_v<IntegralT>);
//    assert(n_bits_in_number + 1 < 8 * sizeof(unsigned));
//    assert(idx_first_bit + n_bits_in_number < 8 * sizeof(IntegralT));
//    const IntegralT mask = static_cast<IntegralT>((1u << n_bits_in_number) - 1) << idx_first_bit;
//    n &= ~mask;
//    n |= n_chunk << n_bits_in_number;
//    //TODO tests!
//}



TEST(KstateOpIntegralRaw, ToBitsRange) {
    const auto r = kstate_op_integral::raw::integral_to_bits_range<uint64_t>(0b1100101110, 10);
    ASSERT_EQ(boost::size(r), 10);
    ASSERT_EQ(*std::next(std::begin(r), 0), false);
    ASSERT_EQ(*std::next(std::begin(r), 1), true);
    ASSERT_EQ(*std::next(std::begin(r), 2), true);
    ASSERT_EQ(*std::next(std::begin(r), 3), true);
    ASSERT_EQ(*std::next(std::begin(r), 4), false);
    ASSERT_EQ(*std::next(std::begin(r), 5), true);
    ASSERT_EQ(*std::next(std::begin(r), 6), false);
    ASSERT_EQ(*std::next(std::begin(r), 7), false);
    ASSERT_EQ(*std::next(std::begin(r), 8), true);
    ASSERT_EQ(*std::next(std::begin(r), 9), true);
}

TEST(KstateOpIntegralRaw, FromBitsRange) {
    using BufferT = typename std::array<bool, 10>;
    //using IteratorT = typename BufferT::iterator;
    //using ConstIteratorT = typename BufferT::const_iterator;
    //using RangeT = typename boost::iterator_range<IteratorT>;
    //using ConstRangeT = typename boost::iterator_range<ConstIteratorT>;
    BufferT a {false, true, true, true, false, true, false, false, true, true};
    const auto n = kstate_op_integral::raw::integral_from_bits_range<uint64_t>(a);
    ASSERT_EQ(n, 0b1100101110);
}
