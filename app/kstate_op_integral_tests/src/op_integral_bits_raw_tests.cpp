#include <kstate_op_integral/op_integral_bits_raw.hpp>

#include <boost/range.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cstdint> //for types like: uint64_t

TEST(KstateOpIntegralRaw, BitsExtractBit) {
    EXPECT_FALSE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 0));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 1));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 2));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 3));
    EXPECT_FALSE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 4));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 5));
    EXPECT_FALSE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 6));
    EXPECT_FALSE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 7));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 8));
    EXPECT_TRUE(kstate_op_integral::raw::extract_bit<uint64_t>(0b1100101110u, 9));
}

TEST(KstateOpIntegralRaw, BitsRotated) {
    EXPECT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110u, 10, 0), 0b1100101110u);
    EXPECT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110u, 10, 1), 0b0110010111u);
    EXPECT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110u, 10, 2), 0b1011001011u);
    EXPECT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110u, 10, 9), 0b1001011101u);
}

TEST(KstateOpIntegralRaw, BitsRefined) {
    EXPECT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110u, false, 3), 0b1100100110u);
    EXPECT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110u, true, 3), 0b1100101110u);
    EXPECT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110u, false, 4), 0b1100101110u);
    EXPECT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110u, true, 4), 0b1100111110u);
}

TEST(KstateOpIntegralRaw, ExtractChunkNumber) {
    EXPECT_EQ((kstate_op_integral::raw::extract_chunk_number<uint64_t, unsigned>(0b1100101110u, 0, 4)), 0b1110u);
    EXPECT_EQ((kstate_op_integral::raw::extract_chunk_number<uint64_t, unsigned>(0b1100101110u, 1, 4)), 0b0111u);
    EXPECT_EQ((kstate_op_integral::raw::extract_chunk_number<uint64_t, unsigned>(0b1100101110u, 2, 4)), 0b1011u);
    EXPECT_EQ((kstate_op_integral::raw::extract_chunk_number<uint64_t, unsigned>(0b1100101110u, 3, 4)), 0b0101u);
    EXPECT_EQ((kstate_op_integral::raw::extract_chunk_number<uint64_t, unsigned>(0b1100101110u, 4, 4)), 0b0010u);
    EXPECT_EQ((kstate_op_integral::raw::extract_chunk_number<uint64_t, unsigned>(0b1100101110u, 5, 4)), 0b1001u);
    EXPECT_EQ((kstate_op_integral::raw::extract_chunk_number<uint64_t, unsigned>(0b1100101110u, 6, 4)), 0b1100u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber0) {
    uint64_t n = 0b0u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 0, 5, 0b10101u);
    EXPECT_EQ(n, 0b10101u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber1) {
    uint64_t n = 0b0u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 0, 5, 0b01010u);
    EXPECT_EQ(n, 0b01010u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber2) {
    uint64_t n = 0b111111111u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 0, 5, 0b10101u);
    EXPECT_EQ(n, 0b111110101u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber3) {
    uint64_t n = 0b111111111u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 0, 5, 0b01010u);
    ASSERT_EQ(n, 0b111101010u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber4) {
    uint64_t n = 0b0u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 1, 5, 0b10101u);
    EXPECT_EQ(n, 0b101010u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber5) {
    uint64_t n = 0b0u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 1, 5, 0b01010u);
    EXPECT_EQ(n, 0b010100u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber6) {
    uint64_t n = 0b111111111u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 1, 5, 0b10101u);
    EXPECT_EQ(n, 0b111101011u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber7) {
    uint64_t n = 0b111111111u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 1, 5, 0b01010u);
    EXPECT_EQ(n, 0b111010101u);
}

TEST(KstateOpIntegralRaw, SetChunkNumber8) {
    //           0b9876543210
    //           0b-11010----
    uint64_t n = 0b1100101110u;
    kstate_op_integral::raw::set_chunk_number<uint64_t, unsigned>(n, 4, 5, 0b11010u);
    ASSERT_EQ(n, 0b1110101110u);
}

TEST(KstateOpIntegralRaw, ToBitsRange0) {
    const auto r = kstate_op_integral::raw::integral_to_bits_range<uint64_t>(0b0u, 0);
    ASSERT_EQ(boost::size(r), 0);
}

TEST(KstateOpIntegralRaw, ToBitsRange1) {
    const auto r = kstate_op_integral::raw::integral_to_bits_range<uint64_t>(0b1u, 1);
    ASSERT_EQ(boost::size(r), 1);
    EXPECT_EQ(*std::next(std::begin(r), 0), true);
}

TEST(KstateOpIntegralRaw, ToBitsRange2) {
    const auto r = kstate_op_integral::raw::integral_to_bits_range<uint64_t>(0b1100101110u, 10);
    ASSERT_EQ(boost::size(r), 10);
    EXPECT_EQ(*std::next(std::begin(r), 0), false);
    EXPECT_EQ(*std::next(std::begin(r), 1), true);
    EXPECT_EQ(*std::next(std::begin(r), 2), true);
    EXPECT_EQ(*std::next(std::begin(r), 3), true);
    EXPECT_EQ(*std::next(std::begin(r), 4), false);
    EXPECT_EQ(*std::next(std::begin(r), 5), true);
    EXPECT_EQ(*std::next(std::begin(r), 6), false);
    EXPECT_EQ(*std::next(std::begin(r), 7), false);
    EXPECT_EQ(*std::next(std::begin(r), 8), true);
    EXPECT_EQ(*std::next(std::begin(r), 9), true);
}

TEST(KstateOpIntegralRaw, FromBitsRange0) {
    using BufferT = typename std::array<bool, 0>;
    BufferT bits_array{};
    const auto n = kstate_op_integral::raw::integral_from_bits_range<uint64_t>(bits_array);
    EXPECT_EQ(n, 0b0u);
}

TEST(KstateOpIntegralRaw, FromBitsRange1) {
    using BufferT = typename std::array<bool, 10>;
    BufferT bits_array{false, true, true, true, false, true, false, false, true, true};
    const auto n = kstate_op_integral::raw::integral_from_bits_range<uint64_t>(bits_array);
    EXPECT_EQ(n, 0b1100101110u);
}

TEST(KstateOpIntegralRaw, ToChunkNumbersRange0) {
    const auto r = kstate_op_integral::raw::integral_to_chunk_numbers_range<uint64_t, unsigned>(0b1011001, 7, 7);
    ASSERT_EQ(boost::size(r), 1);
    EXPECT_EQ(*std::next(std::begin(r), 0), 0b1011001);
}

TEST(KstateOpIntegralRaw, ToChunkNumbersRange1) {
    const auto r = kstate_op_integral::raw::integral_to_chunk_numbers_range<uint64_t, unsigned>(0b110010111000101u, 15, 3);
    ASSERT_EQ(boost::size(r), 5);
    EXPECT_EQ(*std::next(std::begin(r), 0), 0b101u);
    EXPECT_EQ(*std::next(std::begin(r), 1), 0b000u);
    EXPECT_EQ(*std::next(std::begin(r), 2), 0b111u);
    EXPECT_EQ(*std::next(std::begin(r), 3), 0b010u);
    EXPECT_EQ(*std::next(std::begin(r), 4), 0b110u);
}

TEST(KstateOpIntegralRaw, FromChunkNumbersRange0) {
    using BufferT = typename std::array<unsigned, 0>;
    BufferT a;
    const auto n = kstate_op_integral::raw::integral_from_chunk_numbers_range<uint64_t, unsigned>(a, 3);
    EXPECT_EQ(n, 0b0u);
}

TEST(KstateOpIntegralRaw, FromChunkNumbersRange1) {
    using BufferT = typename std::array<unsigned, 1>;
    BufferT a = {0b1011001};
    const auto n = kstate_op_integral::raw::integral_from_chunk_numbers_range<uint64_t, unsigned>(a, 7);
    EXPECT_EQ(n, 0b1011001u);
}

TEST(KstateOpIntegralRaw, FromChunkNumbersRange) {
    using BufferT = typename std::array<unsigned, 5>;
    BufferT a {0b101, 0b000, 0b111, 0b010, 0b110};
    const auto n = kstate_op_integral::raw::integral_from_chunk_numbers_range<uint64_t, unsigned>(a, 3);
    EXPECT_EQ(n, 0b110010111000101u);
}
