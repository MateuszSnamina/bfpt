#include <kstate_op_integral/op_integral_bits_raw.hpp>

#include <boost/range.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cstdint> //for types like: uint64_t

TEST(KstateOpIntegral, BitsExtractBit) {
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

TEST(KstateOpIntegral, BitsRotated) {
    ASSERT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110, 10, 0), 0b1100101110);
    ASSERT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110, 10, 1), 0b0110010111);
    ASSERT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110, 10, 2), 0b1011001011);
    ASSERT_EQ(kstate_op_integral::raw::rotate<uint64_t>(0b1100101110, 10, 9), 0b1001011101);
}

TEST(KstateOpIntegral, BitsRefined) {
    ASSERT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110, false, 3), 0b1100100110);
    ASSERT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110, true, 3), 0b1100101110);
    ASSERT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110, false, 4), 0b1100101110);
    ASSERT_EQ(kstate_op_integral::raw::refine<uint64_t>(0b1100101110, true, 4), 0b1100111110);
}

TEST(KstateOpIntegral, ToBitsRange) {
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

TEST(KstateOpIntegral, FromBitsRange) {
    using BufferT = typename std::array<bool, 10>;
    //using IteratorT = typename BufferT::iterator;
    //using ConstIteratorT = typename BufferT::const_iterator;
    //using RangeT = typename boost::iterator_range<IteratorT>;
    //using ConstRangeT = typename boost::iterator_range<ConstIteratorT>;
    BufferT a {false, true, true, true, false, true, false, false, true, true};
    const auto n = kstate_op_integral::raw::integral_from_bits_range<uint64_t>(a);
    ASSERT_EQ(n, 0b1100101110);
}


