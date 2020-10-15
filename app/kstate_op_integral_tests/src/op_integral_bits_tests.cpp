#include <kstate_op_integral/op_integral_bits.hpp>

#include <gtest/gtest.h>

#include <cstdint> //for types like: uint64_t

TEST(KstateOpBits, ExtractBit) {
    EXPECT_FALSE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 0));
    EXPECT_TRUE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 1));
    EXPECT_TRUE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 2));
    EXPECT_TRUE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 3));
    EXPECT_FALSE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 4));
    EXPECT_TRUE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 5));
    EXPECT_FALSE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 6));
    EXPECT_FALSE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 7));
    EXPECT_TRUE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 8));
    EXPECT_TRUE(kstate_op_integral::extract_bit<uint64_t>(0b1100101110, 9));
}

TEST(ExtensionBoostAdaptorsRotated, FilledRawArray) {
    ASSERT_EQ(kstate_op_integral::rotate<uint64_t>(0b1100101110, 10, 0), 0b1100101110);
    ASSERT_EQ(kstate_op_integral::rotate<uint64_t>(0b1100101110, 10, 1), 0b0110010111);
    ASSERT_EQ(kstate_op_integral::rotate<uint64_t>(0b1100101110, 10, 2), 0b1011001011);
    ASSERT_EQ(kstate_op_integral::rotate<uint64_t>(0b1100101110, 10, 9), 0b1001011101);
}
