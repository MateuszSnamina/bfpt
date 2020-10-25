#include<kstate_op_integral/op_integral_bits_least_replication_shift.hpp>
#include<kstate_op_integral/integral_bits_buffer.hpp>

#include <gtest/gtest.h>

#include <cstdint> //for types like: uint64_t

TEST(KstateOpIntegral, LeastReplicationShiftTest0) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b0{0b010011, 6};
    EXPECT_EQ(kstate_op_integral::n_least_replication_shift(b0, 1), 6);
}

TEST(KstateOpIntegral, LeastReplicationShiftTest1) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b110110, 6};
    EXPECT_EQ(kstate_op_integral::n_least_replication_shift(b1, 1), 3);
}

TEST(KstateOpIntegral, LeastReplicationShiftTest2) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b2{0b010101, 6};
    EXPECT_EQ(kstate_op_integral::n_least_replication_shift(b2, 1), 2);
}

TEST(KstateOpIntegral, LeastReplicationShiftTest3) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b3{0b111111, 6};
    EXPECT_EQ(kstate_op_integral::n_least_replication_shift(b3, 1), 1);
}


TEST(KstateOpIntegral, NormFactor0) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b0{0b010011, 6};
    double expected_norm = 1.0 / std::sqrt(6);
    EXPECT_DOUBLE_EQ(kstate_op_integral::norm_factor(b0, 1), expected_norm);
}

TEST(KstateOpIntegral, NormFactor1) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b110110, 6};
    double expected_norm = 1.0 / std::sqrt(3) / 2;
    EXPECT_DOUBLE_EQ(kstate_op_integral::norm_factor(b1, 1), expected_norm);
}

TEST(KstateOpIntegral, NormFactor2) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b2{0b010101, 6};
    double expected_norm = 1.0 / std::sqrt(2) / 3;
    EXPECT_DOUBLE_EQ(kstate_op_integral::norm_factor(b2, 1), expected_norm);
}

TEST(KstateOpIntegral, NormFactor3) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b3{0b111111, 6};
    double expected_norm = 1.0 / 6.0;
    EXPECT_DOUBLE_EQ(kstate_op_integral::norm_factor(b3, 1), expected_norm);
}

TEST(KstateOpIntegral, IsProlificTest0) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b0{0b010011, 6};
    EXPECT_TRUE(kstate_op_integral::is_prolific(b0, 1, 0));
    EXPECT_TRUE(kstate_op_integral::is_prolific(b0, 1, 1));
    EXPECT_TRUE(kstate_op_integral::is_prolific(b0, 1, 2));
    EXPECT_TRUE(kstate_op_integral::is_prolific(b0, 1, 3));
    EXPECT_TRUE(kstate_op_integral::is_prolific(b0, 1, 4));
    EXPECT_TRUE(kstate_op_integral::is_prolific(b0, 1, 5));
}

TEST(KstateOpIntegral, IsProlificTest1) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b110110, 6};
    EXPECT_TRUE(kstate_op_integral::is_prolific(b1, 1, 0));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b1, 1, 1));
    EXPECT_TRUE(kstate_op_integral::is_prolific(b1, 1, 2));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b1, 1, 3));
    EXPECT_TRUE(kstate_op_integral::is_prolific(b1, 1, 4));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b1, 1, 5));
}

TEST(KstateOpIntegral, IsProlificTest2) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b2{0b010101, 6};
    EXPECT_TRUE(kstate_op_integral::is_prolific(b2, 1, 0));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b2, 1, 1));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b2, 1, 2));
    EXPECT_TRUE(kstate_op_integral::is_prolific(b2, 1, 3));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b2, 1, 4));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b2, 1, 5));
}

TEST(KstateOpIntegral, IsProlificTest3) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b3{0b111111, 6};
    EXPECT_TRUE(kstate_op_integral::is_prolific(b3, 1, 0));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b3, 1, 1));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b3, 1, 2));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b3, 1, 3));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b3, 1, 4));
    EXPECT_FALSE(kstate_op_integral::is_prolific(b3, 1, 5));
}
