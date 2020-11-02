#include <kstate_op_range/op_range_least_replication_shift.hpp>

#include <boost/range/algorithm.hpp>

#include <array>

#include <gtest/gtest.h>

TEST(KstateOpRange, LeastReplicationShiftTest0) {
    const std::array<int, 6> v0 = {11, 12, 13, 14, 15, 16};
    EXPECT_EQ(kstate_op_range::n_least_replication_shift(v0), 6);
}

TEST(KstateOpRange, LeastReplicationShiftTest1) {
    const std::array<int, 6> v1 = {11, 12, 13, 11, 12, 13};
    EXPECT_EQ(kstate_op_range::n_least_replication_shift(v1), 3);
}

TEST(KstateOpRange, LeastReplicationShiftTest2) {
    const std::array<int, 6> v2 = {11, 12, 11, 12, 11, 12};
    EXPECT_EQ(kstate_op_range::n_least_replication_shift(v2), 2);
}

TEST(KstateOpRange, LeastReplicationShiftTest3) {
    const std::array<int, 6> v3 = {11, 11, 11, 11, 11, 11};
    EXPECT_EQ(kstate_op_range::n_least_replication_shift(v3), 1);
}

TEST(KstateOpRange, NormFactor0) {
    const std::array<int, 6> v0 = {11, 12, 13, 14, 15, 16};
    double expected_norm = 1.0 / std::sqrt(6);
    EXPECT_DOUBLE_EQ(kstate_op_range::norm_factor(v0), expected_norm);
}

TEST(KstateOpRange, NormFactor1) {
    const std::array<int, 6> v1 = {11, 12, 13, 11, 12, 13};
    double expected_norm = 1.0 / std::sqrt(3) / 2;
    EXPECT_DOUBLE_EQ(kstate_op_range::norm_factor(v1), expected_norm);
}

TEST(KstateOpRange, NormFactor2) {
    const std::array<int, 6> v2 = {11, 12, 11, 12, 11, 12};
    double expected_norm = 1.0 / std::sqrt(2) / 3;
    EXPECT_DOUBLE_EQ(kstate_op_range::norm_factor(v2), expected_norm);
}

TEST(KstateOpRange, NormFactor3) {
    const std::array<int, 6> v3 = {11, 11, 11, 11, 11, 11};
    double expected_norm = 1.0 / 6.0;
    EXPECT_DOUBLE_EQ(kstate_op_range::norm_factor(v3), expected_norm);
}

TEST(KstateOpRange, IsProlificTest0) {
    const std::array<int, 6> v0 = {11, 12, 13, 14, 15, 16};
    EXPECT_TRUE(kstate_op_range::is_prolific(v0, 0));
    EXPECT_TRUE(kstate_op_range::is_prolific(v0, 1));
    EXPECT_TRUE(kstate_op_range::is_prolific(v0, 2));
    EXPECT_TRUE(kstate_op_range::is_prolific(v0, 3));
    EXPECT_TRUE(kstate_op_range::is_prolific(v0, 4));
    EXPECT_TRUE(kstate_op_range::is_prolific(v0, 5));
}

TEST(KstateOpRange, IsProlificTest1) {
    const std::array<int, 6> v1 = {11, 12, 13, 11, 12, 13};
    EXPECT_TRUE(kstate_op_range::is_prolific(v1, 0));
    EXPECT_FALSE(kstate_op_range::is_prolific(v1, 1));
    EXPECT_TRUE(kstate_op_range::is_prolific(v1, 2));
    EXPECT_FALSE(kstate_op_range::is_prolific(v1, 3));
    EXPECT_TRUE(kstate_op_range::is_prolific(v1, 4));
    EXPECT_FALSE(kstate_op_range::is_prolific(v1, 5));
}

TEST(KstateOpRange, IsProlificTest2) {
    const std::array<int, 6> v2 = {11, 12, 11, 12, 11, 12};
    EXPECT_TRUE(kstate_op_range::is_prolific(v2, 0));
    EXPECT_FALSE(kstate_op_range::is_prolific(v2, 1));
    EXPECT_FALSE(kstate_op_range::is_prolific(v2, 2));
    EXPECT_TRUE(kstate_op_range::is_prolific(v2, 3));
    EXPECT_FALSE(kstate_op_range::is_prolific(v2, 4));
    EXPECT_FALSE(kstate_op_range::is_prolific(v2, 5));
}

TEST(KstateOpRange, IsProlificTest3) {
    const std::array<int, 6> v3 = {11, 11, 11, 11, 11, 11};
    EXPECT_TRUE(kstate_op_range::is_prolific(v3, 0));
    EXPECT_FALSE(kstate_op_range::is_prolific(v3, 1));
    EXPECT_FALSE(kstate_op_range::is_prolific(v3, 2));
    EXPECT_FALSE(kstate_op_range::is_prolific(v3, 3));
    EXPECT_FALSE(kstate_op_range::is_prolific(v3, 4));
    EXPECT_FALSE(kstate_op_range::is_prolific(v3, 5));
}
