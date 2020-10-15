#include<kstate_op_integral/op_integral_least_replication_shift.hpp>
/*
#include <boost/range/algorithm.hpp>

#include <array>

#include <gtest/gtest.h>


TEST(KstateOpIntegral, LeastReplicationShiftTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    EXPECT_EQ(kstate_op_integral::n_least_replication_shift(v1), 6);
}

TEST(KstateOpIntegral, LeastReplicationShiftTest1) {
    const std::array<int, 6> v2 = {11, 12, 13, 11, 12, 13};
    EXPECT_EQ(kstate_op_integral::n_least_replication_shift(v2), 3);
}

TEST(KstateOpIntegral, LeastReplicationShiftTest2) {
    const std::array<int, 6> v3 = {11, 12, 11, 12, 11, 12};
    EXPECT_EQ(kstate_op_integral::n_least_replication_shift(v3), 2);
}

TEST(KstateOpIntegral, LeastReplicationShiftTest3) {
    const std::array<int, 6> v4 = {11, 11, 11, 11, 11, 11};
    EXPECT_EQ(kstate_op_integral::n_least_replication_shift(v4), 1);
}

TEST(KstateOpIntegral, NormFactor0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    double expected_norm = 1.0 / std::sqrt(6);
    EXPECT_DOUBLE_EQ(kstate_op_integral::norm_factor(v1), expected_norm);
}

TEST(KstateOpIntegral, NormFactor1) {
    const std::array<int, 6> v2 = {11, 12, 13, 11, 12, 13};
    double expected_norm = 1.0 / std::sqrt(3) / 2;
    EXPECT_DOUBLE_EQ(kstate_op_integral::norm_factor(v2), expected_norm);
}

TEST(KstateOpIntegral, NormFactor2) {
    const std::array<int, 6> v3 = {11, 12, 11, 12, 11, 12};
    double expected_norm = 1.0 / std::sqrt(2) / 3;
    EXPECT_DOUBLE_EQ(kstate_op_integral::norm_factor(v3), expected_norm);
}

TEST(KstateOpIntegral, NormFactor3) {
    const std::array<int, 6> v4 = {11, 11, 11, 11, 11, 11};
    double expected_norm = 1.0 / 6.0;
    EXPECT_DOUBLE_EQ(kstate_op_integral::norm_factor(v4), expected_norm);
}

TEST(KstateOpIntegral, IsProlificTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    EXPECT_TRUE(kstate_op_integral::is_prolific(v1, 0));
    EXPECT_TRUE(kstate_op_integral::is_prolific(v1, 1));
    EXPECT_TRUE(kstate_op_integral::is_prolific(v1, 2));
    EXPECT_TRUE(kstate_op_integral::is_prolific(v1, 3));
    EXPECT_TRUE(kstate_op_integral::is_prolific(v1, 4));
    EXPECT_TRUE(kstate_op_integral::is_prolific(v1, 5));
}

TEST(KstateOpIntegral, IsProlificTest1) {
    const std::array<int, 6> v2 = {11, 11, 11, 11, 11, 11};
    EXPECT_TRUE(kstate_op_integral::is_prolific(v2, 0));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v2, 1));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v2, 2));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v2, 3));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v2, 4));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v2, 5));
}

TEST(KstateOpIntegral, IsProlificTest2) {
    const std::array<int, 6> v3 = {11, 12, 13, 11, 12, 13};
    EXPECT_TRUE(kstate_op_integral::is_prolific(v3, 0));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v3, 1));
    EXPECT_TRUE(kstate_op_integral::is_prolific(v3, 2));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v3, 3));
    EXPECT_TRUE(kstate_op_integral::is_prolific(v3, 4));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v3, 5));
}

TEST(KstateOpIntegral, IsProlificTest3) {
    const std::array<int, 6> v4 = {11, 12, 11, 12, 11, 12};
    EXPECT_TRUE(kstate_op_integral::is_prolific(v4, 0));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v4, 1));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v4, 2));
    EXPECT_TRUE(kstate_op_integral::is_prolific(v4, 3));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v4, 4));
    EXPECT_FALSE(kstate_op_integral::is_prolific(v4, 5));
}
*/
