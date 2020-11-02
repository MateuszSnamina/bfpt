#include <kstate_op_range/op_range_compare.hpp>

#include <array>

#include <gtest/gtest.h>

TEST(KstateOpRange, CompareEqualityTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const std::array<int, 6> v2 = {13, 14, 15, 16, 11, 12};
    EXPECT_TRUE(kstate_op_range::compare_equality(v1, v1));
    EXPECT_FALSE(kstate_op_range::compare_equality(v1, v2));
}

TEST(KstateOpRange, CompareTranlationalEqualityTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const std::array<int, 6> v2 = {13, 14, 15, 16, 11, 13};
    const std::array<int, 6> v11 = {16, 11, 12, 13, 14, 15};
    const std::array<int, 6> v12 = {15, 16, 11, 12, 13, 14};
    const std::array<int, 6> v13 = {14, 15, 16, 11, 12, 13};
    const std::array<int, 6> v14 = {13, 14, 15, 16, 11, 12};
    const std::array<int, 6> v15 = {12, 13, 14, 15, 16, 11};
    // ---
    ASSERT_FALSE(kstate_op_range::compare_translational_equality(v1, v2));
    ASSERT_TRUE(kstate_op_range::compare_translational_equality(v1, v1));
    ASSERT_EQ(kstate_op_range::compare_translational_equality(v1, v1), 0);
    ASSERT_TRUE(kstate_op_range::compare_translational_equality(v1, v11));
    ASSERT_EQ(kstate_op_range::compare_translational_equality(v1, v11), 1);
    ASSERT_TRUE(kstate_op_range::compare_translational_equality(v1, v12));
    ASSERT_EQ(kstate_op_range::compare_translational_equality(v1, v12), 2);
    ASSERT_TRUE(kstate_op_range::compare_translational_equality(v1, v13));
    ASSERT_EQ(kstate_op_range::compare_translational_equality(v1, v13), 3);
    ASSERT_TRUE(kstate_op_range::compare_translational_equality(v1, v14));
    ASSERT_EQ(kstate_op_range::compare_translational_equality(v1, v14), 4);
    ASSERT_TRUE(kstate_op_range::compare_translational_equality(v1, v15));
    ASSERT_EQ(kstate_op_range::compare_translational_equality(v1, v15), 5);
}
