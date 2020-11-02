#include <kstate_op_range/op_range_unique_shift.hpp>
#include <kstate_op_range/op_range_compare.hpp>

#include <boost/range/algorithm.hpp>

#include <array>

#include <gtest/gtest.h>

TEST(KstateOpRange, NUniqueShiftTest0) {
    const std::array<int, 1> v2 = {11};
    ASSERT_EQ(kstate_op_range::n_unique_shift(v2), 0);
}

TEST(KstateOpRange, NUniqueShiftTest1) {
    const std::array<int, 2> v2 = {11, 12};
    ASSERT_EQ(kstate_op_range::n_unique_shift(v2), 1);
}

TEST(KstateOpRange, NUniqueShiftTest2) {
    const std::array<int, 2> v2 = {12, 11};
    ASSERT_EQ(kstate_op_range::n_unique_shift(v2), 0);
}

TEST(KstateOpRange, NUniqueShiftTest3) {
    const std::array<int, 7> v2 = {12, 11, 14, 13, 14, 14, 13};
    ASSERT_EQ(kstate_op_range::n_unique_shift(v2), 4);
}

TEST(KstateOpRange, MakeUniqueShiftTest0) {
    const std::array<int, 1> v2 = {11};
    const std::array<int, 1> v2u = {11};
    ASSERT_TRUE(kstate_op_range::compare_equality(v2u, kstate_op_range::make_unique_shift(v2)));
}

TEST(KstateOpRange, MakeUniqueShiftTest1) {
    const std::array<int, 2> v2 = {11, 12};
    const std::array<int, 2> v2u = {12, 11};
    ASSERT_TRUE(kstate_op_range::compare_equality(v2u, kstate_op_range::make_unique_shift(v2)));
}

TEST(KstateOpRange, MakeUniqueShiftTest2) {
    const std::array<int, 2> v2 = {12, 11};
    const std::array<int, 2> v2u = {12, 11};
    ASSERT_TRUE(kstate_op_range::compare_equality(v2u, kstate_op_range::make_unique_shift(v2)));
}

TEST(KstateOpRange, MakeUniqueShiftTest3) {
    const std::array<int, 7> v2 = {12, 11, 14, 13, 14, 14, 13};
    const std::array<int, 7> v2u = {14, 14, 13, 12, 11, 14, 13};
    ASSERT_TRUE(kstate_op_range::compare_equality(v2u, kstate_op_range::make_unique_shift(v2)));
}
