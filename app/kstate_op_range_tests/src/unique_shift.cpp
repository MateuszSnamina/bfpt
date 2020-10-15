#include<kstate_op_range/op_unique_shift.hpp>
#include<kstate_op_range/op_compare.hpp>

#include <boost/range/algorithm.hpp>

#include <array>

#include <gtest/gtest.h>

TEST(KstateRangeOp, Test0) {
    const std::array<int, 1> v2 = {11};
    const std::array<int, 1> v2u = {11};
    kstate_op_range::compare_equality(v2u, kstate_op_range::make_unique_shift(v2));
}

TEST(KstateRangeOp, Test1) {
    const std::array<int, 2> v2 = {11, 12};
    const std::array<int, 2> v2u = {12, 11};
    kstate_op_range::compare_equality(v2u, kstate_op_range::make_unique_shift(v2));
}

TEST(KstateRangeOp, Test2) {
    const std::array<int, 2> v2 = {12, 11};
    const std::array<int, 2> v2u = {12, 11};
    kstate_op_range::compare_equality(v2u, kstate_op_range::make_unique_shift(v2));
}

TEST(KstateRangeOp, Test3) {
    const std::array<int, 7> v2 = {12, 11, 14, 13, 14, 14, 13};
    const std::array<int, 7> v2u = {14, 14, 13, 12, 11, 14, 13};
    kstate_op_range::compare_equality(v2u, kstate_op_range::make_unique_shift(v2));
}
