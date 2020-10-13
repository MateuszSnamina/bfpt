#include<kstate_range_op/range_op_least_replication_shift.hpp>

#include <boost/range/algorithm.hpp>

#include <array>

#include <gtest/gtest.h>

/*
// #######################################################################
// ## UniqueKstateRangeOp                                               ##
// #######################################################################

TEST(KstateRangeOp, CtrTest0) {
    const std::array<int, 1> v2 = {11};
    const Kstate k2(kstate_range_op::make_unique_shift(v2), ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

TEST(KstateRangeOp, CtrTest1) {
    const std::array<int, 2> v2 = {11, 12};
    const Kstate k2(kstate_range_op::make_unique_shift(v2), ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(KstateRangeOp, CtrTest2) {
    const std::array<int, 2> v2 = {12, 11};
    const Kstate k2(kstate_range_op::make_unique_shift(v2), ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(KstateRangeOp, CtrTest3) {
    const std::array<int, 7> v2 = {12, 11, 14, 13, 14, 14, 13};
    const Kstate k2(kstate_range_op::make_unique_shift(v2), ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃14∙14∙13∙12∙11∙14∙13⦄");
}
*/
