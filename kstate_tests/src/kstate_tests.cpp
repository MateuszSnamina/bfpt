#include <kstate/kstate.hpp>

#include <boost/range/algorithm.hpp>

#include <array>
#include <list>
#include <vector>

#include <gtest/gtest.h>

using kstate::ctr_from_range;

TEST(DynamicKstate, ConstructorFromRange) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_EQ(k1.n_sites(), 6);
}

TEST(DynamicKstate, CompareTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const std::array<int, 6> v2 = {13, 14, 15, 16, 11, 12};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    EXPECT_TRUE(k1.compare_range(k1.to_range()));
    EXPECT_FALSE(k1.compare_range(k2.to_range()));
}

TEST(DynamicKstate, TranlationalCompareTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const std::array<int, 6> v2 = {13, 14, 15, 16, 11, 12};
    const std::array<int, 6> v11 = {16, 11, 12, 13, 14, 15};
    const std::array<int, 6> v12 = {15, 16, 11, 12, 13, 14};
    const std::array<int, 6> v13 = {14, 15, 16, 11, 12, 13};
    const std::array<int, 6> v14 = {13, 14, 15, 16, 11, 12};
    const std::array<int, 6> v15 = {12, 13, 14, 15, 16, 11};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    const kstate::DynamicKstate<int> k11(v11, ctr_from_range);
    const kstate::DynamicKstate<int> k12(v12, ctr_from_range);
    const kstate::DynamicKstate<int> k13(v13, ctr_from_range);
    const kstate::DynamicKstate<int> k14(v14, ctr_from_range);
    const kstate::DynamicKstate<int> k15(v15, ctr_from_range);
    ASSERT_FALSE(k1.compare_range(v2));
    ASSERT_TRUE(k1.translational_compare_range(v1));
    ASSERT_EQ(*k1.translational_compare_range(v1), 0);
    ASSERT_TRUE(k1.translational_compare_range(v11));
    ASSERT_EQ(*k1.translational_compare_range(v11), 1);
    ASSERT_TRUE(k1.translational_compare_range(v12));
    ASSERT_EQ(*k1.translational_compare_range(v12), 2);
    ASSERT_TRUE(k1.translational_compare_range(v13));
    ASSERT_EQ(*k1.translational_compare_range(v13), 3);
    ASSERT_TRUE(k1.translational_compare_range(v13));
    ASSERT_EQ(*k1.translational_compare_range(v14), 4);
    ASSERT_TRUE(k1.translational_compare_range(v14));
    ASSERT_EQ(*k1.translational_compare_range(v15), 5);
    ASSERT_TRUE(k1.translational_compare_range(v15));
}

TEST(DynamicKstate, LeastReplicationShiftTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_EQ(k1.n_least_replication_shift(), 6);
}

TEST(DynamicKstate, LeastReplicationShiftTest1) {
    const std::array<int, 6> v2 = {11, 12, 13, 11, 12, 13};
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.n_least_replication_shift(), 3);
}

TEST(DynamicKstate, LeastReplicationShiftTest2) {
    const std::array<int, 6> v3 = {11, 12, 11, 12, 11, 12};
    const kstate::DynamicKstate<int> k3(v3, ctr_from_range);
    EXPECT_EQ(k3.n_least_replication_shift(), 2);
}

TEST(DynamicKstate, LeastReplicationShiftTest3) {
    const std::array<int, 6> v4 = {11, 11, 11, 11, 11, 11};
    const kstate::DynamicKstate<int> k4(v4, ctr_from_range);
    EXPECT_EQ(k4.n_least_replication_shift(), 1);
}

TEST(DynamicKstate, IsProlificTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_TRUE(k1.is_prolific(0));
    EXPECT_TRUE(k1.is_prolific(1));
    EXPECT_TRUE(k1.is_prolific(2));
    EXPECT_TRUE(k1.is_prolific(3));
    EXPECT_TRUE(k1.is_prolific(4));
    EXPECT_TRUE(k1.is_prolific(5));
}

TEST(DynamicKstate, IsProlificTest1) {
    const std::array<int, 6> v2 = {11, 11, 11, 11, 11, 11};
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    EXPECT_TRUE(k2.is_prolific(0));
    EXPECT_FALSE(k2.is_prolific(1));
    EXPECT_FALSE(k2.is_prolific(2));
    EXPECT_FALSE(k2.is_prolific(3));
    EXPECT_FALSE(k2.is_prolific(4));
    EXPECT_FALSE(k2.is_prolific(5));
}

TEST(DynamicKstate, IsProlificTest2) {
    const std::array<int, 6> v3 = {11, 12, 13, 11, 12, 13};
    const kstate::DynamicKstate<int> k3(v3, ctr_from_range);
    EXPECT_TRUE(k3.is_prolific(0));
    EXPECT_FALSE(k3.is_prolific(1));
    EXPECT_TRUE(k3.is_prolific(2));
    EXPECT_FALSE(k3.is_prolific(3));
    EXPECT_TRUE(k3.is_prolific(4));
    EXPECT_FALSE(k3.is_prolific(5));
}

TEST(DynamicKstate, IsProlificTest3) {
    const std::array<int, 6> v4 = {11, 12, 11, 12, 11, 12};
    const kstate::DynamicKstate<int> k4(v4, ctr_from_range);
    EXPECT_TRUE(k4.is_prolific(0));
    EXPECT_FALSE(k4.is_prolific(1));
    EXPECT_FALSE(k4.is_prolific(2));
    EXPECT_TRUE(k4.is_prolific(3));
    EXPECT_FALSE(k4.is_prolific(4));
    EXPECT_FALSE(k4.is_prolific(5));
}

TEST(DynamicKstate, ToStrTest0) {
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_EQ(k1.to_str(), "⦃11∙12∙13∙14∙15∙16⦄");
}

TEST(DynamicKstate, ToStrTest1) {
    const std::array<int, 1> v2 = {11};
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

// #######################################################################
// ## UniqueDynamicKstate                                               ##
// #######################################################################

TEST(DynamicUniqueKstate, CtrTest0) {
    const std::array<int, 1> v2 = {11};
    const kstate::DynamicUniqueKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

TEST(DynamicUniqueKstate, CtrTest1) {
    const std::array<int, 2> v2 = {11, 12};
    const kstate::DynamicUniqueKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(DynamicUniqueKstate, CtrTest2) {
    const std::array<int, 2> v2 = {12, 11};
    const kstate::DynamicUniqueKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(DynamicUniqueKstate, CtrTest3) {
    const std::array<int, 7> v2 = {12, 11, 14, 13, 14, 14, 13};
    const kstate::DynamicUniqueKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃14∙14∙13∙12∙11∙14∙13⦄");
}

// #######################################################################
// ## KstateUniqueView                                                  ##
// #######################################################################

// TEST(KstateUniqueView, OfFilledDynamicKstate) {
//     const std::array<int, 6> v1[7] = {14, 13, 14, 14, 13, 12, 11};
//     const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
//     const auto k2 = kstate::KstateUniqueView(k1);
//     EXPECT_EQ(k1.to_str(), "⦃14∙13∙14∙14∙13∙12∙11⦄");
//     EXPECT_EQ(k2.to_str(), "⦃14∙14∙13∙12∙11∙14∙13⦄");
// }
