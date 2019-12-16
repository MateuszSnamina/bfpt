#include <kstate/kstate.hpp>

#include <boost/range/algorithm.hpp>

#include <array>
#include <list>
#include <vector>

#include <gtest/gtest.h>

using kstate::ctr_from_range;

TEST(DynamicKstate, ConstructorFromRange) {
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_EQ(k1.n_sites(), 6);
}

TEST(DynamicKstate, CompareTest0) {
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const int v2[6] = {13, 14, 15, 16, 11, 12};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    EXPECT_TRUE(k1.compare(k1));
    EXPECT_FALSE(k1.compare(k2));
}

TEST(DynamicKstate, TranlationalCompareTest0) {
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const int v2[6] = {13, 14, 15, 16, 11, 12};
    const int v11[6] = {16, 11, 12, 13, 14, 15};
    const int v12[6] = {15, 16, 11, 12, 13, 14};
    const int v13[6] = {14, 15, 16, 11, 12, 13};
    const int v14[6] = {13, 14, 15, 16, 11, 12};
    const int v15[6] = {12, 13, 14, 15, 16, 11};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    const kstate::DynamicKstate<int> k11(v11, ctr_from_range);
    const kstate::DynamicKstate<int> k12(v12, ctr_from_range);
    const kstate::DynamicKstate<int> k13(v13, ctr_from_range);
    const kstate::DynamicKstate<int> k14(v14, ctr_from_range);
    const kstate::DynamicKstate<int> k15(v15, ctr_from_range);
    ASSERT_FALSE(k1.compare(k2));
    ASSERT_TRUE(k1.tranlational_compare(k1));
    ASSERT_EQ(*k1.tranlational_compare(k1), 0);
    ASSERT_TRUE(k1.tranlational_compare(k11));
    ASSERT_EQ(*k1.tranlational_compare(k11), 1);
    ASSERT_TRUE(k1.tranlational_compare(k12));
    ASSERT_EQ(*k1.tranlational_compare(k12), 2);
    ASSERT_TRUE(k1.tranlational_compare(k13));
    ASSERT_EQ(*k1.tranlational_compare(k13), 3);
    ASSERT_TRUE(k1.tranlational_compare(k13));
    ASSERT_EQ(*k1.tranlational_compare(k14), 4);
    ASSERT_TRUE(k1.tranlational_compare(k14));
    ASSERT_EQ(*k1.tranlational_compare(k15), 5);
    ASSERT_TRUE(k1.tranlational_compare(k15));
}

TEST(DynamicKstate, LeastReplicationShiftTest0) {
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_EQ(k1.n_least_replication_shift(), 6);
}

TEST(DynamicKstate, LeastReplicationShiftTest1) {
    const int v2[6] = {11, 12, 13, 11, 12, 13};
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.n_least_replication_shift(), 3);
}

TEST(DynamicKstate, LeastReplicationShiftTest2) {
    const int v3[6] = {11, 12, 11, 12, 11, 12};
    const kstate::DynamicKstate<int> k3(v3, ctr_from_range);
    EXPECT_EQ(k3.n_least_replication_shift(), 2);
}

TEST(DynamicKstate, LeastReplicationShiftTest3) {
    const int v4[6] = {11, 11, 11, 11, 11, 11};
    const kstate::DynamicKstate<int> k4(v4, ctr_from_range);
    EXPECT_EQ(k4.n_least_replication_shift(), 1);
}

TEST(DynamicKstate, IsProlificTest0) {
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_TRUE(k1.is_prolific(0));
    EXPECT_TRUE(k1.is_prolific(1));
    EXPECT_TRUE(k1.is_prolific(2));
    EXPECT_TRUE(k1.is_prolific(3));
    EXPECT_TRUE(k1.is_prolific(4));
    EXPECT_TRUE(k1.is_prolific(5));
}

TEST(DynamicKstate, IsProlificTest1) {
    const int v2[6] = {11, 11, 11, 11, 11, 11};
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    EXPECT_TRUE(k2.is_prolific(0));
    EXPECT_FALSE(k2.is_prolific(1));
    EXPECT_FALSE(k2.is_prolific(2));
    EXPECT_FALSE(k2.is_prolific(3));
    EXPECT_FALSE(k2.is_prolific(4));
    EXPECT_FALSE(k2.is_prolific(5));
}

TEST(DynamicKstate, IsProlificTest2) {
    const int v3[6] = {11, 12, 13, 11, 12, 13};
    const kstate::DynamicKstate<int> k3(v3, ctr_from_range);
    EXPECT_TRUE(k3.is_prolific(0));
    EXPECT_FALSE(k3.is_prolific(1));
    EXPECT_TRUE(k3.is_prolific(2));
    EXPECT_FALSE(k3.is_prolific(3));
    EXPECT_TRUE(k3.is_prolific(4));
    EXPECT_FALSE(k3.is_prolific(5));
}

TEST(DynamicKstate, IsProlificTest3) {
    const int v4[6] = {11, 12, 11, 12, 11, 12};
    const kstate::DynamicKstate<int> k4(v4, ctr_from_range);
    EXPECT_TRUE(k4.is_prolific(0));
    EXPECT_FALSE(k4.is_prolific(1));
    EXPECT_FALSE(k4.is_prolific(2));
    EXPECT_TRUE(k4.is_prolific(3));
    EXPECT_FALSE(k4.is_prolific(4));
    EXPECT_FALSE(k4.is_prolific(5));
}

TEST(DynamicKstate, ToStrTest0) {
    const int v1[6] = {11, 12, 13, 14, 15, 16};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    EXPECT_EQ(k1.to_str(), "⦃11∙12∙13∙14∙15∙16⦄");
}

TEST(DynamicKstate, ToStrTest1) {
    const int v2[1] = {11};
    const kstate::DynamicKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

// #######################################################################
// ## UniqueDynamicKstate                                               ##
// #######################################################################

TEST(DynamicUniqueKstate, CtrTest0) {
    const int v2[1] = {11};
    const kstate::DynamicUniqueKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

TEST(DynamicUniqueKstate, CtrTest1) {
    const int v2[2] = {11, 12};
    const kstate::DynamicUniqueKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(DynamicUniqueKstate, CtrTest2) {
    const int v2[2] = {12, 11};
    const kstate::DynamicUniqueKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(DynamicUniqueKstate, CtrTest3) {
    const int v2[7] = {12, 11, 14, 13, 14, 14, 13};
    const kstate::DynamicUniqueKstate<int> k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃14∙14∙13∙12∙11∙14∙13⦄");
}

// #######################################################################
// ## KstateUniqueView                                                  ##
// #######################################################################

TEST(KstateUniqueView, OfFilledDynamicKstate) {
    const int v1[7] = {14, 13, 14, 14, 13, 12, 11};
    const kstate::DynamicKstate<int> k1(v1, ctr_from_range);
    const auto k2 = kstate::KstateUniqueView(k1);
    EXPECT_EQ(k1.to_str(), "⦃14∙13∙14∙14∙13∙12∙11⦄");
    EXPECT_EQ(k2.to_str(), "⦃14∙14∙13∙12∙11∙14∙13⦄");
}
