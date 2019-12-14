#include <kstate/kstate.hpp>

#include <array>
#include <boost/range/algorithm.hpp>
#include <list>
#include <vector>

#include <gtest/gtest.h>

using kstate::from_range;

TEST(SimpleKstate, ConstructorFromRange) {
  const int v1[6] = {11, 12, 13, 14, 15, 16};
  const kstate::SimpleKstate<int> k1(v1, from_range);
  EXPECT_TRUE(boost::equal(k1.to_range(), v1));
  EXPECT_EQ(k1.n_sites(), 6);
}

TEST(SimpleKstate, CompareTest0) {
  const int v1[6] = {11, 12, 13, 14, 15, 16};
  const int v2[6] = {13, 14, 15, 16, 11, 12};
  const kstate::SimpleKstate<int> k1(v1, from_range);
  const kstate::SimpleKstate<int> k2(v2, from_range);
  EXPECT_TRUE(k1.compare(k1));
  EXPECT_FALSE(k1.compare(k2));
}

TEST(SimpleKstate, TranlationalCompareTest0) {
  const int v1[6] = {11, 12, 13, 14, 15, 16};
  const int v2[6] = {13, 14, 15, 16, 11, 12};
  const int v11[6] = {16, 11, 12, 13, 14, 15};
  const int v12[6] = {15, 16, 11, 12, 13, 14};
  const int v13[6] = {14, 15, 16, 11, 12, 13};
  const int v14[6] = {13, 14, 15, 16, 11, 12};
  const int v15[6] = {12, 13, 14, 15, 16, 11};
  const kstate::SimpleKstate<int> k1(v1, from_range);
  const kstate::SimpleKstate<int> k2(v2, from_range);
  const kstate::SimpleKstate<int> k11(v11, from_range);
  const kstate::SimpleKstate<int> k12(v12, from_range);
  const kstate::SimpleKstate<int> k13(v13, from_range);
  const kstate::SimpleKstate<int> k14(v14, from_range);
  const kstate::SimpleKstate<int> k15(v15, from_range);
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

TEST(SimpleKstate, LeastReplicationShiftTest0) {
  const int v1[6] = {11, 12, 13, 14, 15, 16};
  const kstate::SimpleKstate<int> k1(v1, from_range);
  EXPECT_EQ(k1.n_least_replication_shift(), 6);
}

TEST(SimpleKstate, LeastReplicationShiftTest1) {
  const int v2[6] = {11, 12, 13, 11, 12, 13};
  const kstate::SimpleKstate<int> k2(v2, from_range);
  EXPECT_EQ(k2.n_least_replication_shift(), 3);
}

TEST(SimpleKstate, LeastReplicationShiftTest2) {
  const int v3[6] = {11, 12, 11, 12, 11, 12};
  const kstate::SimpleKstate<int> k3(v3, from_range);
  EXPECT_EQ(k3.n_least_replication_shift(), 2);
}

TEST(SimpleKstate, LeastReplicationShiftTest3) {
  const int v4[6] = {11, 11, 11, 11, 11, 11};
  const kstate::SimpleKstate<int> k4(v4, from_range);
  EXPECT_EQ(k4.n_least_replication_shift(), 1);
}

TEST(SimpleKstate, IsProlificTest0) {
  const int v1[6] = {11, 12, 13, 14, 15, 16};
  const kstate::SimpleKstate<int> k1(v1, from_range);
  EXPECT_TRUE(k1.is_prolific(0));
  EXPECT_TRUE(k1.is_prolific(1));
  EXPECT_TRUE(k1.is_prolific(2));
  EXPECT_TRUE(k1.is_prolific(3));
  EXPECT_TRUE(k1.is_prolific(4));
  EXPECT_TRUE(k1.is_prolific(5));
}

TEST(SimpleKstate, IsProlificTest1) {
  const int v2[6] = {11, 11, 11, 11, 11, 11};
  const kstate::SimpleKstate<int> k2(v2, from_range);
  EXPECT_TRUE(k2.is_prolific(0));
  EXPECT_FALSE(k2.is_prolific(1));
  EXPECT_FALSE(k2.is_prolific(2));
  EXPECT_FALSE(k2.is_prolific(3));
  EXPECT_FALSE(k2.is_prolific(4));
  EXPECT_FALSE(k2.is_prolific(5));
}

TEST(SimpleKstate, IsProlificTest2) {
  const int v3[6] = {11, 12, 13, 11, 12, 13};
  const kstate::SimpleKstate<int> k3(v3, from_range);
  EXPECT_TRUE(k3.is_prolific(0));
  EXPECT_FALSE(k3.is_prolific(1));
  EXPECT_TRUE(k3.is_prolific(2));
  EXPECT_FALSE(k3.is_prolific(3));
  EXPECT_TRUE(k3.is_prolific(4));
  EXPECT_FALSE(k3.is_prolific(5));
}

TEST(SimpleKstate, IsProlificTest3) {
  const int v4[6] = {11, 12, 11, 12, 11, 12};
  const kstate::SimpleKstate<int> k4(v4, from_range);
  EXPECT_TRUE(k4.is_prolific(0));
  EXPECT_FALSE(k4.is_prolific(1));
  EXPECT_FALSE(k4.is_prolific(2));
  EXPECT_TRUE(k4.is_prolific(3));
  EXPECT_FALSE(k4.is_prolific(4));
  EXPECT_FALSE(k4.is_prolific(5));
}

TEST(SimpleKstate, ToStrTest0) {
  const int v1[6] = {11, 12, 13, 14, 15, 16};
  const kstate::SimpleKstate<int> k1(v1, from_range);
  EXPECT_EQ(k1.to_str(), "⦃11∙12∙13∙14∙15∙16⦄");
}

TEST(SimpleKstate, ToStrTest1) {
  const int v2[1] = {11};
  const kstate::SimpleKstate<int> k2(v2, from_range);
  EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

TEST(UniqueKstate, CtrTest0) {
  const int v2[1] = {11};
  const kstate::SimpleUniqueKstate<int> k2(v2, from_range);
  EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

TEST(UniqueKstate, CtrTest1) {
  const int v2[2] = {11, 12};
  const kstate::SimpleUniqueKstate<int> k2(v2, from_range);
  EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(UniqueKstate, CtrTest2) {
  const int v2[2] = {12, 11};
  const kstate::SimpleUniqueKstate<int> k2(v2, from_range);
  EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(UniqueKstate, CtrTest3) {
  const int v2[7] = {12, 11, 14, 13, 14, 14, 13};
  const kstate::SimpleUniqueKstate<int> k2(v2, from_range);
  EXPECT_EQ(k2.to_str(), "⦃14∙14∙13∙12∙11∙14∙13⦄");
}