#include <kstate/kstate.hpp>

#include <array>
#include <list>
#include <vector>

#include <boost/range/algorithm.hpp>

#include <gtest/gtest.h>

TEST(Kstate, compare) {
  int v1[6] = {11, 12, 13, 14, 15, 16};
  int v2[6] = {13, 14, 15, 16, 11, 12};

  const auto k1 = kstate::SimpleKstate<int>(v1);
  const auto k2 = kstate::SimpleKstate<int>(v2);
  EXPECT_TRUE(k1.compare(k1));
  EXPECT_FALSE(k1.compare(k2));
}

TEST(Kstate, tranlational_compare) {
  int v1[6] = {11, 12, 13, 14, 15, 16};
  int v2[6] = {13, 14, 15, 16, 11, 12};
  int v11[6] = {16, 11, 12, 13, 14, 15};
  int v12[6] = {15, 16, 11, 12, 13, 14};
  int v13[6] = {14, 15, 16, 11, 12, 13};
  int v14[6] = {13, 14, 15, 16, 11, 12};
  int v15[6] = {12, 13, 14, 15, 16, 11};

  const auto k1 = kstate::SimpleKstate<int>(v1);
  const auto k2 = kstate::SimpleKstate<int>(v2);
  const auto k11 = kstate::SimpleKstate<int>(v11);
  const auto k12 = kstate::SimpleKstate<int>(v12);
  const auto k13 = kstate::SimpleKstate<int>(v13);
  const auto k14 = kstate::SimpleKstate<int>(v14);
  const auto k15 = kstate::SimpleKstate<int>(v15);

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

TEST(Kstate, LeastReplicationShift) {
  int v1[6] = {11, 12, 13, 14, 15, 16};
  const auto k1 = kstate::SimpleKstate<int>(v1);
  EXPECT_EQ(k1.n_least_replication_shift(), 6);

  int v2[6] = {11, 12, 13, 11, 12, 13};
  const auto k2 = kstate::SimpleKstate<int>(v2);
  EXPECT_EQ(k2.n_least_replication_shift(), 3);

  int v3[6] = {11, 12, 11, 12, 11, 12};
  const auto k3 = kstate::SimpleKstate<int>(v3);
  EXPECT_EQ(k3.n_least_replication_shift(), 2);

  int v4[6] = {11, 11, 11, 11, 11, 11};
  const auto k4 = kstate::SimpleKstate<int>(v4);
  EXPECT_EQ(k4.n_least_replication_shift(), 1);
}