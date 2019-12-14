#include <kstate/basis.hpp>
#include <kstate/kstate.hpp>

#include <array>
#include <boost/range/algorithm.hpp>
#include <list>
#include <memory>
#include <vector>

#include <gtest/gtest.h>

using kstate::from_range;

TEST(Basis, Empty) {
  kstate::VecMap<kstate::SimpleKstate<int>> basis;
  ASSERT_EQ(basis.size(), 0);
  ASSERT_TRUE(basis.vec_index().begin() == basis.vec_index().end());
  ASSERT_TRUE(basis.map_index().begin() == basis.map_index().end());
  const int v10[3] = {11, 12, 113};
  const auto k10 = std::make_shared<kstate::SimpleKstate<int>>(v10, from_range);
  EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k10->to_range()));
  //EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v10)); //fix me!
}

TEST(Basis, OneElement) {
  kstate::VecMap<kstate::SimpleKstate<int>> basis;
  const int v1[3] = {11, 12, 13};
  const auto k1 = std::make_shared<kstate::SimpleKstate<int>>(v1, from_range);
  basis.add_element(k1);
  ASSERT_EQ(basis.size(), 1);
  ASSERT_TRUE(std::next(basis.vec_index().begin(), 1) ==
              basis.vec_index().end());
  ASSERT_TRUE(std::next(basis.map_index().begin(), 1) ==
              basis.map_index().end());
  EXPECT_TRUE(*basis.vec_index().begin() == k1);
  EXPECT_TRUE(*basis.map_index().begin() == k1);
}

TEST(Basis, TwoDifferentElementsTest0) {
  kstate::VecMap<kstate::SimpleKstate<int>> basis;
  const int v1[3] = {11, 12, 13};
  const int v2[3] = {13, 14, 15};
  const auto k1 = std::make_shared<kstate::SimpleKstate<int>>(v1, from_range);
  const auto k2 = std::make_shared<kstate::SimpleKstate<int>>(v2, from_range);
  basis.add_element(k1);
  basis.add_element(k2);
  ASSERT_EQ(basis.size(), 2);
  ASSERT_TRUE(std::next(basis.vec_index().begin(), 2) ==
              basis.vec_index().end());
  ASSERT_TRUE(std::next(basis.map_index().begin(), 2) ==
              basis.map_index().end());
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 0) == k1);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 1) == k2);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k1);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 1) == k2);
}

TEST(Basis, TwoDifferentElementsTest1) {
  kstate::VecMap<kstate::SimpleKstate<int>> basis;
  const int v1[3] = {11, 12, 13};
  const int v2[3] = {13, 14, 15};
  const auto k1 = std::make_shared<kstate::SimpleKstate<int>>(v1, from_range);
  const auto k2 = std::make_shared<kstate::SimpleKstate<int>>(v2, from_range);
  basis.add_element(k2);
  basis.add_element(k1);
  ASSERT_EQ(basis.size(), 2);
  ASSERT_TRUE(std::next(basis.vec_index().begin(), 2) ==
              basis.vec_index().end());
  ASSERT_TRUE(std::next(basis.map_index().begin(), 2) ==
              basis.map_index().end());
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 0) == k2);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 1) == k1);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k1);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 1) == k2);
}

TEST(Basis, TwoSameElementsTest0) {
  kstate::VecMap<kstate::SimpleKstate<int>> basis;
  const int v1[3] = {11, 12, 13};
  const int v2[3] = {11, 12, 13};
  const auto k1 = std::make_shared<kstate::SimpleKstate<int>>(v1, from_range);
  const auto k2 = std::make_shared<kstate::SimpleKstate<int>>(v2, from_range);
  basis.add_element(k1);
  basis.add_element(k2);
  ASSERT_EQ(basis.size(), 1);
  ASSERT_TRUE(std::next(basis.vec_index().begin(), 1) ==
              basis.vec_index().end());
  ASSERT_TRUE(std::next(basis.map_index().begin(), 1) ==
              basis.map_index().end());
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 0) == k1);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k1);
  EXPECT_FALSE(*std::next(basis.vec_index().begin(), 0) == k2);
  EXPECT_FALSE(*std::next(basis.map_index().begin(), 0) == k2);
}

TEST(Basis, BigTest) {
  kstate::VecMap<kstate::SimpleKstate<int>> basis;

  const int v0[3] = {7, 12, 13};
  const int v1[3] = {11, 12, 13};
  const int v2[3] = {13, 14, 15};
  const int v3[3] = {13, 14, 15};  // replica of v2
  const int v4[3] = {13, 14, 15};  // replica of v2
  const int v5[3] = {1, 20, 15};
  const int v6[3] = {1, 14, 15};
  const int v7[3] = {1, 20, 15};  // replica of v5
  const int v8[3] = {3, 20, 15};
  const int v9[3] = {3, 20, 18};
  const int v10[3] = {3, 21, 15};
  const int v11[3] = {3, 21, 10};
  const auto k0 = std::make_shared<kstate::SimpleKstate<int>>(v0, from_range);
  const auto k1 = std::make_shared<kstate::SimpleKstate<int>>(v1, from_range);
  const auto k2 = std::make_shared<kstate::SimpleKstate<int>>(v2, from_range);
  const auto k3 = std::make_shared<kstate::SimpleKstate<int>>(v3, from_range);
  const auto k4 = std::make_shared<kstate::SimpleKstate<int>>(v4, from_range);
  const auto k5 = std::make_shared<kstate::SimpleKstate<int>>(v5, from_range);
  const auto k6 = std::make_shared<kstate::SimpleKstate<int>>(v6, from_range);
  const auto k7 = std::make_shared<kstate::SimpleKstate<int>>(v7, from_range);
  const auto k8 = std::make_shared<kstate::SimpleKstate<int>>(v8, from_range);
  const auto k9 = std::make_shared<kstate::SimpleKstate<int>>(v9, from_range);
  const auto k10 = std::make_shared<kstate::SimpleKstate<int>>(v10, from_range);
  const auto k11 = std::make_shared<kstate::SimpleKstate<int>>(v11, from_range);

  basis.add_element(k0);
  basis.add_element(k1);
  basis.add_element(k2);
  basis.add_element(k3);
  basis.add_element(k4);
  basis.add_element(k5);
  basis.add_element(k6);
  basis.add_element(k7);
  basis.add_element(k8);
  basis.add_element(k9);
  basis.add_element(k10);
  basis.add_element(k11);

  ASSERT_EQ(basis.size(), 9);
  ASSERT_TRUE(std::next(basis.vec_index().begin(), 9) ==
              basis.vec_index().end());
  ASSERT_TRUE(std::next(basis.map_index().begin(), 9) ==
              basis.map_index().end());

  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 0) == k0);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 1) == k1);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 2) == k2);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 3) == k5);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 4) == k6);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 5) == k8);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 6) == k9);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 7) == k10);
  EXPECT_TRUE(*std::next(basis.vec_index().begin(), 8) == k11);

  EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k6);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 1) == k5);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 2) == k8);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 3) == k9);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 4) == k11);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 5) == k10);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 6) == k0);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 7) == k1);
  EXPECT_TRUE(*std::next(basis.map_index().begin(), 8) == k2);

  // std::cout << "size:" << basis.size() << std::endl;
  // for (const auto& _ : basis.vec_index()) {
  //   std::cout << "vec_index: " << _->to_str() << std::endl;
  // }
  // for (const auto& _ : basis.map_index()) {
  //   std::cout << "map_index: " << _->to_str() << std::endl;
  // }
}
