#include <kstate/basis.hpp>
#include <kstate/kstate.hpp>

#include <gtest/gtest.h>
#include <array>
#include <boost/range/algorithm.hpp>
#include <list>
#include <memory>
#include <vector>

TEST(Basis, Fake0) {
  kstate::VecMap<kstate::SimpleKstate<int>> basis;

  int v1[3] = {11, 12, 13};
  int v2[3] = {13, 14, 15};
  int v3[3] = {13, 14, 15};
  int v4[3] = {13, 14, 15};
  int v5[3] = {1, 20, 15};
  int v6[3] = {1, 14, 15};
  int v7[3] = {1, 20, 15};
  int v8[3] = {3, 20, 15};

  const auto k1 = std::make_shared<kstate::SimpleKstate<int>>(v1, 42l);
  const auto k2 = std::make_shared<kstate::SimpleKstate<int>>(v2, 42l);
  const auto k3 = std::make_shared<kstate::SimpleKstate<int>>(v3, 42l);
  const auto k4 = std::make_shared<kstate::SimpleKstate<int>>(v4, 42l);
  const auto k5 = std::make_shared<kstate::SimpleKstate<int>>(v5, 42l);
  const auto k6 = std::make_shared<kstate::SimpleKstate<int>>(v6, 42l);
  const auto k7 = std::make_shared<kstate::SimpleKstate<int>>(v7, 42l);
  const auto k8 = std::make_shared<kstate::SimpleKstate<int>>(v8, 42l);

  std::cout << k1->to_str() << std::endl;
  basis.add_element(k1);
  basis.add_element(k2);
  basis.add_element(k3);
  basis.add_element(k4);
  basis.add_element(k5);
  basis.add_element(k6);
  basis.add_element(k7);
  basis.add_element(k8);

  //   EXPECT_TRUE(k1.compare(k1));
  //   EXPECT_FALSE(k1.compare(k2));

  std::cout << "SIZE:" << basis.size() << std::endl;
  for (const auto& _ : basis.vec_index()) {
    std::cout << _->to_str() << std::endl;
  }
  std::cout << "SIZE:" << basis.size() << std::endl;
  for (const auto& _ : basis.map_index()) {
    std::cout << _->to_str() << std::endl;
  }
}
