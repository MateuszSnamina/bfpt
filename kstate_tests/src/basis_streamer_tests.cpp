#include <kstate/basis_streamer.hpp>

#include <gtest/gtest.h>

using kstate::ctr_from_range;

TEST(BasisStreamer, TESTTRY) {
    // not yet a test...
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
    const auto k0 = std::make_shared<kstate::DynamicKstate<int>>(v0, ctr_from_range);
    const auto k1 = std::make_shared<kstate::DynamicKstate<int>>(v1, ctr_from_range);
    const auto k2 = std::make_shared<kstate::DynamicKstate<int>>(v2, ctr_from_range);
    const auto k3 = std::make_shared<kstate::DynamicKstate<int>>(v3, ctr_from_range);
    const auto k4 = std::make_shared<kstate::DynamicKstate<int>>(v4, ctr_from_range);
    const auto k5 = std::make_shared<kstate::DynamicKstate<int>>(v5, ctr_from_range);
    const auto k6 = std::make_shared<kstate::DynamicKstate<int>>(v6, ctr_from_range);
    const auto k7 = std::make_shared<kstate::DynamicKstate<int>>(v7, ctr_from_range);
    const auto k8 = std::make_shared<kstate::DynamicKstate<int>>(v8, ctr_from_range);
    const auto k9 = std::make_shared<kstate::DynamicKstate<int>>(v9, ctr_from_range);
    const auto k10 = std::make_shared<kstate::DynamicKstate<int>>(v10, ctr_from_range);
    const auto k11 = std::make_shared<kstate::DynamicKstate<int>>(v11, ctr_from_range);
    kstate::Basis<kstate::DynamicKstate<int>> basis(3);
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
    kstate::BasisStreamer bs(std::cout);
    bs.stream(basis);
}