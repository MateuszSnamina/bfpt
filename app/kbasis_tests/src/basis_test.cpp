#include <kbasis_tests/site_state_trait_for_int.hpp>

#include <kstate_trait/trait_site_state.hpp>
#include <kstate_trait/trait_kstate.hpp>

#include <kbasis/basis.hpp>

#include <kstate_impl/kstate_concrete_stl.hpp>

#include <kstate_op_range/op_unique_shift.hpp>

#include <boost/range/algorithm.hpp>

#include <memory>

#include <gtest/gtest.h>

using kstate_impl::ctr_from_range;

TEST(Basis, Empty) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using KStateTrait = kstate_trait::TraitKstate<Kstate>;
    kbasis::Basis<KStateTrait> basis(3);
    // test basis size:
    ASSERT_EQ(basis.size(), 0);
    ASSERT_TRUE(basis.vec_index().begin() == basis.vec_index().end());
    ASSERT_TRUE(basis.map_index().begin() == basis.map_index().end());//TODO remove
    // test not-finding a not existing element:
    const int v100[3] = {11, 12, 113};
    const auto k100 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v100, ctr_from_range);
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k100->to_range()));
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v100));
}

TEST(Basis, OneElement) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using KStateTrait = kstate_trait::TraitKstate<Kstate>;
    kbasis::Basis<KStateTrait> basis(3);
    const int v1[3] = {11, 12, 13};
    const auto k1 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v1, ctr_from_range);
    basis.add_element(k1);
    // test basis size:
    ASSERT_EQ(basis.size(), 1);
    ASSERT_TRUE(basis.vec_index().begin() + 1 == basis.vec_index().end());
    ASSERT_TRUE(std::next(basis.map_index().begin(), 1) == basis.map_index().end());
    // test element access:
    EXPECT_TRUE(*basis.vec_index().begin() == k1);
    EXPECT_TRUE(*basis.map_index().begin() == k1);
    // test finding the existing element:
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k1->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v1));
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k1->to_range()), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v1), 0);
    const int v10[3] = {11, 12, 113};
    const auto k10 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v10, ctr_from_range);
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k10->to_range()));
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v10));
    // test not-finding a not existing element:
    const int v100[3] = {11, 12, 113};
    const auto k100 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v100, ctr_from_range);
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k100->to_range()));
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v100));
}

TEST(Basis, TwoDifferentElementsTest0) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using KStateTrait = kstate_trait::TraitKstate<Kstate>;
    kbasis::Basis<KStateTrait> basis(3);
    const int v1[3] = {11, 12, 13};
    const int v2[3] = {13, 14, 15};
    const auto k1 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v1, ctr_from_range);
    const auto k2 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v2, ctr_from_range);
    basis.add_element(k1);
    basis.add_element(k2);
    // test basis size:
    ASSERT_EQ(basis.size(), 2);
    ASSERT_TRUE(basis.vec_index().begin() + 2 == basis.vec_index().end());
    ASSERT_TRUE(std::next(basis.map_index().begin(), 2) == basis.map_index().end());
    // test element access:
    EXPECT_TRUE(*(basis.vec_index().begin() + 0) == k1);
    EXPECT_TRUE(*(basis.vec_index().begin() + 1) == k2);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k1);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 1) == k2);
    // test finding the existing elements:
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k1->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v1));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k2->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v2));
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k1->to_range()), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v1), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k2->to_range()), 1);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v2), 1);
    // test not-finding a not existing element:
    const int v100[3] = {11, 12, 113};
    const auto k100 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v100, ctr_from_range);
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k100->to_range()));
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v100));
}

TEST(Basis, TwoDifferentElementsTest1) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using KStateTrait = kstate_trait::TraitKstate<Kstate>;
    kbasis::Basis<KStateTrait> basis(3);
    const int v1[3] = {11, 12, 13};
    const int v2[3] = {13, 14, 15};
    const auto k1 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v1, ctr_from_range);
    const auto k2 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v2, ctr_from_range);
    basis.add_element(k2);
    basis.add_element(k1);
    // test basis size:
    ASSERT_EQ(basis.size(), 2);
    ASSERT_TRUE(basis.vec_index().begin() + 2 == basis.vec_index().end());
    ASSERT_TRUE(std::next(basis.map_index().begin(), 2) == basis.map_index().end());
    // test element access:
    EXPECT_TRUE(*(basis.vec_index().begin() + 0) == k2);
    EXPECT_TRUE(*(basis.vec_index().begin() + 1) == k1);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k1);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 1) == k2);
    // test finding the existing elements:
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k1->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v1));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k2->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v2));
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k1->to_range()), 1);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v1), 1);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k2->to_range()), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v2), 0);
    // test not-finding a not existing element:
    const int v100[3] = {11, 12, 113};
    const auto k100 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v100, ctr_from_range);
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k100->to_range()));
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v100));
}

TEST(Basis, TwoSameElementsTest0) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using KStateTrait = kstate_trait::TraitKstate<Kstate>;
    kbasis::Basis<KStateTrait> basis(3);
    const int v1[3] = {11, 12, 13};
    const int v2[3] = {11, 12, 13};
    const auto k1 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v1, ctr_from_range);
    const auto k2 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v2, ctr_from_range);
    basis.add_element(k1);
    basis.add_element(k2);
    // test basis size:
    ASSERT_EQ(basis.size(), 1);
    ASSERT_TRUE(basis.vec_index().begin() + 1 == basis.vec_index().end());
    ASSERT_TRUE(std::next(basis.map_index().begin(), 1) == basis.map_index().end());
    // test element access:
    EXPECT_TRUE(*(basis.vec_index().begin() + 0) == k1);
    EXPECT_FALSE(*(basis.vec_index().begin() + 0) == k2);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k1);
    EXPECT_FALSE(*std::next(basis.map_index().begin(), 0) == k2);
    // test finding the existing elements:
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k1->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v1));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k2->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v2));
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k1->to_range()), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v1), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k2->to_range()), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v2), 0);
    // test not-finding a not existing element:
    const int v100[3] = {11, 12, 113};
    const auto k100 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v100, ctr_from_range);
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k100->to_range()));
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v100));
}

TEST(Basis, BigTest) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using KStateTrait = kstate_trait::TraitKstate<Kstate>;
    kbasis::Basis<KStateTrait> basis(3);
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
    const auto k0 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v0, ctr_from_range);
    const auto k1 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v1, ctr_from_range);
    const auto k2 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v2, ctr_from_range);
    const auto k3 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v3, ctr_from_range);
    const auto k4 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v4, ctr_from_range);
    const auto k5 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v5, ctr_from_range);
    const auto k6 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v6, ctr_from_range);
    const auto k7 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v7, ctr_from_range);
    const auto k8 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v8, ctr_from_range);
    const auto k9 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v9, ctr_from_range);
    const auto k10 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v10, ctr_from_range);
    const auto k11 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v11, ctr_from_range);
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
    // test basis size:
    ASSERT_EQ(basis.size(), 9);
    ASSERT_TRUE(basis.vec_index().begin() + 9 == basis.vec_index().end());
    ASSERT_TRUE(std::next(basis.map_index().begin(), 9) == basis.map_index().end());
    // test element access:
    EXPECT_TRUE(*(basis.vec_index().begin() + 0) == k0);
    EXPECT_TRUE(*(basis.vec_index().begin() + 1) == k1);
    EXPECT_TRUE(*(basis.vec_index().begin() + 2) == k2);
    EXPECT_TRUE(*(basis.vec_index().begin() + 3) == k5);
    EXPECT_TRUE(*(basis.vec_index().begin() + 4) == k6);
    EXPECT_TRUE(*(basis.vec_index().begin() + 5) == k8);
    EXPECT_TRUE(*(basis.vec_index().begin() + 6) == k9);
    EXPECT_TRUE(*(basis.vec_index().begin() + 7) == k10);
    EXPECT_TRUE(*(basis.vec_index().begin() + 8) == k11);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k6);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 1) == k5);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 2) == k8);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 3) == k9);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 4) == k11);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 5) == k10);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 6) == k0);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 7) == k1);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 8) == k2);
    // test finding the existing elements:
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k0->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v0));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k1->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v1));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k2->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v2));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k3->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v3));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k4->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v4));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k5->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v5));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k6->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v6));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k7->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v7));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k8->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v8));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k9->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v9));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k10->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v10));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k11->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v11));
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k0->to_range()), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v0), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k1->to_range()), 1);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v1), 1);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k2->to_range()), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v2), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k3->to_range()), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v3), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k4->to_range()), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v4), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k5->to_range()), 3);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v5), 3);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k6->to_range()), 4);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v6), 4);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k7->to_range()), 3);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v7), 3);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k8->to_range()), 5);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v8), 5);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k9->to_range()), 6);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v9), 6);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k10->to_range()), 7);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v10), 7);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k11->to_range()), 8);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(v11), 8);
    // test not-finding a not existing element:
    const int v100[3] = {11, 12, 113};
    const auto k100 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(v100, ctr_from_range);
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k100->to_range()));
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v100));

    // print:
    // std::cout << "size:" << basis.size() << std::endl;
    // for (const auto& _ : basis.vec_index()) {
    //   std::cout << "vec_index: " << _->to_str() << std::endl;
    // }
    // for (const auto& _ : basis.map_index()) {
    //   std::cout << "map_index: " << _->to_str() << std::endl;
    // }
}

TEST(Basis, BigTestWithUniqueStates) {
    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
    using Kstate = kstate_impl::DynamicKstate<SiteStateTrait>;
    using KStateTrait = kstate_trait::TraitKstate<Kstate>;
    kbasis::Basis<KStateTrait> basis(3);
    const int v0[3] = {7, 12, 13};
    const int v1[3] = {11, 12, 13};
    const int v2[3] = {13, 14, 15};
    const int v3[3] = {13, 14, 15};  // (equivalent) replica of v2
    const int v4[3] = {15, 13, 14};  // (equivalent) replica of v2
    const int v5[3] = {1, 20, 15};
    const int v6[3] = {1, 14, 15};
    const int v7[3] = {20, 15, 1};  // (equivalent) replica of v5
    const int v8[3] = {3, 20, 15};
    const auto k0 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v0), ctr_from_range);
    const auto k1 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v1), ctr_from_range);
    const auto k2 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v2), ctr_from_range);
    const auto k3 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v3), ctr_from_range);
    const auto k4 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v4), ctr_from_range);
    const auto k5 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v5), ctr_from_range);
    const auto k6 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v6), ctr_from_range);
    const auto k7 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v7), ctr_from_range);
    const auto k8 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v8), ctr_from_range);
    basis.add_element(k0);
    basis.add_element(k1);
    basis.add_element(k2);
    basis.add_element(k3);
    basis.add_element(k4);
    basis.add_element(k5);
    basis.add_element(k6);
    basis.add_element(k7);
    basis.add_element(k8);
    // test basis size:
    ASSERT_EQ(basis.size(), 6);
    ASSERT_TRUE(basis.vec_index().begin() + 6 == basis.vec_index().end());
    ASSERT_TRUE(std::next(basis.map_index().begin(), 6) == basis.map_index().end());
    // test element access:
    EXPECT_TRUE(*(basis.vec_index().begin() + 0) == k0);
    EXPECT_TRUE(*(basis.vec_index().begin() + 1) == k1);
    EXPECT_TRUE(*(basis.vec_index().begin() + 2) == k2);
    EXPECT_TRUE(*(basis.vec_index().begin() + 3) == k5);
    EXPECT_TRUE(*(basis.vec_index().begin() + 4) == k6);
    EXPECT_TRUE(*(basis.vec_index().begin() + 5) == k8);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 0) == k0);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 1) == k1);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 2) == k6);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 3) == k2);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 4) == k5);
    EXPECT_TRUE(*std::next(basis.map_index().begin(), 5) == k8);
    // test finding the existing elements:
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k0->to_range()));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v0 | extension::boost::adaptors::rotated(0)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v0 | extension::boost::adaptors::rotated(1)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v0 | extension::boost::adaptors::rotated(2)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k1->to_range()));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v1 | extension::boost::adaptors::rotated(0)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v1 | extension::boost::adaptors::rotated(1)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v1 | extension::boost::adaptors::rotated(2)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k2->to_range()));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v2 | extension::boost::adaptors::rotated(0)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v2 | extension::boost::adaptors::rotated(1)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v2 | extension::boost::adaptors::rotated(2)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k3->to_range()));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v3 | extension::boost::adaptors::rotated(0)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v3 | extension::boost::adaptors::rotated(1)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v3 | extension::boost::adaptors::rotated(2)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k4->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v4 | extension::boost::adaptors::rotated(0)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v4 | extension::boost::adaptors::rotated(1)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v4 | extension::boost::adaptors::rotated(2)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k5->to_range()));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v5 | extension::boost::adaptors::rotated(0)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v5 | extension::boost::adaptors::rotated(1)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v5 | extension::boost::adaptors::rotated(2)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k6->to_range()));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v6 | extension::boost::adaptors::rotated(0)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v6 | extension::boost::adaptors::rotated(1)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v6 | extension::boost::adaptors::rotated(2)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k7->to_range()));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v7 | extension::boost::adaptors::rotated(0)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v7 | extension::boost::adaptors::rotated(1)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v7 | extension::boost::adaptors::rotated(2)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(k8->to_range()));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v8 | extension::boost::adaptors::rotated(0)));
    ASSERT_TRUE(basis.find_element_and_get_its_ra_index(v8 | extension::boost::adaptors::rotated(1)));
    ASSERT_FALSE(basis.find_element_and_get_its_ra_index(v8 | extension::boost::adaptors::rotated(2)));
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k0->to_range()), 0);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k1->to_range()), 1);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k2->to_range()), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k3->to_range()), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k4->to_range()), 2);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k5->to_range()), 3);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k6->to_range()), 4);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k7->to_range()), 3);
    EXPECT_EQ(*basis.find_element_and_get_its_ra_index(k8->to_range()), 5);
    // test not-finding a not existing element:
    const int v100[3] = {11, 12, 113};
    const auto k100 = std::make_shared<kstate_impl::DynamicKstate<kstate_trait::TraitSiteState<int>>>(kstate_op_range::make_unique_shift(v100), ctr_from_range);
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(k100->to_range()));
    EXPECT_FALSE(basis.find_element_and_get_its_ra_index(v100));

    //print:
    // std::cout << "size:" << basis.size() << std::endl;
    // for (const auto& _ : basis.vec_index()) {
    //   std::cout << "vec_index: " << _->to_str() << std::endl;
    // }
    // for (const auto& _ : basis.map_index()) {
    //   std::cout << "map_index: " << _->to_str() << std::endl;
    // }
}
