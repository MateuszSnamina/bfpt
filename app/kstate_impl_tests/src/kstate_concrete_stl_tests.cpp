#include <kstate_impl_tests/site_state_trait_for_int.hpp>

#include <kstate_op_range/op_range_unique_shift.hpp>

#include <kstate_impl/kstate_concrete_stl.hpp>

#include <boost/range/algorithm.hpp>

#include <array>
#include <list>
#include <vector>

#include <gtest/gtest.h>

using kstate_impl::ctr_from_range;

TEST(DynamicStlKstate, ConstructorFromRange) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const KstateT k1(v1, ctr_from_range);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_EQ(k1.n_sites(), 6);
}

TEST(DynamicStlKstate, CompareEqualityTest0) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const std::array<int, 6> v2 = {13, 14, 15, 16, 11, 12};
    const KstateT k1(v1, ctr_from_range);
    const KstateT k2(v2, ctr_from_range);
    EXPECT_TRUE(k1.compare_equality_range(k1.to_range()));
    EXPECT_FALSE(k1.compare_equality_range(k2.to_range()));
}

TEST(DynamicStlKstate, CompareTranlationalEqualityTest0) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const std::array<int, 6> v2 = {13, 14, 15, 16, 11, 12};
    const std::array<int, 6> v11 = {16, 11, 12, 13, 14, 15};
    const std::array<int, 6> v12 = {15, 16, 11, 12, 13, 14};
    const std::array<int, 6> v13 = {14, 15, 16, 11, 12, 13};
    const std::array<int, 6> v14 = {13, 14, 15, 16, 11, 12};
    const std::array<int, 6> v15 = {12, 13, 14, 15, 16, 11};
    const KstateT k1(v1, ctr_from_range);
    const KstateT k2(v2, ctr_from_range);
    const KstateT k11(v11, ctr_from_range);
    const KstateT k12(v12, ctr_from_range);
    const KstateT k13(v13, ctr_from_range);
    const KstateT k14(v14, ctr_from_range);
    const KstateT k15(v15, ctr_from_range);
    //ASSERT_FALSE(k1.compare_translational_equality_range(v2)); //TODO rethink!
    ASSERT_TRUE(k1.compare_translational_equality_range(v1));
    ASSERT_EQ(*k1.compare_translational_equality_range(v1), 0);
    ASSERT_TRUE(k1.compare_translational_equality_range(v11));
    ASSERT_EQ(*k1.compare_translational_equality_range(v11), 1);
    ASSERT_TRUE(k1.compare_translational_equality_range(v12));
    ASSERT_EQ(*k1.compare_translational_equality_range(v12), 2);
    ASSERT_TRUE(k1.compare_translational_equality_range(v13));
    ASSERT_EQ(*k1.compare_translational_equality_range(v13), 3);
    ASSERT_TRUE(k1.compare_translational_equality_range(v14));
    ASSERT_EQ(*k1.compare_translational_equality_range(v14), 4);
    ASSERT_TRUE(k1.compare_translational_equality_range(v15));
    ASSERT_EQ(*k1.compare_translational_equality_range(v15), 5);
}

TEST(DynamicStlKstate, LeastReplicationShiftTest0) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const KstateT k1(v1, ctr_from_range);
    EXPECT_EQ(k1.n_least_replication_shift(), 6);
}

TEST(DynamicStlKstate, LeastReplicationShiftTest1) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v2 = {11, 12, 13, 11, 12, 13};
    const KstateT k2(v2, ctr_from_range);
    EXPECT_EQ(k2.n_least_replication_shift(), 3);
}

TEST(DynamicStlKstate, LeastReplicationShiftTest2) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v3 = {11, 12, 11, 12, 11, 12};
    const KstateT k3(v3, ctr_from_range);
    EXPECT_EQ(k3.n_least_replication_shift(), 2);
}

TEST(DynamicStlKstate, LeastReplicationShiftTest3) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v4 = {11, 11, 11, 11, 11, 11};
    const KstateT k4(v4, ctr_from_range);
    EXPECT_EQ(k4.n_least_replication_shift(), 1);
}

//---
TEST(DynamicStlKstate, NormFactor0) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const KstateT k1(v1, ctr_from_range);
    double expected_norm = 1.0 / std::sqrt(6);
    EXPECT_DOUBLE_EQ(k1.norm_factor(), expected_norm);
}

TEST(DynamicStlKstate, NormFactor1) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v2 = {11, 12, 13, 11, 12, 13};
    const KstateT k2(v2, ctr_from_range);
    double expected_norm = 1.0 / std::sqrt(3) / 2;
    EXPECT_DOUBLE_EQ(k2.norm_factor(), expected_norm);
}

TEST(DynamicStlKstate, NormFactor2) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v3 = {11, 12, 11, 12, 11, 12};
    const KstateT k3(v3, ctr_from_range);
    double expected_norm = 1.0 / std::sqrt(2) / 3;
    EXPECT_DOUBLE_EQ(k3.norm_factor(), expected_norm);
}

TEST(DynamicStlKstate, NormFactor3) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v4 = {11, 11, 11, 11, 11, 11};
    const KstateT k4(v4, ctr_from_range);
    double expected_norm = 1.0 / 6.0;
    EXPECT_DOUBLE_EQ(k4.norm_factor(), expected_norm);
}

TEST(DynamicStlKstate, IsProlificTest0) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const KstateT k1(v1, ctr_from_range);
    EXPECT_TRUE(k1.is_prolific(0));
    EXPECT_TRUE(k1.is_prolific(1));
    EXPECT_TRUE(k1.is_prolific(2));
    EXPECT_TRUE(k1.is_prolific(3));
    EXPECT_TRUE(k1.is_prolific(4));
    EXPECT_TRUE(k1.is_prolific(5));
}

TEST(DynamicStlKstate, IsProlificTest1) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v2 = {11, 11, 11, 11, 11, 11};
    const KstateT k2(v2, ctr_from_range);
    EXPECT_TRUE(k2.is_prolific(0));
    EXPECT_FALSE(k2.is_prolific(1));
    EXPECT_FALSE(k2.is_prolific(2));
    EXPECT_FALSE(k2.is_prolific(3));
    EXPECT_FALSE(k2.is_prolific(4));
    EXPECT_FALSE(k2.is_prolific(5));
}

TEST(DynamicStlKstate, IsProlificTest2) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v3 = {11, 12, 13, 11, 12, 13};
    const KstateT k3(v3, ctr_from_range);
    EXPECT_TRUE(k3.is_prolific(0));
    EXPECT_FALSE(k3.is_prolific(1));
    EXPECT_TRUE(k3.is_prolific(2));
    EXPECT_FALSE(k3.is_prolific(3));
    EXPECT_TRUE(k3.is_prolific(4));
    EXPECT_FALSE(k3.is_prolific(5));
}

TEST(DynamicStlKstate, IsProlificTest3) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v4 = {11, 12, 11, 12, 11, 12};
    const KstateT k4(v4, ctr_from_range);
    EXPECT_TRUE(k4.is_prolific(0));
    EXPECT_FALSE(k4.is_prolific(1));
    EXPECT_FALSE(k4.is_prolific(2));
    EXPECT_TRUE(k4.is_prolific(3));
    EXPECT_FALSE(k4.is_prolific(4));
    EXPECT_FALSE(k4.is_prolific(5));
}

TEST(DynamicStlKstate, ToStrTest0) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
    const KstateT k1(v1, ctr_from_range);
    EXPECT_EQ(k1.to_str(), "⦃11∙12∙13∙14∙15∙16⦄");
}

TEST(DynamicStlKstate, ToStrTest1) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 1> v2 = {11};
    const KstateT k2(v2, ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

// #######################################################################
// ## UniqueDynamicStlKstate                                            ##
// #######################################################################

TEST(DynamicUniqueStlKstate, CtrTest0) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 1> v2 = {11};
    const KstateT k2(kstate_op_range::make_unique_shift(v2), ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃11⦄");
}

TEST(DynamicUniqueStlKstate, CtrTest1) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 2> v2 = {11, 12};
    const KstateT k2(kstate_op_range::make_unique_shift(v2), ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(DynamicUniqueStlKstate, CtrTest2) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 2> v2 = {12, 11};
    const KstateT k2(kstate_op_range::make_unique_shift(v2), ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
}

TEST(DynamicUniqueStlKstate, CtrTest3) {
    using SiteStateTraitT = kstate_trait::TraitSiteState<int>;
    using KstateT = kstate_impl::DynamicStlKstate<SiteStateTraitT>;
    const std::array<int, 7> v2 = {12, 11, 14, 13, 14, 14, 13};
    const KstateT k2(kstate_op_range::make_unique_shift(v2), ctr_from_range);
    EXPECT_EQ(k2.to_str(), "⦃14∙14∙13∙12∙11∙14∙13⦄");
}
