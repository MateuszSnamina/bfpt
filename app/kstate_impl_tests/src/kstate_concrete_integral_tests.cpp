#include <kstate_impl_tests/site_state_trait_for_monostar.hpp>

#include <kstate_impl/kstate_concrete_integral.hpp>

//#include <boost/range/algorithm.hpp>

#include <array>

#include <gtest/gtest.h>

TEST(DynamicIntegralKstate, IntegralToSiteStateRange) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> integral_bits{0b01101, 5u};
    const auto r = kstate_impl::helpers::integral_bits_to_two_level_site_state_range<MonostarSiteStateTrait, uint64_t>(integral_bits);
    EXPECT_EQ(boost::size(r), 5);
    EXPECT_TRUE((*std::next(std::begin(r), 0))._is_excited);
    EXPECT_FALSE((*std::next(std::begin(r), 1))._is_excited);
    EXPECT_TRUE((*std::next(std::begin(r), 2))._is_excited);
    EXPECT_TRUE((*std::next(std::begin(r), 3))._is_excited);
    EXPECT_FALSE((*std::next(std::begin(r), 4))._is_excited);
    EXPECT_TRUE(std::next(std::begin(r), 0)->_is_excited);
    EXPECT_FALSE(std::next(std::begin(r), 1)->_is_excited);
    EXPECT_TRUE(std::next(std::begin(r), 2)->_is_excited);
    EXPECT_TRUE(std::next(std::begin(r), 3)->_is_excited);
    EXPECT_FALSE(std::next(std::begin(r), 4)->_is_excited);
}

TEST(DynamicIntegralKstate, SiteStateRangeToIntegral) {
    std::array<MonostarSiteState, 5> siste_state_array{es, gs, gs, es, es};
    kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> integra_bits =
            kstate_impl::helpers::integral_bits_from_two_level_site_state_range<MonostarSiteStateTrait, uint64_t>(siste_state_array);
    EXPECT_EQ(integra_bits.get_n_all_bits() , 5);
    EXPECT_EQ(integra_bits.get_number() , 0b11001);
}

// ****************************
// ** TESTS  ******************
// ****************************

TEST(DynamicIntegralKstate, ConstructorFromRange) {
    using Kstate = kstate_impl::DynamicTwoLevelIntegral64Kstate<MonostarSiteStateTrait>;
    const std::array<MonostarSiteState, 6> v1 = {gs, gs, gs, gs, gs, gs};
    const Kstate k1(v1, ctr_from_range);
//    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
//    EXPECT_EQ(k1.n_sites(), 6);
}

//TEST(DynamicIntegralKstate, CompareEqualityTest0) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
//    const std::array<int, 6> v2 = {13, 14, 15, 16, 11, 12};
//    const Kstate k1(v1, ctr_from_range);
//    const Kstate k2(v2, ctr_from_range);
//    EXPECT_TRUE(k1.compare_equality_range(k1.to_range()));
//    EXPECT_FALSE(k1.compare_equality_range(k2.to_range()));
//    // ---
//    EXPECT_TRUE(k1.compare_equality_any_range(k1.to_range()));
//    EXPECT_FALSE(k1.compare_equality_any_range(k2.to_range()));
//    // ---
//    EXPECT_TRUE(k1.compare_equality_kstate(k1));
//    EXPECT_FALSE(k1.compare_equality_kstate(k2));
//}

//TEST(DynamicIntegralKstate, CompareTranlationalEqualityTest0) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
//    const std::array<int, 6> v2 = {13, 14, 15, 16, 11, 12};
//    const std::array<int, 6> v11 = {16, 11, 12, 13, 14, 15};
//    const std::array<int, 6> v12 = {15, 16, 11, 12, 13, 14};
//    const std::array<int, 6> v13 = {14, 15, 16, 11, 12, 13};
//    const std::array<int, 6> v14 = {13, 14, 15, 16, 11, 12};
//    const std::array<int, 6> v15 = {12, 13, 14, 15, 16, 11};
//    const Kstate k1(v1, ctr_from_range);
//    const Kstate k2(v2, ctr_from_range);
//    const Kstate k11(v11, ctr_from_range);
//    const Kstate k12(v12, ctr_from_range);
//    const Kstate k13(v13, ctr_from_range);
//    const Kstate k14(v14, ctr_from_range);
//    const Kstate k15(v15, ctr_from_range);
//    // ---
//    //ASSERT_FALSE(k1.compare_translational_equality_range(v2)); //TODO rethink!
//    ASSERT_TRUE(k1.compare_translational_equality_range(v1));
//    ASSERT_EQ(*k1.compare_translational_equality_range(v1), 0);
//    ASSERT_TRUE(k1.compare_translational_equality_range(v11));
//    ASSERT_EQ(*k1.compare_translational_equality_range(v11), 1);
//    ASSERT_TRUE(k1.compare_translational_equality_range(v12));
//    ASSERT_EQ(*k1.compare_translational_equality_range(v12), 2);
//    ASSERT_TRUE(k1.compare_translational_equality_range(v13));
//    ASSERT_EQ(*k1.compare_translational_equality_range(v13), 3);
//    ASSERT_TRUE(k1.compare_translational_equality_range(v14));
//    ASSERT_EQ(*k1.compare_translational_equality_range(v14), 4);
//    ASSERT_TRUE(k1.compare_translational_equality_range(v15));
//    ASSERT_EQ(*k1.compare_translational_equality_range(v15), 5);
//    // ---
//    //ASSERT_FALSE(k1.compare_translational_equality_range(v2)); //TODO rethink!
//    ASSERT_TRUE(k1.compare_translational_equality_any_range(v1));
//    ASSERT_EQ(*k1.compare_translational_equality_any_range(v1), 0);
//    ASSERT_TRUE(k1.compare_translational_equality_any_range(v11));
//    ASSERT_EQ(*k1.compare_translational_equality_any_range(v11), 1);
//    ASSERT_TRUE(k1.compare_translational_equality_any_range(v12));
//    ASSERT_EQ(*k1.compare_translational_equality_any_range(v12), 2);
//    ASSERT_TRUE(k1.compare_translational_equality_any_range(v13));
//    ASSERT_EQ(*k1.compare_translational_equality_any_range(v13), 3);
//    ASSERT_TRUE(k1.compare_translational_equality_any_range(v14));
//    ASSERT_EQ(*k1.compare_translational_equality_any_range(v14), 4);
//    ASSERT_TRUE(k1.compare_translational_equality_any_range(v15));
//    ASSERT_EQ(*k1.compare_translational_equality_any_range(v15), 5);
//    // ---
//    //ASSERT_FALSE(k1.compare_translational_equality_range(k2)); //TODO rethink!
//    ASSERT_TRUE(k1.compare_translational_equality_kstate(k1));
//    ASSERT_EQ(*k1.compare_translational_equality_kstate(k1), 0);
//    ASSERT_TRUE(k1.compare_translational_equality_kstate(k11));
//    ASSERT_EQ(*k1.compare_translational_equality_kstate(k11), 1);
//    ASSERT_TRUE(k1.compare_translational_equality_kstate(k12));
//    ASSERT_EQ(*k1.compare_translational_equality_kstate(k12), 2);
//    ASSERT_TRUE(k1.compare_translational_equality_kstate(k13));
//    ASSERT_EQ(*k1.compare_translational_equality_kstate(k13), 3);
//    ASSERT_TRUE(k1.compare_translational_equality_kstate(k14));
//    ASSERT_EQ(*k1.compare_translational_equality_kstate(k14), 4);
//    ASSERT_TRUE(k1.compare_translational_equality_kstate(k15));
//    ASSERT_EQ(*k1.compare_translational_equality_kstate(k15), 5);
//}

//TEST(DynamicIntegralKstate, LeastReplicationShiftTest0) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
//    const Kstate k1(v1, ctr_from_range);
//    EXPECT_EQ(k1.n_least_replication_shift(), 6);
//    EXPECT_EQ((k1.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::n_least_replication_shift()), 6);
//}

//TEST(DynamicIntegralKstate, LeastReplicationShiftTest1) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v2 = {11, 12, 13, 11, 12, 13};
//    const Kstate k2(v2, ctr_from_range);
//    EXPECT_EQ(k2.n_least_replication_shift(), 3);
//    EXPECT_EQ((k2.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::n_least_replication_shift()), 3);
//}

//TEST(DynamicIntegralKstate, LeastReplicationShiftTest2) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v3 = {11, 12, 11, 12, 11, 12};
//    const Kstate k3(v3, ctr_from_range);
//    EXPECT_EQ(k3.n_least_replication_shift(), 2);
//    EXPECT_EQ((k3.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::n_least_replication_shift()), 2);
//}

//TEST(DynamicIntegralKstate, LeastReplicationShiftTest3) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v4 = {11, 11, 11, 11, 11, 11};
//    const Kstate k4(v4, ctr_from_range);
//    EXPECT_EQ(k4.n_least_replication_shift(), 1);
//    EXPECT_EQ((k4.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::n_least_replication_shift()), 1);
//}

////---
//TEST(DynamicIntegralKstate, NormFactor0) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
//    const Kstate k1(v1, ctr_from_range);
//    double expected_norm = 1.0 / std::sqrt(6);
//    EXPECT_DOUBLE_EQ(k1.norm_factor(), expected_norm);
//    // EXPECT_DOUBLE_EQ((k1.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::norm_factor()), expected_norm);
//}

//TEST(DynamicIntegralKstate, NormFactor1) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v2 = {11, 12, 13, 11, 12, 13};
//    const Kstate k2(v2, ctr_from_range);
//    double expected_norm = 1.0 / std::sqrt(3) / 2;
//    EXPECT_DOUBLE_EQ(k2.norm_factor(), expected_norm);
//    // EXPECT_DOUBLE_EQ((k2.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::norm_factor()), expected_norm);
//}

//TEST(DynamicIntegralKstate, NormFactor2) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v3 = {11, 12, 11, 12, 11, 12};
//    const Kstate k3(v3, ctr_from_range);
//    double expected_norm = 1.0 / std::sqrt(2) / 3;
//    EXPECT_DOUBLE_EQ(k3.norm_factor(), expected_norm);
//    // EXPECT_DOUBLE_EQ((k3.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::norm_factor()), expected_norm);
//}

//TEST(DynamicIntegralKstate, NormFactor3) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v4 = {11, 11, 11, 11, 11, 11};
//    const Kstate k4(v4, ctr_from_range);
//    double expected_norm = 1.0 / 6.0;
//    EXPECT_DOUBLE_EQ(k4.norm_factor(), expected_norm);
//    // EXPECT_DOUBLE_EQ((k4.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::norm_factor()), expected_norm);
//}

//TEST(DynamicIntegralKstate, IsProlificTest0) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
//    const Kstate k1(v1, ctr_from_range);
//    EXPECT_TRUE(k1.is_prolific(0));
//    EXPECT_TRUE(k1.is_prolific(1));
//    EXPECT_TRUE(k1.is_prolific(2));
//    EXPECT_TRUE(k1.is_prolific(3));
//    EXPECT_TRUE(k1.is_prolific(4));
//    EXPECT_TRUE(k1.is_prolific(5));
//    EXPECT_TRUE((k1.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(0)));
//    EXPECT_TRUE((k1.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(1)));
//    EXPECT_TRUE((k1.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(2)));
//    EXPECT_TRUE((k1.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(3)));
//    EXPECT_TRUE((k1.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(4)));
//    EXPECT_TRUE((k1.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(5)));
//}

//TEST(DynamicIntegralKstate, IsProlificTest1) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v2 = {11, 11, 11, 11, 11, 11};
//    const Kstate k2(v2, ctr_from_range);
//    EXPECT_TRUE(k2.is_prolific(0));
//    EXPECT_FALSE(k2.is_prolific(1));
//    EXPECT_FALSE(k2.is_prolific(2));
//    EXPECT_FALSE(k2.is_prolific(3));
//    EXPECT_FALSE(k2.is_prolific(4));
//    EXPECT_FALSE(k2.is_prolific(5));
//    EXPECT_TRUE((k2.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(0)));
//    EXPECT_FALSE((k2.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(1)));
//    EXPECT_FALSE((k2.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(2)));
//    EXPECT_FALSE((k2.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(3)));
//    EXPECT_FALSE((k2.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(4)));
//    EXPECT_FALSE((k2.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(5)));
//}

//TEST(DynamicIntegralKstate, IsProlificTest2) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v3 = {11, 12, 13, 11, 12, 13};
//    const Kstate k3(v3, ctr_from_range);
//    EXPECT_TRUE(k3.is_prolific(0));
//    EXPECT_FALSE(k3.is_prolific(1));
//    EXPECT_TRUE(k3.is_prolific(2));
//    EXPECT_FALSE(k3.is_prolific(3));
//    EXPECT_TRUE(k3.is_prolific(4));
//    EXPECT_FALSE(k3.is_prolific(5));
//    EXPECT_TRUE((k3.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(0)));
//    EXPECT_FALSE((k3.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(1)));
//    EXPECT_TRUE((k3.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(2)));
//    EXPECT_FALSE((k3.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(3)));
//    EXPECT_TRUE((k3.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(4)));
//    EXPECT_FALSE((k3.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(5)));
//}

//TEST(DynamicIntegralKstate, IsProlificTest3) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v4 = {11, 12, 11, 12, 11, 12};
//    const Kstate k4(v4, ctr_from_range);
//    EXPECT_TRUE(k4.is_prolific(0));
//    EXPECT_FALSE(k4.is_prolific(1));
//    EXPECT_FALSE(k4.is_prolific(2));
//    EXPECT_TRUE(k4.is_prolific(3));
//    EXPECT_FALSE(k4.is_prolific(4));
//    EXPECT_FALSE(k4.is_prolific(5));
//    EXPECT_TRUE((k4.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(0)));
//    EXPECT_FALSE((k4.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(1)));
//    EXPECT_FALSE((k4.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(2)));
//    EXPECT_TRUE((k4.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(3)));
//    EXPECT_FALSE((k4.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(4)));
//    EXPECT_FALSE((k4.template Kstate<kstate_trait::TraitSiteState<int>, boost::random_access_traversal_tag>::is_prolific(5)));
//}

//TEST(DynamicIntegralKstate, ToStrTest0) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 6> v1 = {11, 12, 13, 14, 15, 16};
//    const Kstate k1(v1, ctr_from_range);
//    EXPECT_EQ(k1.to_str(), "⦃11∙12∙13∙14∙15∙16⦄");
//}

//TEST(DynamicIntegralKstate, ToStrTest1) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 1> v2 = {11};
//    const Kstate k2(v2, ctr_from_range);
//    EXPECT_EQ(k2.to_str(), "⦃11⦄");
//}

//// #######################################################################
//// ## UniqueDynamicIntegralKstate                                       ##
//// #######################################################################

//TEST(DynamicUniqueKstate, CtrTest0) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 1> v2 = {11};
//    const Kstate k2(kstate_op_range::make_unique_shift(v2), ctr_from_range);
//    EXPECT_EQ(k2.to_str(), "⦃11⦄");
//}

//TEST(DynamicUniqueKstate, CtrTest1) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 2> v2 = {11, 12};
//    const Kstate k2(kstate_op_range::make_unique_shift(v2), ctr_from_range);
//    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
//}

//TEST(DynamicUniqueKstate, CtrTest2) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 2> v2 = {12, 11};
//    const Kstate k2(kstate_op_range::make_unique_shift(v2), ctr_from_range);
//    EXPECT_EQ(k2.to_str(), "⦃12∙11⦄");
//}

//TEST(DynamicUniqueKstate, CtrTest3) {
//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    using Kstate = kstate_impl::DynamicIntegralKstate<SiteStateTrait>;
//    const std::array<int, 7> v2 = {12, 11, 14, 13, 14, 14, 13};
//    const Kstate k2(kstate_op_range::make_unique_shift(v2), ctr_from_range);
//    EXPECT_EQ(k2.to_str(), "⦃14∙14∙13∙12∙11∙14∙13⦄");
//}
