#include <kstate_impl_tests/site_state_trait_for_monostar.hpp>
#include <kstate_impl_tests/site_state_trait_for_int.hpp>
#include <kstate_impl_tests/site_state_trait_for_my_site.h>

#include <kstate_impl/kstate_concrete_integral.hpp>

#include <array>
#include <iterator>

#include <gtest/gtest.h>

TEST(DynamicIntegralKstate, IntegralToSiteStateRange) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> integral_bits{0b01101, 5u};
    const auto r = kstate_impl::helpers::integral_bits_to_site_state_range<MonostarSiteStateTrait, uint64_t>(integral_bits, 1u);
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
        kstate_impl::helpers::integral_bits_from_site_state_range<MonostarSiteStateTrait, uint64_t>(siste_state_array, 1u);
    EXPECT_EQ(integra_bits.get_n_all_bits(), 5);
    EXPECT_EQ(integra_bits.get_number(), 0b11001);
}

// ***************************************************
// ** DynamicIntegralKstate -- one bit per site     **
// ***************************************************

TEST(DynamicIntegralKstate, ConstructorFromRange0) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MonostarSiteStateTrait, 1u>;
    const std::array<MonostarSiteState, 0> v1;
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 0);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

TEST(DynamicIntegralKstate, ConstructorFromRange1) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MonostarSiteStateTrait, 1u>;
    const std::array<MonostarSiteState, 6> v1 = {gs, gs, gs, gs, gs, gs};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 0), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 1), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 2), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 3), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 4), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 5), gs);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

TEST(DynamicIntegralKstate, ConstructorFromRange2) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MonostarSiteStateTrait, 1u>;
    const std::array<MonostarSiteState, 6> v1 = {gs, es, es, gs, es, gs};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 0), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 1), es);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 2), es);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 3), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 4), es);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 5), gs);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

// ***********************************************************
// ** DynamicIntegralKstate -- three bits per site, part 1  **
// ***********************************************************

TEST(DynamicIntegralKstate, ConstructorFromRange3) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MonostarSiteStateTrait, 3u>;
    const std::array<MonostarSiteState, 0> v1;
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 0);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

TEST(DynamicIntegralKstate, ConstructorFromRange4) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MonostarSiteStateTrait, 3u>;
    const std::array<MonostarSiteState, 6> v1 = {gs, gs, gs, gs, gs, gs};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 0), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 1), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 2), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 3), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 4), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 5), gs);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

TEST(DynamicIntegralKstate, ConstructorFromRange5) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MonostarSiteStateTrait, 3u>;
    const std::array<MonostarSiteState, 6> v1 = {gs, es, es, gs, es, gs};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 0), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 1), es);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 2), es);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 3), gs);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 4), es);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 5), gs);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

TEST(DynamicIntegralKstate, ConstructorFromRange6) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 0> v1;
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 0);
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

// ***********************************************************
// ** DynamicIntegralKstate -- three bits per site, part 2  **
// ***********************************************************

TEST(DynamicIntegralKstate, ConstructorFromRange7) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 5> v1 = {MySiteState(12), MySiteState(12), MySiteState(12), MySiteState(12), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 5);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 0), MySiteState(12));
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 1), MySiteState(12));
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 2), MySiteState(12));
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 3), MySiteState(12));
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 4), MySiteState(12));
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

TEST(DynamicIntegralKstate, ConstructorFromRange8) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 5> v1 = {MySiteState(14), MySiteState(12), MySiteState(13), MySiteState(13), MySiteState(14)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 5);
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 0), MySiteState(14));
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 1), MySiteState(12));
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 2), MySiteState(13));
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 3), MySiteState(13));
    EXPECT_EQ(*std::next(std::begin(k1.to_range()), 4), MySiteState(14));
    EXPECT_TRUE(boost::equal(k1.to_range(), v1));
}

TEST(DynamicIntegralKstate, CompareEqualityTest0) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(11), MySiteState(12), MySiteState(13), MySiteState(14), MySiteState(14), MySiteState(11)};
    const std::array<MySiteState, 6> v2 = {MySiteState(13), MySiteState(14), MySiteState(14), MySiteState(11), MySiteState(11), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    const Kstate k2(v2, ctr_from_range);
    EXPECT_TRUE(k1.compare_equality_range(k1.to_range()));
    EXPECT_FALSE(k1.compare_equality_range(k2.to_range()));
}

TEST(DynamicIntegralKstate, LeastReplicationShiftTest0) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_EQ(k1.n_least_replication_shift(), 1);
}

TEST(DynamicIntegralKstate, LeastReplicationShiftTest1) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(13), MySiteState(14), MySiteState(12), MySiteState(13)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_EQ(k1.n_least_replication_shift(), 3);
}

TEST(DynamicIntegralKstate, LeastReplicationShiftTest2) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_EQ(k1.n_least_replication_shift(), 2);
}

TEST(DynamicIntegralKstate, LeastReplicationShiftTest3) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_EQ(k1.n_least_replication_shift(), 6);
}

TEST(DynamicIntegralKstate, NormFactor0) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    const double expected_norm = 1.0 / 6;
    EXPECT_DOUBLE_EQ(k1.norm_factor(), expected_norm);
}

TEST(DynamicIntegralKstate, NormFactor1) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(13), MySiteState(14), MySiteState(12), MySiteState(13)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    const double expected_norm = 1.0 / std::sqrt(3) / 2;
    EXPECT_DOUBLE_EQ(k1.norm_factor(), expected_norm);
}

TEST(DynamicIntegralKstate, NormFactor2) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    const double expected_norm = 1.0 / std::sqrt(2) / 3;
    EXPECT_DOUBLE_EQ(k1.norm_factor(), expected_norm);
}

TEST(DynamicIntegralKstate, NormFactor3) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    const double expected_norm = 1.0 / std::sqrt(6);
    EXPECT_DOUBLE_EQ(k1.norm_factor(), expected_norm);
}

TEST(DynamicIntegralKstate, IsProlificTest0) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14), MySiteState(14)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_TRUE(k1.is_prolific(0));
    EXPECT_FALSE(k1.is_prolific(1));
    EXPECT_FALSE(k1.is_prolific(2));
    EXPECT_FALSE(k1.is_prolific(3));
    EXPECT_FALSE(k1.is_prolific(4));
    EXPECT_FALSE(k1.is_prolific(5));
}

TEST(DynamicIntegralKstate, IsProlificTest1) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(13), MySiteState(14), MySiteState(12), MySiteState(13)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_TRUE(k1.is_prolific(0));
    EXPECT_FALSE(k1.is_prolific(1));
    EXPECT_TRUE(k1.is_prolific(2));
    EXPECT_FALSE(k1.is_prolific(3));
    EXPECT_TRUE(k1.is_prolific(4));
    EXPECT_FALSE(k1.is_prolific(5));
}

TEST(DynamicIntegralKstate, IsProlificTest2) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_TRUE(k1.is_prolific(0));
    EXPECT_FALSE(k1.is_prolific(1));
    EXPECT_FALSE(k1.is_prolific(2));
    EXPECT_TRUE(k1.is_prolific(3));
    EXPECT_FALSE(k1.is_prolific(4));
    EXPECT_FALSE(k1.is_prolific(5));
}

TEST(DynamicIntegralKstate, IsProlificTest3) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    ASSERT_EQ(k1.n_sites(), 6);
    ASSERT_TRUE(boost::equal(k1.to_range(), v1));
    EXPECT_TRUE(k1.is_prolific(0));
    EXPECT_TRUE(k1.is_prolific(1));
    EXPECT_TRUE(k1.is_prolific(2));
    EXPECT_TRUE(k1.is_prolific(3));
    EXPECT_TRUE(k1.is_prolific(4));
    EXPECT_TRUE(k1.is_prolific(5));
}

TEST(DynamicIntegralKstate, ToStrTest0) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 1> v1 = {MySiteState(11)};
    const Kstate k1(v1, ctr_from_range);
    //MySiteState(10) => 'K'
    //MySiteState(11) => 'L'
    //MySiteState(12) => 'M'
    //MySiteState(13) => 'N'
    //MySiteState(14) => 'O'
    //MySiteState(15) => 'P'
    EXPECT_EQ(k1.to_str(), "⦃L⦄");
}

TEST(DynamicIntegralKstate, ToStrTest1) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    const std::array<MySiteState, 6> v1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const Kstate k1(v1, ctr_from_range);
    EXPECT_EQ(k1.to_str(), "⦃O∙M∙L∙L∙M∙M⦄");
}

//// #######################################################################
//// ## UniqueDynamicIntegralKstate                                       ##
//// #######################################################################

//TEST(UniqueDynamicIntegralKstate, CtrTest0) {
//    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
//    const std::array<MySiteState, 1> a1 = {MySiteState(11)};
//    const Kstate k1(a1, ctr_from_range);
//    //const Kstate k2(kstate_op_range::make_unique_shift(a2), ctr_from_range);
//    EXPECT_EQ(k1.to_str(), "⦃L⦄");
//}

//TEST(UniqueDynamicIntegralKstate, CtrTest1) {
//    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
//    const std::array<MySiteState, 2> a1 = {MySiteState(11), MySiteState(12)};
//    const Kstate k1(a1, ctr_from_range);
//    //    const Kstate k2(kstate_op_range::make_unique_shift(a1), ctr_from_range);
//    EXPECT_EQ(k1.to_str(), "⦃L∙M⦄");
//}

//TEST(UniqueDynamicIntegralKstate, CtrTest2) {
//    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
//    const std::array<MySiteState, 2> a1 = {MySiteState(12), MySiteState(11)};
//    const Kstate k1(a1, ctr_from_range);
//    //    const Kstate k1(kstate_op_range::make_unique_shift(a1), ctr_from_range);
//    EXPECT_EQ(k1.to_str(), "⦃L∙M⦄");
//}

//TEST(UniqueDynamicUniqueKstate, CtrTest3) {
//    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
//    const std::array<MySiteState, 7> a1 = {MySiteState(12), MySiteState(11), MySiteState(14), MySiteState(13), MySiteState(14), MySiteState(14), MySiteState(13)};
//    const Kstate k1(a1, ctr_from_range);//    using SiteStateTrait = kstate_trait::TraitSiteState<int>;
//    //    const Kstate k2(kstate_op_range::make_unique_shift(a1), ctr_from_range);
//    EXPECT_EQ(k1.to_str(), "⦃N∙M∙L∙O∙N∙O∙O⦄");
//}
