#include <kstate_impl_tests/site_state_trait_for_my_site.h>

#include <kstate_impl/kstate_concrete_integral.hpp>

#include <array>
#include <iterator>

#include <gtest/gtest.h>

TEST(DynamicIntegralKstateTrait, Test0) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    using KstateTrait = kstate_trait::TraitKstate<Kstate>;
    const std::array<MySiteState, 6> array_1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // before unique shift: {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // after  unique shift: {MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12), MySiteState(14)};
    const auto kstate_1 = KstateTrait::from_range(array_1);
    ASSERT_TRUE(boost::equal(kstate_1.to_range(), array_1));
    const auto view_1 = KstateTrait::to_view(kstate_1);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 0), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 1), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 2), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 3), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 4), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 5), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1), 1u);
}

TEST(DynamicIntegralKstateTrait, Test1) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    using KstateTrait = kstate_trait::TraitKstate<Kstate>;
    const std::array<MySiteState, 6> array_1 = {MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(11), MySiteState(14), MySiteState(12)};
    // before unique shift: {MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(11), MySiteState(14), MySiteState(12)};
    // after  unique shift: {MySiteState(11), MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(12), MySiteState(14)};
    const auto kstate_1 = KstateTrait::from_range(array_1);
    ASSERT_TRUE(boost::equal(kstate_1.to_range(), array_1));
    const auto view_1 = KstateTrait::to_view(kstate_1);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 0), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 1), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 2), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 3), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 4), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 5), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1), 3u);
}

TEST(DynamicIntegralKstateTrait, Test2) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    using KstateTrait = kstate_trait::TraitKstate<Kstate>;
    const std::array<MySiteState, 6> array_1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const std::array<MySiteState, 6> array_2 = {MySiteState(15), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const std::array<MySiteState, 6> array_3 = {MySiteState(14), MySiteState(13), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const std::array<MySiteState, 6> array_4 = {MySiteState(13), MySiteState(13), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(13)};
    const auto kstate_1 = KstateTrait::from_range(array_1);
    const auto kstate_2 = KstateTrait::from_range(array_2);
    const auto kstate_3 = KstateTrait::from_range(array_3);
    const auto kstate_4 = KstateTrait::from_range(array_4);
    ASSERT_TRUE(boost::equal(kstate_1.to_range(), array_1));
    ASSERT_TRUE(boost::equal(kstate_2.to_range(), array_2));
    ASSERT_TRUE(boost::equal(kstate_3.to_range(), array_3));
    ASSERT_TRUE(boost::equal(kstate_4.to_range(), array_4));
    const auto view_1 = KstateTrait::to_view(kstate_1);
    const auto view_2 = KstateTrait::to_view(kstate_2);
    const auto view_3 = KstateTrait::to_view(kstate_3);
    const auto view_4 = KstateTrait::to_view(kstate_4);
    ASSERT_TRUE(KstateTrait::view_compare_equality(view_1, view_1));
    ASSERT_FALSE(KstateTrait::view_compare_equality(view_1, view_2));
    ASSERT_FALSE(KstateTrait::view_compare_equality(view_1, view_3));
    ASSERT_FALSE(KstateTrait::view_compare_equality(view_1, view_4));
    ASSERT_FALSE(KstateTrait::view_compare_less(view_1, view_1));
    ASSERT_TRUE(KstateTrait::view_compare_less(view_1, view_2));
    ASSERT_FALSE(KstateTrait::view_compare_less(view_2, view_1));
    ASSERT_TRUE(KstateTrait::view_compare_less(view_1, view_3));
    ASSERT_FALSE(KstateTrait::view_compare_less(view_3, view_1));
    ASSERT_TRUE(KstateTrait::view_compare_less(view_1, view_4));
    ASSERT_FALSE(KstateTrait::view_compare_less(view_4, view_1));
}

TEST(DynamicIntegralKstateTrait, Test3) {
    using Kstate = kstate_impl::DynamicIntegral64Kstate<MySiteStateTrait, 3u>;
    using KstateTrait = kstate_trait::TraitKstate<Kstate>;
    const std::array<MySiteState, 6> array_1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const auto kstate_1 = KstateTrait::from_range(array_1);
    ASSERT_TRUE(boost::equal(kstate_1.to_range(), array_1));
    const auto view_1 = KstateTrait::to_view(kstate_1);
    ASSERT_TRUE(KstateTrait::view_compare_equality(view_1, view_1));
    const auto kstate_1_re = KstateTrait::from_view(view_1);
    ASSERT_TRUE(boost::equal(kstate_1_re.to_range(), array_1));
    const auto view_1_re = KstateTrait::to_view(kstate_1_re);
    ASSERT_TRUE(KstateTrait::view_compare_equality(view_1, view_1_re));
}
