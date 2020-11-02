#include <kstate_impl_tests/site_state_trait_for_my_site.h>

#include <kstate_impl/kstate_concrete_stl.hpp>

#include <array>
#include <iterator>

#include <gtest/gtest.h>

TEST(DynamicStlKstateTrait, Test0) {
    using Kstate = kstate_impl::DynamicStlKstate<MySiteStateTrait>;
    using KstateTrait = kstate_trait::TraitKstate<Kstate>;
    const std::array<MySiteState, 6> array_1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // before unique shift: {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // after unique shift:  {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    const auto kstate_1 = KstateTrait::from_range(array_1);
    ASSERT_TRUE(boost::equal(kstate_1.to_range(), array_1));
    const auto view_1 = KstateTrait::to_view(kstate_1);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 0), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 1), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 2), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 3), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 4), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 5), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1), 0u);
}

TEST(DynamicStlKstateTrait, Test1) {
    using Kstate = kstate_impl::DynamicStlKstate<MySiteStateTrait>;
    using KstateTrait = kstate_trait::TraitKstate<Kstate>;
    const std::array<MySiteState, 6> array_1 = {MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(11), MySiteState(14), MySiteState(12)};
    // before unique shift: {MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(11), MySiteState(14), MySiteState(12)};
    // after unique shift:  {MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(12), MySiteState(14), MySiteState(11)};
    const auto kstate_1 = KstateTrait::from_range(array_1);
    ASSERT_TRUE(boost::equal(kstate_1.to_range(), array_1));
    const auto view_1 = KstateTrait::to_view(kstate_1);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 0), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 1), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 2), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 3), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 4), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 5), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1), 4u);
}

TEST(DynamicStlKstateTrait, Test2) {
    using Kstate = kstate_impl::DynamicStlKstate<MySiteStateTrait>;
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

    ASSERT_FALSE(KstateTrait::view_compare_less(view_1, view_4));
    ASSERT_TRUE(KstateTrait::view_compare_less(view_4, view_1));
}

TEST(DynamicStlKstateTrait, Test3) {
    using Kstate = kstate_impl::DynamicStlKstate<MySiteStateTrait>;
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

TEST(DynamicStlKstateTrait, Test4) {
    using Kstate = kstate_impl::DynamicStlKstate<MySiteStateTrait>;
    using KstateTrait = kstate_trait::TraitKstate<Kstate>;
    const std::array<MySiteState, 6> array_1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // before unique shift: {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // after unique shift:  {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // after rotation(2):   {MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12), MySiteState(14), MySiteState(12)};
    const auto kstate_1 = KstateTrait::from_range(array_1);
    ASSERT_TRUE(boost::equal(kstate_1.to_range(), array_1));
    const auto view_1 = KstateTrait::to_view(kstate_1);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 0), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 1), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 2), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 3), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 4), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 5), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1), 0u);
    const kstate_view_amend_spec::RotateHolder rotation_holder{2};
    const auto view_1_rotated = KstateTrait::rotated_view(view_1, rotation_holder);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_rotated, 0), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_rotated, 1), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_rotated, 2), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_rotated, 3), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_rotated, 4), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_rotated, 5), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1_rotated), 4u);
}

TEST(DynamicStlKstateTrait, Test5) {
    using Kstate = kstate_impl::DynamicStlKstate<MySiteStateTrait>;
    using KstateTrait = kstate_trait::TraitKstate<Kstate>;
    const std::array<MySiteState, 6> array_1 = {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // before unique shift:         {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12)};
    // after unique shift:          {MySiteState(12), MySiteState(11), MySiteState(11), MySiteState(12), MySiteState(12), MySiteState(14)};
    // after refined(3, State{15}): {MySiteState(14), MySiteState(12), MySiteState(11), MySiteState(15), MySiteState(12), MySiteState(12)};
    const auto kstate_1 = KstateTrait::from_range(array_1);
    ASSERT_TRUE(boost::equal(kstate_1.to_range(), array_1));
    const auto view_1 = KstateTrait::to_view(kstate_1);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 0), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 1), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 2), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 3), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 4), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1, 5), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1), 0u);
    const kstate_view_amend_spec::RefinedHolder<MySiteState> refined_holder{3, MySiteState(15)};  // MySiteState(12) => idx = 5;
    const auto view_1_refined = KstateTrait::refined_view(view_1, refined_holder);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined, 0), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined, 1), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined, 2), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined, 3), MySiteState(15));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined, 4), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined, 5), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1_refined), 3u);
    const kstate_view_amend_spec::RotateHolder rotated_holder{4};
    const auto view_1_refined_rotated = KstateTrait::rotated_view(view_1_refined, rotated_holder);
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined_rotated, 0), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined_rotated, 1), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined_rotated, 2), MySiteState(14));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined_rotated, 3), MySiteState(12));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined_rotated, 4), MySiteState(11));
    ASSERT_EQ(KstateTrait::view_n_th_site_state(view_1_refined_rotated, 5), MySiteState(15));
    ASSERT_EQ(KstateTrait::view_n_unique_shift(view_1_refined_rotated), 5u);
}
