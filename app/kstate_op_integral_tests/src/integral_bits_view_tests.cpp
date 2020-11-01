#include <kstate_op_integral/integral_bits_view.hpp>
#include <kstate_op_integral/integral_bits_buffer.hpp>

#include <kstate_op_integral_tests/site_state_trait_for_my_site.h>

#include <boost/range.hpp>

#include <gtest/gtest.h>

#include <array>
#include <cstdint> //for types like: uint64_t

TEST(IntegralBitsRotatedView, Test0) {
    using IntegralBitsBufferT = kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t>;
    const IntegralBitsBufferT b{0b110010111, 9};
    const kstate_view_amend_spec::RotateHolder h{0};
    const kstate_op_integral::IntegralBitsRotatedView<IntegralBitsBufferT, 3u> rotated_view{b, h};
    EXPECT_EQ(rotated_view.get_n_all_bits(), 9);
    EXPECT_EQ(rotated_view.get_number(), 0b110010111);
}

TEST(IntegralBitsRotatedView, Test1) {
    using IntegralBitsBufferT = kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t>;
    const IntegralBitsBufferT b{0b110010111, 9};
    const kstate_view_amend_spec::RotateHolder h{1};
    const kstate_op_integral::IntegralBitsRotatedView<IntegralBitsBufferT, 3u> rotated_view{b, h};
    EXPECT_EQ(rotated_view.get_n_all_bits(), 9);
    EXPECT_EQ(rotated_view.get_number(), 0b111110010);
}

TEST(IntegralBitsRotatedView, Test2) {
    using IntegralBitsBufferT = kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t>;
    const IntegralBitsBufferT buffer{0b110010111, 9};
    const kstate_view_amend_spec::RotateHolder holder_1{1};
    const kstate_view_amend_spec::RotateHolder holder_2{1};
    using ViewT1 = kstate_op_integral::IntegralBitsRotatedView<IntegralBitsBufferT, 3u>;
    using ViewT2 = kstate_op_integral::IntegralBitsRotatedView<ViewT1, 3u>;
    const ViewT1 rotated_view_1{buffer, holder_1};
    const ViewT2 rotated_view_2{rotated_view_1, holder_2};
    EXPECT_EQ(rotated_view_2.get_n_all_bits(), 9);
    EXPECT_EQ(rotated_view_2.get_number(), 0b010111110);
}

TEST(IntegralBitsRefinedView, Test1) {
    using IntegralBitsBufferT = kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t>;
    const IntegralBitsBufferT b{0b110010111, 9};
    const kstate_view_amend_spec::RefinedHolder<MySiteState> h{0, MySiteState(15)};// MySiteState(12) => idx = 5;
    const kstate_op_integral::IntegralBitsRefinedView<IntegralBitsBufferT, MySiteStateTrait, 3u> refined_view{b, h};
    EXPECT_EQ(refined_view.get_n_all_bits(), 9);
    EXPECT_EQ(refined_view.get_number(), 0b110010101);
}

TEST(IntegralBitsRefinedView, Test2) {
    using IntegralBitsBufferT = kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t>;
    const IntegralBitsBufferT b{0b110010111, 9};
    const kstate_view_amend_spec::RefinedHolder<MySiteState> h{1, MySiteState(15)};// MySiteState(12) => idx = 5;
    const kstate_op_integral::IntegralBitsRefinedView<IntegralBitsBufferT, MySiteStateTrait, 3u> refined_view{b, h};
    EXPECT_EQ(refined_view.get_n_all_bits(), 9);
    EXPECT_EQ(refined_view.get_number(), 0b110101111);
}

TEST(IntegralBitsRefinedView, Test3) {
    using IntegralBitsBufferT = kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t>;
    const IntegralBitsBufferT b{0b110010111, 9};
    const kstate_view_amend_spec::RefinedHolder<MySiteState> h{2, MySiteState(15)};// MySiteState(12) => idx = 5;
    const kstate_op_integral::IntegralBitsRefinedView<IntegralBitsBufferT, MySiteStateTrait, 3u> refined_view{b, h};
    EXPECT_EQ(refined_view.get_n_all_bits(), 9);
    EXPECT_EQ(refined_view.get_number(), 0b101010111);
}
