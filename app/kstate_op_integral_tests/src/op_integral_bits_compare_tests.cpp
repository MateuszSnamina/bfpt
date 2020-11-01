#include<kstate_op_integral/integral_bits_buffer.hpp>
#include<kstate_op_integral/op_integral_bits_compare.hpp>

#include <gtest/gtest.h>

#include <cstdint> //for types like: uint64_t

TEST(KstateOpIntegral, CompareEqualityTest0) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1100101110, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b2{0b1100101100, 10};
    EXPECT_TRUE(kstate_op_integral::compare_equality(b1, b1));
    EXPECT_FALSE(kstate_op_integral::compare_equality(b1, b2));
}

TEST(KstateOpIntegral, CompareTranlationalEqualityTest0) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1100101110, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b2{0b1100101100, 10};

    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b11{0b1001011101, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b12{0b0010111011, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b13{0b0101110110, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b14{0b1011101100, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b15{0b0111011001, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b16{0b1110110010, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b17{0b1101100101, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b18{0b1011001011, 10};
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b19{0b0110010111, 10};
    // ---
    ASSERT_FALSE(kstate_op_integral::compare_translational_equality(b1, b2));
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b1));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b1), 0);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b11));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b11), 1);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b12));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b12), 2);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b13));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b13), 3);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b14));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b14), 4);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b15));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b15), 5);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b16));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b16), 6);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b17));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b17), 7);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b18));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b18), 8);
    ASSERT_TRUE(kstate_op_integral::compare_translational_equality(b1, b19));
    ASSERT_EQ(*kstate_op_integral::compare_translational_equality(b1, b19), 9);
}
