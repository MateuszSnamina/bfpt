#include <kstate_op_integral/op_integral_bits_unique_shift.hpp>
#include <kstate_op_integral/integral_bits_buffer.hpp>

#include <gtest/gtest.h>

#include <cstdint>  //for types like: uint64_t

TEST(KstateOpIntegral, NUniqueShiftTest0) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1, 1};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 0);
}

TEST(KstateOpIntegral, NUniqueShiftTest1) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b01, 2};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 1);
}

TEST(KstateOpIntegral, NUniqueShiftTest2) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b10, 2};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 0);
}

TEST(KstateOpIntegral, NUniqueShiftTest3) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1101101001, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 1);
}

TEST(KstateOpIntegral, NUniqueShiftTest4) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1101101000, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 0);
}

TEST(KstateOpIntegral, NUniqueShiftTest5) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b0001011011, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 5);
}

TEST(KstateOpIntegral, NUniqueShiftTest6) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1100001100, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 4);
}

TEST(KstateOpIntegral, NUniqueShiftTest7) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b0110100110, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 9);
}

TEST(KstateOpIntegral, NUniqueShiftTest8) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b011110, 6};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 1), 5);
}

TEST(KstateOpIntegral, NUniqueShiftTest9) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b011110, 6};  // 1=b01 3=b11 2=b10
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1, 2), 2);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest0) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1, 1};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b1);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest1) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b01, 2};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b10);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest2) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b10, 2};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b10);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest3) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1101101001, 10};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b1110110100);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest4) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1101101000, 10};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b1101101000);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest5) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b0001011011, 10};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b1101100010);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest6) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b1100001100, 10};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b1100110000);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest7) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b0110100110, 10};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b1101001100);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest8) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b011110, 6};
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 1), 0b111100);
}

TEST(KstateOpIntegral, MakeUniqueShiftTest9) {
    const kstate_op_integral::IntegralBitsDynamicBuffer<uint64_t> b1{0b011110, 6};  // 1=b01 3=b11 2=b10
    ASSERT_EQ(kstate_op_integral::make_unique_shift(b1, 2), 0b111001);
}
