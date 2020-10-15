#include<kstate_op_integral/op_integral_unique_shift.hpp>
//#include<kstate_op_integral/op_integral_compare.hpp>

#include <gtest/gtest.h>

#include <cstdint> //for types like: uint64_t

TEST(KstateOpIntegral, NUniqueShiftTest0) {
    const kstate_op_integral::IntegralBitsDynamic<uint64_t> b1{0b1, 1};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1), 0);
}

TEST(KstateOpIntegral, NUniqueShiftTest1) {
    const kstate_op_integral::IntegralBitsDynamic<uint64_t> b1{0b01, 2};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1), 1);
}

TEST(KstateOpIntegral, NUniqueShiftTest2) {
    const kstate_op_integral::IntegralBitsDynamic<uint64_t> b1{0b10, 2};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1), 0);
}

TEST(KstateOpIntegral, NUniqueShiftTest3) {
    const kstate_op_integral::IntegralBitsDynamic<uint64_t> b1{0b1101101001, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1), 1);
}

TEST(KstateOpIntegral, NUniqueShiftTest4) {
    const kstate_op_integral::IntegralBitsDynamic<uint64_t> b1{0b1101101000, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1), 0);
}

TEST(KstateOpIntegral, NUniqueShiftTest5) {
    const kstate_op_integral::IntegralBitsDynamic<uint64_t> b1{0b0001011011, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1), 5);
}

TEST(KstateOpIntegral, NUniqueShiftTest6) {
    const kstate_op_integral::IntegralBitsDynamic<uint64_t> b1{0b1100001100, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1), 4);
}

TEST(KstateOpIntegral, NUniqueShiftTest7) {
    const kstate_op_integral::IntegralBitsDynamic<uint64_t> b1{0b0110100110, 10};
    ASSERT_EQ(kstate_op_integral::n_unique_shift(b1), 9);
}

/*
TEST(KstateOpIntegral, MakeUniqueShiftTest0) {
    const std::array<int, 1> v2 = {11};
    const std::array<int, 1> v2u = {11};
    kstate_op_integral::compare_equality(v2u, kstate_op_integral::make_unique_shift(v2));
}

TEST(KstateOpIntegral, MakeUniqueShiftTest1) {
    const std::array<int, 2> v2 = {11, 12};
    const std::array<int, 2> v2u = {12, 11};
    kstate_op_integral::compare_equality(v2u, kstate_op_integral::make_unique_shift(v2));
}

TEST(KstateOpIntegral, MakeUniqueShiftTest2) {
    const std::array<int, 2> v2 = {12, 11};
    const std::array<int, 2> v2u = {12, 11};
    kstate_op_integral::compare_equality(v2u, kstate_op_integral::make_unique_shift(v2));
}

TEST(KstateOpIntegral, MakeUniqueShiftTest3) {
    const std::array<int, 7> v2 = {12, 11, 14, 13, 14, 14, 13};
    const std::array<int, 7> v2u = {14, 14, 13, 12, 11, 14, 13};
    kstate_op_integral::compare_equality(v2u, kstate_op_integral::make_unique_shift(v2));
}
*/
