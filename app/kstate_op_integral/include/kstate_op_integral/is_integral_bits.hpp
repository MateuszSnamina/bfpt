#pragma once

#include <type_traits>

namespace kstate_op_integral {

template <typename T>
struct IsIntegralBits : std::false_type {};

}  // namespace kstate_op_integral
