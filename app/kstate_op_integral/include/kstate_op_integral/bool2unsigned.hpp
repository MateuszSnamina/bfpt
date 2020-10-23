#pragma once

namespace kstate_op_integral {

constexpr bool bool_from_unsiged(unsigned u) {
    return (u ? true : false);
}

constexpr unsigned bool_to_unsigned(bool b) {
    return (b ? 1u : 0u);
}

}
