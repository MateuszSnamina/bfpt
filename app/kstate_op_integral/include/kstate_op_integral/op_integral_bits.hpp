#pragma once

#include <boost/range.hpp>
#include <boost/range/any_range.hpp>

//#include <cstdint> for types like: uint64_t
#include <iterator>
#include <type_traits>
#include <memory>
#include <cassert>

// #######################################################################
// ## init_vector_from_range (a helper function)                        ##
// #######################################################################

namespace kstate_op_integral {

template <typename T>
bool extract_bit(T n, unsigned idx) {
    static_assert(std::is_integral_v<T>);
    assert(idx < 8 * sizeof(T));
    return static_cast<bool>(n & (static_cast<T>(1u) << idx));
}

} // end of namespace kstate_op_integral

