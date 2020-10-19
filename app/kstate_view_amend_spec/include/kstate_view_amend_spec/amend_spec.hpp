#pragma once

#include<cstddef>

// #######################################################################
// ##  rotated                                                          ##
// #######################################################################

namespace kstate_view_amend_spec {

class RotateHolder {
   public:
    RotateHolder(size_t n) : n(n){};
    const size_t n;
};

inline RotateHolder rotated(size_t n) {
    return RotateHolder(n);
}

} // namespace kstate_view_amend_spec

// #######################################################################
// ##  doubled                                                          ##
// #######################################################################

namespace kstate_view_amend_spec {

class Doubler {};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
static Doubler doubled{};
#pragma GCC diagnostic pop

} // namespace kstate_view_amend_spec

// #######################################################################
// ##  refined                                                          ##
// #######################################################################

namespace kstate_view_amend_spec {

template <typename T>
struct RefinedHolder {
    RefinedHolder(size_t n, T v)
        : n(n),
          v{v} {};
    const size_t n;
    const T v[1];
};

template <typename T>
RefinedHolder<T> refined(size_t n, T v) {
    return RefinedHolder<T>(n, v);
}

} // namespace kstate_view_amend_spec
