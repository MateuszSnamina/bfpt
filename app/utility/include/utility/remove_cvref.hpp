#ifndef KSTATE_REMOVE_CVREF_HPP
#define KSTATE_REMOVE_CVREF_HPP

#include <type_traits>

// #######################################################################
// ##  remove_cvref (is added in to C++20)                              ##
// #######################################################################

namespace utility {

template <class T>
struct remove_cvref {
    typedef std::remove_cv_t<std::remove_reference_t<T>> type;
};

template <class T>
using remove_cvref_t = typename remove_cvref<T>::type;

}  // namespace utility

#endif  // REMOVE_CVREF_HPP
