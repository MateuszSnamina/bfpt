#pragma once

#define KSTATE_IS_BASE_OF_TEMPLATE_HPP

#include <type_traits>

// #######################################################################
// ## is_base_of_template                                               ##
// #######################################################################

// Ref:
// - https://stackoverflow.com/questions/51910808/detect-is-base-of-with-base-class-template for discussion.
// - https://stackoverflow.com/questions/34672441/stdis-base-of-for-template-classes

namespace utility {

template <template <typename...> class C, typename...Ts>
std::true_type is_base_of_template_impl(const C<Ts...>*);

template <template <typename...> class C>
std::false_type is_base_of_template_impl(...);

template <typename T, template <typename...> class C>
using is_base_of_template = decltype(is_base_of_template_impl<C>(std::declval<T*>()));

template <typename T, template <typename...> class C>
inline constexpr bool is_base_of_template_v = is_base_of_template<T, C>::value;

}
