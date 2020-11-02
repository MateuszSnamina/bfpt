#ifndef UTILITY_ALMOST_EQUAL_HPP
#define UTILITY_ALMOST_EQUAL_HPP

#include <type_traits>
#include <numeric>

#include <cmath>
#include <set>

// #######################################################################
// ## almost_equal                                                      ##
// #######################################################################

namespace utility {

// Ref: https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
almost_equal(T x, T y, T scale = 0.0, int ulp = 1e3) {
    return std::fabs(x - y) <= std::numeric_limits<T>::epsilon() * (std::fabs(scale) + std::fabs(x) + std::fabs(y)) * ulp || std::fabs(x - y) < std::numeric_limits<T>::min();
}

template <class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, std::set<T>>::type
remove_almost_equal_numbers(const std::set<double>& numbers, T scale = 0.0, int ulp = 1e3) {
    std::set<T> result;
    for (const auto& number : numbers) {
        if (result.empty() || !utility::almost_equal(*(result.crbegin()), number, scale, ulp)) {
            result.insert(number);
        }
    }
    return result;
}

}  // end of namespace utility

#endif  // UTILITY_ALMOST_EQUAL_HPP
