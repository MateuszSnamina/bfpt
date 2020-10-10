#pragma once

#include<functional>
#include<set>
#include<optional>
#include<utility>

// #######################################################################
// ## find_zero_in_given_range                                          ##
// #######################################################################

/*
 * Finds x₀ such as
 * `min(a,b) ≤ x₀ ≤ max(a,b)` ⋀ `f(x₀) = 0`.
 *
 * The function works if (f(a)>0 ⋀ f(b)<0) ⋁ (f(a)<0 ⋀ f(b)>0)
 * and assumes there is exactly one solution of the problem.
 */

namespace monostar_hamiltonians::root_finder {

std::optional<double> find_zero_in_given_range(
        const std::function<double(double)>& fn,
        std::pair<double, double> range);

} // end of namespace monostar_hamiltonians::root_finder

// #######################################################################
// ## find_zero_in_subranges                                            ##
// #######################################################################

/*
 * Repeats the find_zero_in_given_range function
 * for n_subranges subranges of the given range.
 */

namespace monostar_hamiltonians::root_finder {

std::set<double> find_zero_in_subranges(
        const std::function<double(double)>& fn,
        const std::pair<double, double>& range,
        unsigned n_subranges = 100);

} // end of namespace monostar_hamiltonians::root_finder
