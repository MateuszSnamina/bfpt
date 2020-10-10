#pragma once

#include<functional>
#include<set>
#include<utility>

// #######################################################################
// ## is_minimum                                                        ##
// #######################################################################

namespace monostar_hamiltonians::function_analyzer {

double is_minimum(const std::function<double(double)>& fn, double x);

} // end of namespace monostar_hamiltonians

// #######################################################################
// ## find_all_extrema                                                  ##
// #######################################################################

/*
 * Finds all x₀ such as
 * `min(a,b) ≤ x₀ ≤ max(a,b)` ⋀ `x₀ - extrememum of f(x)`.
 *
 * The function uses the provided `[d/dx](fn)` function.
 *
 * Args:
 * `fn_prim` for `[d/dx](fn)` function.
 * `range` for `{a, b}`,
 * `n_subranges` is an implemenation parameter.
 *
 */

namespace monostar_hamiltonians::function_analyzer {

std::set<double> find_all_extrema(const std::function<double(double)>& fn_prim,
                                  const std::pair<double, double>& range,
                                  unsigned n_subranges = 100);

} // end of namespace monostar_hamiltonians::function_analyzer

// #######################################################################
// ## find_all_{local,global}_minima                                    ##
// #######################################################################

namespace monostar_hamiltonians::function_analyzer {

std::set<double> find_all_local_minima(
        const std::function<double(double)>& fn,
        const std::function<double(double)>& fn_prim, // d/dx(fn)
        const std::pair<double, double>& range,
        unsigned n_subranges = 100);
std::set<double> find_all_global_minima(
        const std::function<double(double)>& fn,
        const std::function<double(double)>& fn_prim, // d/dx(fn)
        const std::pair<double, double>& range,
        unsigned n_subranges = 10000);
std::set<double> find_all_global_minima_periodic_2_pi(
        const std::function<double(double)>& fn,
        const std::function<double(double)>& fn_prim, // d/dx(fn)
        unsigned n_subranges = 10000);

} // end of namespace monostar_hamiltonians::function_analyzer
