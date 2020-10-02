#include<model_monostar/simple_numerical_function_analyzer.hpp>

#include<model_monostar/simple_numerical_root_finder.hpp>

#include<utility/almost_equal.hpp>

#include<limits>

#include<cmath>
#include<cassert>

#include<iostream> //TODO remove (debug)

// #######################################################################
// ## M_PI                                                              ##
// #######################################################################

// Ref: https://stackoverflow.com/questions/1727881/how-to-use-the-pi-constant-in-c
#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

// #######################################################################
// ## is_minimum                                                        ##
// #######################################################################

double is_minimum(const std::function<double(double)>& fn, double x) {
    const double fnx = fn(x);
    for(double epsilon = 10 * std::numeric_limits<double>::epsilon(); epsilon < std::numeric_limits<double>::max() / 1000; epsilon *=2) {
        const double a = x + epsilon;
        const double b = x - epsilon;
        const double fna = fn(a);
        const double fnb = fn(b);
        if ((fna > fnx && !utility::almost_equal(fna, fnx, 1.0)) &&
                (fnb > fnx && !utility::almost_equal(fnb, fnx, 1.0))) {
            return true;
        } else if ((fna < fnx && !utility::almost_equal(fna, fnx, 1.0)) ||
                   (fnb < fnx && !utility::almost_equal(fnb, fnx, 1.0))) {
            return false;
        } else {
            continue;
        }
    }
    return false;
}

// #######################################################################
// ## find_all_extrema                                                  ##
// #######################################################################

std::set<double> find_all_extrema(const std::function<double(double)>& fn_prim,
                                  const std::pair<double, double>& range,
                                  unsigned n_subranges) {
    const std::set<double> all_extrema_raw = find_zero_in_subranges(fn_prim, range, n_subranges);
    const std::set<double> all_extrema_refined = utility::remove_almost_equal_numbers(all_extrema_raw, (range.second - range.first) / n_subranges);
    return all_extrema_refined;
}

// #######################################################################
// ## find_all_{local,global}_minima                                    ##
// #######################################################################

std::set<double> find_all_local_minima(
        const std::function<double(double)>& fn,
        const std::function<double(double)>& fn_prim, // d/dx(fn)
        const std::pair<double, double>& range,
        unsigned n_subranges) {
    const std::set<double> all_extrema = find_all_extrema(fn_prim, range, n_subranges);
    std::set<double> all_local_minima;
    for (const auto& extremum : all_extrema) {
        if  (is_minimum(fn, extremum)) {
            all_local_minima.insert(extremum);
        }
    }
    return all_local_minima;
}

std::set<double> find_all_global_minima(
        const std::function<double(double)>& fn,
        const std::function<double(double)>& fn_prim, // d/dx(fn)
        const std::pair<double, double>& range,
        unsigned n_subranges) {
    const std::set<double> all_local_minima = find_all_local_minima(fn, fn_prim, range, n_subranges);
    if (all_local_minima.empty()) {
        return {};
    } else {
        std::set<double> results;
        std::set<double> all_fn_local_minima = [&all_local_minima, &fn](){
            std::set<double> all_fn_local_minima_builder;
            for (const auto local_minimum : all_local_minima) {
                all_fn_local_minima_builder.insert(fn(local_minimum));
            }
            return all_fn_local_minima_builder;
        }();
        double fn_global_minimum = *all_fn_local_minima.cbegin();
        for (const auto local_minimum : all_local_minima) {
            double fn_local_minimum = fn(local_minimum);
            if (utility::almost_equal(fn_global_minimum, fn_local_minimum, 1.0)) {
                results.insert(local_minimum);
            }
        }
        return results;
    }
}

std::set<double> find_all_global_minima_periodic_2_pi(
        const std::function<double(double)>& fn,
        const std::function<double(double)>& fn_prim, // d/dx(fn)
        unsigned n_subranges) {
    const std::set<double> all_local_minima = [&fn, &fn_prim, &n_subranges]() {
        std::set<double> all_local_minima_builder = find_all_local_minima(fn, fn_prim, {-M_PI, +M_PI}, n_subranges);
        if (is_minimum(fn, -M_PI)) {
            all_local_minima_builder.insert(-M_PI);
        }
        if (is_minimum(fn, +M_PI)) {
            all_local_minima_builder.insert(+M_PI);
        }
        return all_local_minima_builder;
    }();

    if (all_local_minima.empty()) {
        return {};
    } else {
        std::set<double> results;
        std::set<double> all_fn_local_minima = [&all_local_minima, &fn](){
            std::set<double> all_fn_local_minima_builder;
            for (const auto local_minimum : all_local_minima) {
                all_fn_local_minima_builder.insert(fn(local_minimum));
            }
            return all_fn_local_minima_builder;
        }();
        double fn_global_minimum = *all_fn_local_minima.cbegin();
        for (const auto local_minimum : all_local_minima) {
            double fn_local_minimum = fn(local_minimum);
            if (utility::almost_equal(fn_global_minimum, fn_local_minimum, 1.0)) {
                results.insert(local_minimum);
            }
        }
        if (utility::almost_equal(*results.cbegin(), -M_PI) && utility::almost_equal(*results.crbegin(), +M_PI)) {
            results.erase(*results.cbegin());
        }
        return results;
    }
}
