#include<monostar_hamiltonians/simple_numerical_root_finder.hpp>

#include<utility/almost_equal.hpp>

#include<cmath>
#include<cassert>

// #######################################################################
// ## find_zero_in_given_range                                          ##
// #######################################################################

namespace monostar_hamiltonians::root_finder {

std::optional<double> find_zero_in_given_range(
        const std::function<double(double)>& fn,
        std::pair<double, double> range) {
    double& a = range.first;
    double& b = range.second;
    const double problem_scale = std::abs(b - a);
    double fna = fn(a);
    double fnb = fn(b);
    if (fna == 0) {
        return a;
    } else if (fnb == 0) {
        return b;
    } else if (std::copysign(1.0, fna) == std::copysign(1.0, fnb)) {
        return std::nullopt;
    } else if ((std::copysign(1.0, fna) == +1 && std::copysign(1.0, fnb) == -1) ||
               (std::copysign(1.0, fna) == -1 && std::copysign(1.0, fnb) == +1) ){
        assert(fna != 0);
        assert(fnb != 0);
        while (!utility::almost_equal(a,b, problem_scale)) {
            assert(std::copysign(1.0, fna) == -1 || std::copysign(1.0, fna) == +1);
            assert(std::copysign(1.0, fnb) == -1 || std::copysign(1.0, fnb) == +1);
            assert(std::copysign(1.0, fna) != std::copysign(1.0, fnb));
            const double c = (b + a) / 2;
            const double fnc = fn(c);
            if (fnc == 0) {
                a = c;
                b = c;
            } else if (std::copysign(1.0, fnc) == std::copysign(1.0, fnb)) {
                assert((std::copysign(1.0, fna) == +1 && std::copysign(1.0, fnc) == -1) ||
                       (std::copysign(1.0, fna) == -1 && std::copysign(1.0, fnc) == +1));
                b = c;
                fnb = fnc;
            } else if (std::copysign(1.0, fna) == std::copysign(1.0, fnc)) {
                assert((std::copysign(1.0, fnc) == +1 && std::copysign(1.0, fnb) == -1) ||
                       (std::copysign(1.0, fnc) == -1 && std::copysign(1.0, fnb) == +1));
                a = c;
                fna = fnc;
            } else {
                // fn(c) is NaN.
                return std::nullopt;
            }
        } // end of while loop
        return (a + b) / 2;
    } else {
        // at least one of: fn(a), fn(b) is NaN.
        return std::nullopt;
    }
}

} // end of namespace monostar_hamiltonians::root_finder

// #######################################################################
// ## find_zero_in_subranges                                            ##
// #######################################################################

namespace monostar_hamiltonians::root_finder {

std::set<double> find_zero_in_subranges(
        const std::function<double(double)>& fn,
        const std::pair<double, double>& range,
        unsigned n_subranges) {
    assert(n_subranges > 0);
    const std::optional<std::pair<double, double>> sanitized_range = [&range]() -> std::optional<std::pair<double, double>> {
        const double& a = range.first;
        const double& b = range.second;
        if (a < b) {
            return std::make_pair(a, b);
        } else if (a < b) {
            return std::make_pair(b, a);
        } else {
            return std::nullopt;
        }
    }();
    if (!sanitized_range) {
        return {};
    }
    const double& a = (*sanitized_range).first;
    const double& b = (*sanitized_range).second;
    const double problem_scale = (b - a);
    std::set<double> result;
    for (unsigned n_subrange = 0; n_subrange < n_subranges; n_subrange++) {
        const double subrange_size = (b - a) / n_subranges;
        const std::pair<double, double> subrange = {a + n_subrange * subrange_size, a + (n_subrange + 1) * subrange_size};
        if (const auto result_subrange = find_zero_in_given_range(fn, subrange)) {
            result.insert(*result_subrange);
        }
    }
    result = utility::remove_almost_equal_numbers(result, problem_scale);
    return result;
}

} // end of namespace monostar_hamiltonians::root_finder
