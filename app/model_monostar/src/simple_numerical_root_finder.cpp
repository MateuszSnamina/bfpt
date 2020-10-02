#include<model_monostar/simple_numerical_root_finder.hpp>

#include<cmath>
#include<cassert>

// #######################################################################
// ## RootFinder                                                        ##
// #######################################################################

std::optional<double> find_zero_in_given_range(
        const std::function<double(double)>& fn,
        std::pair<double, double> range) {
    double& a = range.first;
    double& b = range.second;
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
        while ( std::abs(b-a) > 5 * std::numeric_limits<double>::epsilon() ) {
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
    std::set<double> result;
    for (unsigned n_subrange = 0; n_subrange < n_subranges; n_subrange++) {
        const double& a = (*sanitized_range).first;
        const double& b = (*sanitized_range).second;
        const double subrange_size = (b - a) / n_subranges;
        const std::pair<double, double> subrange = {a + n_subrange * subrange_size, a + (n_subrange + 1) * subrange_size};
        if (const auto result_subrange = find_zero_in_given_range(fn, subrange)){
            result.insert(*result_subrange);
        }
    }
    return result;
}
