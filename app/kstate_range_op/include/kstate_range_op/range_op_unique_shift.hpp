#pragma once

#include <extensions/adaptors.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <iterator>

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_range_op {

template <typename ForwardRange>
size_t n_unique_shift(const ForwardRange& rng) noexcept {
    using Difference = typename boost::range_difference<ForwardRange>::type;
    const Difference d = std::distance(std::begin(rng), std::end(rng));
    size_t i = 0;
    for (size_t _ = 1; boost::numeric_cast<Difference>(_) < d; _++) {
        const bool result_cmp = boost::lexicographical_compare(
            rng | extension::boost::adaptors::rotated(i),
            rng | extension::boost::adaptors::rotated(_));
        if (result_cmp) {
            i = _;
        }
    }
    return i;
}

}  // namespace kstate_range_op

// #######################################################################
// ## make_unique_shift                                                 ##
// #######################################################################

namespace kstate_range_op {

template <typename ForwardRange>
extension::boost::adaptors::RotatedRangeType<ForwardRange> make_unique_shift(const ForwardRange& rng) noexcept {
    return rng | extension::boost::adaptors::rotated(n_unique_shift(rng));
}

}  // namespace kstate_range_op
