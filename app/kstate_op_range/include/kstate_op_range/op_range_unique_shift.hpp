#pragma once

#include <kstate_op_range/op_range_raw_adaptors.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/numeric/conversion/cast.hpp>

#include <iterator>

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_range {

template <typename ForwardRange>
size_t n_unique_shift(const ForwardRange& rng) noexcept {
    using raw::adaptors::operator|;
    using Difference = typename boost::range_difference<ForwardRange>::type;
    const Difference d = std::distance(std::begin(rng), std::end(rng));
    size_t i = 0;
    for (size_t _ = 1; boost::numeric_cast<Difference>(_) < d; _++) {
        const bool result_cmp = boost::lexicographical_compare(
            rng | kstate_view_amend_spec::rotated(i),
            rng | kstate_view_amend_spec::rotated(_));
        if (result_cmp) {
            i = _;
        }
    }
    return i;
}

}  // namespace kstate_op_range

// #######################################################################
// ## make_unique_shift                                                 ##
// #######################################################################

namespace kstate_op_range {

template <typename ForwardRange>
raw::RotatedRangeType<ForwardRange> make_unique_shift(const ForwardRange& rng) noexcept {
    using raw::adaptors::operator|;
    return rng | kstate_view_amend_spec::rotated(n_unique_shift(rng));
}

}  // namespace kstate_op_range
