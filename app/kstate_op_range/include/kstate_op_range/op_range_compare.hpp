#pragma once

#include <kstate_op_range/op_range_raw_adaptors.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/algorithm/search.hpp>
#include <boost/range.hpp>

#include <iterator>
#include <optional>
#include <cmath>

// #######################################################################
// ## compare_XXX                                                       ##
// #######################################################################

namespace kstate_op_range {

template <typename ForwardRange1, typename ForwardRange2>
bool compare_less(const ForwardRange1& rng1, const ForwardRange2& rng2) noexcept {
    assert(boost::size(rng1) == boost::size(rng2));
    return boost::lexicographical_compare(rng1, rng2);
}

template <typename ForwardRange1, typename ForwardRange2>
bool compare_equality(const ForwardRange1& rng1, const ForwardRange2& rng2) noexcept {
    assert(boost::size(rng1) == boost::size(rng2));
    return boost::range::equal(rng1, rng2);
}

template <typename ForwardRange1, typename ForwardRange2>
std::optional<size_t> compare_translational_equality(const ForwardRange1& rng1, const ForwardRange2& rng2) noexcept {
    using raw::adaptors::operator|;
    assert(boost::size(rng1) == boost::size(rng2));
    const auto rng2d = rng2 | kstate_view_amend_spec::doubled;
    const auto it = boost::range::search(rng2d, rng1);
    return it == std::end(rng2d)
            ? std::optional<size_t>()
            : static_cast<size_t>(std::distance(std::begin(rng2d), it));
}

}  // namespace kstate_op_range
