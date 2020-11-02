#pragma once

#include <kstate_op_range/op_range_raw_adaptors.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/algorithm/search.hpp>

#include <iterator>

#include <cmath>
#include <cassert>

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_range {

template <typename ForwardRange>
size_t n_least_replication_shift(const ForwardRange& rng) {
    // OLD IMPLEMENTATION:
    //    using raw::adaptors::operator|;
    //    assert(boost::size(rng) > 0);
    //    const auto rngdr = rng | kstate_view_amend_spec::doubled | kstate_view_amend_spec::rotated(1);
    //    const auto it = boost::range::search(rngdr, rng);
    //    const auto _ = std::distance(std::begin(rngdr), it);
    //    assert(_ >= 0);
    //    return static_cast<size_t>(_ + 1);
    // NEW IMPLEMENTATION:
    using raw::adaptors::operator|;
    [[maybe_unused]] const unsigned rng_size = boost::size(rng);
    assert(rng_size > 0);
    for (unsigned i = 1; i < rng_size; i++) {
        if (boost::range::equal(rng, rng | kstate_view_amend_spec::rotated(i))) {
            return i;
        }
    }
    return rng_size;
}

}  // namespace kstate_op_range

// #######################################################################
// ## norm_factor                                                       ##
// #######################################################################

namespace kstate_op_range {

template <typename ForwardRange>
double norm_factor(const ForwardRange& rng) noexcept {
    const size_t n_sites = boost::size(rng);
    return std::sqrt(n_least_replication_shift(rng)) / n_sites;
    // The result is equal to 1 / std::sqrt(n_least_replication_shift) / n_replicas;
    // where n_replicas = n_sites / n_least_replication_shift
}

}  // namespace kstate_op_range

// #######################################################################
// ## is_prolific                                                       ##
// #######################################################################

namespace kstate_op_range {

template <typename ForwardRange>
bool is_prolific(const ForwardRange& rng, int n_k) noexcept {
    const size_t n_sites = boost::size(rng);
    return !((n_least_replication_shift(rng) * n_k) % n_sites);
}

}  // namespace kstate_op_range
