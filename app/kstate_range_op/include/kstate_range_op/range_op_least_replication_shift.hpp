#pragma once

#include <extensions/adaptors.hpp>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/algorithm/search.hpp>

#include <iterator>

#include <cmath>
#include <cassert>

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_range_op {

template <typename ForwardRange>
size_t n_least_replication_shift(const ForwardRange& rng) {
    assert(boost::size(rng) > 0);
    const auto rngdr = rng | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    const auto it = boost::range::search(rngdr, rng);
    const auto _ = std::distance(std::begin(rngdr), it);
    assert(_ >= 0);
    return static_cast<size_t>(_ + 1);
}

}  // namespace kstate_range_op

// #######################################################################
// ## norm_factor                                                       ##
// #######################################################################

namespace kstate_range_op {

template <typename ForwardRange>
double norm_factor(const ForwardRange& rng) noexcept {
    const size_t n_sites = boost::size(rng);
    return std::sqrt(n_least_replication_shift(rng)) / n_sites;
    // The result is equal to 1 / std::sqrt(n_least_replication_shift) / n_replicas;
    // where n_replicas = n_sites / n_least_replication_shift
}

}  // namespace kstate_range_op

// #######################################################################
// ## is_prolific                                                       ##
// #######################################################################

namespace kstate_range_op {

template <typename ForwardRange>
bool is_prolific(const ForwardRange& rng, int n_k) noexcept {
    const size_t n_sites = boost::size(rng);
    return !((n_least_replication_shift(rng) * n_k) % n_sites);
}

}  // namespace kstate_range_op
