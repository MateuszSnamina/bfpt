#pragma once

#include<kstate_op_integral/integral_bits.hpp>
#include<kstate_op_integral/op_integral_bits.hpp>

#include <cmath>
#include <cassert>

/*

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
size_t n_least_replication_shift(IntegralBitsT integral_bits) {
    assert(boost::size(integral_bits) > 0);
    const auto rngdr = integral_bits | extension::boost::adaptors::doubled | extension::boost::adaptors::rotated(1);
    const auto it = boost::range::search(rngdr, integral_bits);
    const auto _ = std::distance(std::begin(rngdr), it);
    assert(_ >= 0);
    return static_cast<size_t>(_ + 1);
}

}  // namespace kstate_op_integral

// #######################################################################
// ## norm_factor                                                       ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
double norm_factor(IntegralBitsT integral_bits) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const size_t n_sites = boost::size(integral_bits);
    return std::sqrt(n_least_replication_shift(integral_bits)) / n_sites;
    // The result is equal to 1 / std::sqrt(n_least_replication_shift) / n_replicas;
    // where n_replicas = n_sites / n_least_replication_shift
}

}  // namespace kstate_op_integral

// #######################################################################
// ## is_prolific                                                       ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
bool is_prolific(IntegralBitsT integral_bits, int n_k) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const size_t n_sites = boost::size(integral_bits);
    return !((n_least_replication_shift(integral_bits) * n_k) % n_sites);
}

}  // namespace kstate_op_integral
*/
