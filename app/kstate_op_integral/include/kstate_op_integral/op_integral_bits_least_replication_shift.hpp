#pragma once

#include<kstate_op_integral/is_integral_bits.hpp>
#include<kstate_op_integral/op_integral_bits_raw.hpp>

#include <cmath>
#include <cassert>

// #######################################################################
// ## n_unique_shift                                                    ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
size_t n_least_replication_shift(IntegralBitsT integral_bits, unsigned n_bits_per_site) {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    assert(n_all_bits % n_bits_per_site == 0);
    const auto n_buffer = integral_bits.get_number();
    auto n_rotated_buffer = integral_bits.get_number();
    size_t result = 1;
    for (size_t _ = n_bits_per_site; _ < n_all_bits; _ += n_bits_per_site) {
        n_rotated_buffer = ::kstate_op_integral::raw::rotate(n_rotated_buffer, n_all_bits, n_bits_per_site);
        if (n_buffer == n_rotated_buffer) {
            break;
        }
        result++;
    }
    return result;
}

}  // namespace kstate_op_integral

// #######################################################################
// ## norm_factor                                                       ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
double norm_factor(IntegralBitsT integral_bits, unsigned n_bits_per_site) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    assert(n_all_bits % n_bits_per_site == 0);
    const auto n_sites = n_all_bits / n_bits_per_site;
    return std::sqrt(n_least_replication_shift(integral_bits, n_bits_per_site)) / n_sites;
    // The result is equal to 1 / std::sqrt(n_least_replication_shift) / n_replicas;
    // where n_replicas = n_sites / n_least_replication_shift
}

}  // namespace kstate_op_integral

// #######################################################################
// ## is_prolific                                                       ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
bool is_prolific(IntegralBitsT integral_bits, unsigned n_bits_per_site, int n_k) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    assert(n_all_bits % n_bits_per_site == 0);
    const auto n_sites = n_all_bits / n_bits_per_site;
    return !((n_least_replication_shift(integral_bits, n_bits_per_site) * n_k) % n_sites);
}

}  // namespace kstate_op_integral
