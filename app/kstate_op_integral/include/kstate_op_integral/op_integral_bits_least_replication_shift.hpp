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
size_t n_least_replication_shift(IntegralBitsT integral_bits) {
    // LETS ASSUME A SITE PER BIT. HOLDS FOR TWO-LEVEL SYSTEM ONLY!!
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    const auto n_buffer = integral_bits.get_number();
    auto n2_buffer = integral_bits.get_number();
    for (size_t _ = 1; _ < n_all_bits; _++) {
        n2_buffer = ::kstate_op_integral::raw::rotate(n2_buffer, n_all_bits, 1);
        if (n_buffer == n2_buffer) {
            return _;
        }
    }
    return n_all_bits;
}

}  // namespace kstate_op_integral

// #######################################################################
// ## norm_factor                                                       ##
// #######################################################################

namespace kstate_op_integral {

template <typename IntegralBitsT>
double norm_factor(IntegralBitsT integral_bits) noexcept {
    static_assert(IsIntegralBits<IntegralBitsT>::value);
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    const auto n_sites = n_all_bits; // LETS ASSUME A SITE PER BIT. HOLDS FOR TWO-LEVEL SYSTEM ONLY!!
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
    const auto n_all_bits = integral_bits.get_n_all_bits();
    assert(n_all_bits < 8 * sizeof(typename IntegralBitsT::IntegralT));
    const auto n_sites = n_all_bits; // LETS ASSUME A SITE PER BIT. HOLDS FOR TWO-LEVEL SYSTEM ONLY!!
    return !((n_least_replication_shift(integral_bits) * n_k) % n_sites);
}

}  // namespace kstate_op_integral
