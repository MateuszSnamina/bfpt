#pragma once

#include <kstate_trait/trait_site_state.hpp>

#include <kstate_op_integral/is_integral_bits.hpp>
#include <kstate_op_integral/op_integral_bits_raw.hpp>
#include <kstate_op_integral/bool2unsigned.hpp>

#include <kstate_view_amend_spec/amend_spec.hpp>

#include <type_traits>
#include <cassert>


// #######################################################################
// ## IntegralBitsView                                                  ##
// #######################################################################

namespace kstate_op_integral {

template<typename _IntegralBitsT>
struct IntegralBitsView {
    static_assert(IsIntegralBits<_IntegralBitsT>::value);
    using IntegralBitsT = _IntegralBitsT;
    using IntegralT = typename IntegralBitsT::IntegralT;
    const IntegralBitsT& integral_bits;
    const IntegralT& get_number_ref() const noexcept {
        return integral_bits.get_number_ref();
    }
    IntegralT get_number() const noexcept {
        return integral_bits.get_number();
    }
    unsigned char get_n_all_bits() const noexcept {
        return integral_bits.get_n_all_bits();
    };
};

template<typename _IntegralBitsT>
struct IsIntegralBits<IntegralBitsView<_IntegralBitsT>> : std::true_type {};

}

// #######################################################################
// ## IntegralBitsRotatedView                                           ##
// #######################################################################

namespace kstate_op_integral {

template<typename _IntegralBitsT, unsigned char _n_bits_per_site>//TODO use _n_bits_per_site (now we assume it is equal to 1u).
struct IntegralBitsRotatedView {
    static_assert(IsIntegralBits<_IntegralBitsT>::value);
    static_assert(_n_bits_per_site >= 0u);
    using IntegralBitsT = _IntegralBitsT;
    using IntegralT = typename IntegralBitsT::IntegralT;
    static constexpr unsigned char n_bits_per_site = _n_bits_per_site;
    const IntegralBitsT& integral_bits;
    const kstate_view_amend_spec::RotateHolder& h;
    IntegralT get_number() const noexcept {
        const IntegralT underlying_value = integral_bits.get_number();
        const unsigned char n_all_bits = integral_bits.get_n_all_bits();
        const unsigned char idx_pivot_chunk_number = h.n;
        const unsigned char n_bits_in_chunk_number = n_bits_per_site;
        const IntegralT refined_value = kstate_op_integral::raw::rotate_chunk_number(
                    underlying_value,
                    n_all_bits,
                    idx_pivot_chunk_number,
                    n_bits_in_chunk_number);
        return refined_value;
    }
    unsigned char get_n_all_bits() const noexcept {
        return integral_bits.get_n_all_bits();
    };
};

template<typename _IntegralBitsT, unsigned char _n_bits_per_site>
struct IsIntegralBits<IntegralBitsRotatedView<_IntegralBitsT, _n_bits_per_site>> : std::true_type {};

}

// #######################################################################
// ## IntegralBitsRefinedView                                           ##
// #######################################################################

namespace kstate_op_integral {

template<typename _IntegralBitsT, typename _SiteStateTraitT, unsigned char _n_bits_per_site>
struct IntegralBitsRefinedView {
    static_assert(IsIntegralBits<_IntegralBitsT>::value);
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(_SiteStateTraitT::site_basis_dim() == 2);
    static_assert(_n_bits_per_site >= 1u);
    using IntegralBitsT = _IntegralBitsT;
    using SiteStateTraitT = _SiteStateTraitT;
    using IntegralT = typename IntegralBitsT::IntegralT;
    using SiteStateT = typename SiteStateTraitT::SiteStateT;
    static constexpr unsigned char n_bits_per_site = _n_bits_per_site;
    const IntegralBitsT& integral_bits;
    const kstate_view_amend_spec::RefinedHolder<typename SiteStateTraitT::SiteStateT>& h;
    IntegralT get_number() const noexcept {
        const IntegralT underlying_value = integral_bits.get_number();
        const unsigned new_chunk_number = SiteStateTraitT::get_index(h.v[0]);
        const unsigned char idx_first_bit = n_bits_per_site * h.n;
        const unsigned char n_bits_in_chunk_number = n_bits_per_site;
        const IntegralT refined_value = kstate_op_integral::raw::refine_chunk_number<IntegralT, unsigned>(
                    underlying_value,
                    idx_first_bit,
                    n_bits_in_chunk_number,
                    new_chunk_number);
        return refined_value;
    }
    unsigned char get_n_all_bits() const noexcept {
        return integral_bits.get_n_all_bits();
    };
};

template<typename _IntegralBitsT, typename _SiteStateTraitT, unsigned char _n_bits_per_site>
struct IsIntegralBits<IntegralBitsRefinedView<_IntegralBitsT, _SiteStateTraitT, _n_bits_per_site>> : std::true_type {};

}
