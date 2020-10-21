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

template<typename _IntegralBitsT>
struct IntegralBitsRotatedView {
    static_assert(IsIntegralBits<_IntegralBitsT>::value);
    using IntegralBitsT = _IntegralBitsT;
    using IntegralT = typename IntegralBitsT::IntegralT;
    const IntegralBitsT& integral_bits;
    const kstate_view_amend_spec::RotateHolder& h;
    IntegralT get_number() const noexcept {
        return kstate_op_integral::raw::rotate(integral_bits.get_number(), integral_bits.get_n_all_bits(), h.n);
    }
    unsigned char get_n_all_bits() const noexcept {
        return integral_bits.get_n_all_bits();
    };
};

template<typename _IntegralBitsT>
struct IsIntegralBits<IntegralBitsRotatedView<_IntegralBitsT>> : std::true_type {};

}

// #######################################################################
// ## IntegralBitsRefinedView                                           ##
// #######################################################################

namespace kstate_op_integral {

template<typename _IntegralBitsT, typename _SiteStateTraitT>
struct IntegralBitsRefinedView {
    static_assert(IsIntegralBits<_IntegralBitsT>::value);
    static_assert(kstate_trait::IsTraitSiteState<_SiteStateTraitT>::value);
    static_assert(_SiteStateTraitT::is_site_state_trait);
    static_assert(_SiteStateTraitT::site_basis_dim() == 2);
    using IntegralBitsT = _IntegralBitsT;
    using SiteStateTraitT = _SiteStateTraitT;
    using IntegralT = typename IntegralBitsT::IntegralT;
    using SiteStateT = typename SiteStateTraitT::SiteStateT;
    const IntegralBitsT& integral_bits;
    const kstate_view_amend_spec::RefinedHolder<typename SiteStateTraitT::SiteStateT>& h;
    IntegralT get_number() const noexcept {
        const unsigned new_value_unsigned = SiteStateTraitT::get_index(h.v[1]);
        const bool new_value_bool = bool_from_unsiged(new_value_unsigned);
        return kstate_op_integral::raw::refine(integral_bits.get_number(), new_value_bool, h.n) ;
    }
    unsigned char get_n_all_bits() const noexcept {
        return integral_bits.get_n_all_bits();
    };
};

template<typename _IntegralBitsT, typename _SiteStateTraitT>
struct IsIntegralBits<IntegralBitsRefinedView<_IntegralBitsT, _SiteStateTraitT>> : std::true_type {};

}


//TODO remove!?
//template<typename ViewT>
//static auto refined_view(const ViewT& v, const kstate_view_amend_spec::RefinedHolder<typename SiteStateTraitT::SiteStateT>& h) noexcept {
//    return kstate_op_range::raw::refined(v, h);
//}
//template<typename ViewT>
//static auto rotated_view(const ViewT& v, const kstate_view_amend_spec::RotateHolder& h) noexcept {
//    return kstate_op_range::raw::rotated(v, h);





