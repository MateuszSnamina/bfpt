#pragma once

#include <so_system/so_site_state.hpp>

#include <kstate_impl/kstate_concrete_stl.hpp>
#include <kstate_impl/kstate_concrete_integral.hpp>

#include <iostream>

// #######################################################################
// ## SoKstate                                                          ##
// #######################################################################

namespace so_system {

using SoDynamicStlKstate = kstate_impl::DynamicStlKstate<SoSiteStateTrait>;

template <size_t N>
using StaticStlSoKstate = kstate_impl::StaticStlKstate<SoSiteStateTrait, N>;

using SoDynamicIntegral64Kstate = kstate_impl::DynamicIntegral64Kstate<SoSiteStateTrait, 2u>;

//template<size_t N>
//using StaticIntegralSoKstate = kstate_impl::StaticTwoLevelIntegral64Kstate<SoSiteStateTrait, N, 1u>;
// ^^^^ TODO: UNCOMMENT WHEN IT WILL BE IMPLEMENTED

/*
 *  Select the type you want to use:
 */

//using SoKstate = SoDynamicStlKstate;
//using SoKstate = StaticStlSoKstate<20>;
using SoKstate = SoDynamicIntegral64Kstate;

}  // namespace so_system

// #######################################################################
// ## SoKstateTrait - trait                                             ##
// #######################################################################

namespace so_system {

using SoKstateTrait = kstate_trait::TraitKstate<SoKstate>;

}

// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace so_system {

SoKstate classical_gs_kstate(const unsigned n_sites);

SoKstate classical_es_kstate(const unsigned n_sites);

}  // namespace so_system

// #######################################################################
// ## SoKstate - printing                                               ##
// #######################################################################

namespace so_system {

std::ostream& operator<<(std::ostream& stream, const SoKstate& state);

}  // namespace so_system
