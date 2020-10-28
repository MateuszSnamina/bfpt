#pragma once

#include <monostar_system/monostar_site_state.hpp>

#include <kstate_impl/kstate_concrete_stl.hpp>
#include <kstate_impl/kstate_concrete_integral.hpp>

#include <iostream>

// #######################################################################
// ## DynamicMonostarKstate                                             ##
// #######################################################################

namespace monostar_system {

using MonostarDynamicStlKstate = kstate_impl::DynamicStlKstate<MonostarSiteStateTrait>;

template<size_t N>
using StaticStlMonostarKstate = kstate_impl::StaticStlKstate<MonostarSiteStateTrait, N>;

using MonostarDynamicIntegral64Kstate = kstate_impl::DynamicTwoLevelIntegral64Kstate<MonostarSiteStateTrait>;

//template<size_t N>
//using StaticIntegralMonostarKstate = kstate_impl::StaticTwoLevelIntegral64Kstate<MonostarSiteStateTrait, N>;
// ^^^^ TODO: UNCOMMENT WHEN IT WILL BE IMPLEMENTED

/*
 *  Select the type you want to use:
 */

//using MonostarKstate = MonostarDynamicStlKstate;
//using MonostarKstate = StaticStlMonostarKstate<20>;
using MonostarKstate = MonostarDynamicIntegral64Kstate;

}  // namespace monostar_system


// #######################################################################
// ## DynamicMonostarKstate - trait                                     ##
// #######################################################################

namespace monostar_system {

using MonostarKstateTrait = kstate_trait::TraitKstate<MonostarKstate>;

}


// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace monostar_system {

MonostarKstate classical_gs_kstate(const unsigned n_sites);

MonostarKstate classical_es_kstate(const unsigned n_sites);

}  // namespace monostar_system


// #######################################################################
// ## DynamicMonostarKstate - printing                                  ##
// #######################################################################

namespace monostar_system {

std::ostream& operator<<(std::ostream& stream, const MonostarKstate& state);

}  // namespace monostar_system
