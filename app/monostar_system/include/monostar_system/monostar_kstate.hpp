#pragma once

#include <monostar_system/monostar_site_state.hpp>

#include <kstate_impl/kstate_concrete_stl.hpp>

#include <iostream>

// #######################################################################
// ## DynamicMonostarKstate                                             ##
// #######################################################################

namespace monostar_system {

using DynamicMonostarKstate = kstate_impl::DynamicStlKstate<MonostarSiteStateTrait>;

template<size_t N>
using StaticMonostarKstate = kstate_impl::StaticStlKstate<MonostarSiteStateTrait, N>;

using Static20MonostarKstate = StaticMonostarKstate<20>;

}  // namespace monostar_system


// #######################################################################
// ## DynamicMonostarKstate - trait                                     ##
// #######################################################################

namespace monostar_system {

using DynamicMonostarKstateTrait = kstate_trait::TraitKstate<DynamicMonostarKstate>;

}


// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace monostar_system {

DynamicMonostarKstate classical_gs_kstate(const unsigned n_sites);

DynamicMonostarKstate classical_es_kstate(const unsigned n_sites);

}  // namespace monostar_system


// #######################################################################
// ## DynamicMonostarKstate - printing                                  ##
// #######################################################################

namespace monostar_system {

std::ostream& operator<<(std::ostream& stream, const DynamicMonostarKstate& state);

}  // namespace monostar_system
