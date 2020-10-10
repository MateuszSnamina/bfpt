#pragma once

#include <monostar_system/monostar_site_state.hpp>

#include <kstate/kstate_concrete.hpp>

#include <iostream>

// #######################################################################
// ## DynamicMonostarKstate                                             ##
// #######################################################################

namespace monostar_system {

using DynamicMonostarKstate = kstate::DynamicKstate<kstate::TraitSiteState<MonostarSiteState>>;

}  // namespace monostar_system


// #######################################################################
// ## DynamicMonostarKstate - trait                                     ##
// #######################################################################

namespace monostar_system {

using DynamicMonostarKstateTrait = kstate::TraitKstate<DynamicMonostarKstate>;

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
