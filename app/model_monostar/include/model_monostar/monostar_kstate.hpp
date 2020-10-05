#ifndef MODEL_MONOSTAR_MONOSTAR_KSTATE_HPP
#define MODEL_MONOSTAR_MONOSTAR_KSTATE_HPP

#include <model_monostar/monostar_site_state.hpp>

#include <kstate/kstate_concrete.hpp>

#include <iostream>

// #######################################################################
// ## DynamicMonostarKstate                                             ##
// #######################################################################

namespace model_monostar {

using DynamicMonostarKstate = kstate::DynamicKstate<kstate::TraitSiteState<MonostarSiteState>>;

}  // namespace model_monostar


// #######################################################################
// ## DynamicMonostarKstate - trait                                     ##
// #######################################################################

namespace model_monostar {

using DynamicMonostarKstateTrait = kstate::TraitKstate<DynamicMonostarKstate>;

}


// #######################################################################
// ## classical_gs_kstate, classical_es_kstate                          ##
// #######################################################################

namespace model_monostar {

DynamicMonostarKstate classical_gs_kstate(const unsigned n_sites);

DynamicMonostarKstate classical_es_kstate(const unsigned n_sites);

}  // namespace model_monostar


// #######################################################################
// ## DynamicMonostarKstate - printing                                  ##
// #######################################################################

namespace model_monostar {

std::ostream& operator<<(std::ostream& stream, const DynamicMonostarKstate& state);

}  // namespace model_monostar

#endif
