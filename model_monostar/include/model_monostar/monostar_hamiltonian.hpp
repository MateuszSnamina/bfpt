#ifndef MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP
#define MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP

#include<bfpt_common/hamiltonian_12.hpp>
#include<model_monostar/monostar_site_state.hpp>

#include<map>

// #######################################################################
// ## Helper function for preparing Hamiltonian12                       ##
// #######################################################################

namespace model_monostar {

std::multimap<bfpt_common::SiteStatePair<MonostarSiteState>, bfpt_common::CoupleInfo<MonostarSiteState>>
prepare_half_off_diag_info_for_af(double J);

std::map<bfpt_common::SiteStatePair<MonostarSiteState>, double>
prepare_diag_info(double J);

bfpt_common::Hamiltonian12<MonostarSiteState>
prepare_hamiltonian_12(double J_classical = 1.0, double J_quantum = 1.0);

} // end of namespace model_monostar

#endif
