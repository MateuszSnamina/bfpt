#ifndef MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP
#define MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_HPP

#include<bfpt_common/hamiltonian_kernel.hpp>
#include<model_monostar/monostar_site_state.hpp>

#include<map>

// #######################################################################
// ## prepare_hamiltonian_12_{af,fm}                                    ##
// #######################################################################

namespace model_monostar {

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_12_af(double J_classical = 1.0, double J_quantum = 1.0);

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_12_fm(double J_classical = 1.0, double J_quantum = 1.0);

} // end of namespace model_monostar

#endif
