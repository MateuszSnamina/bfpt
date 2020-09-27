#ifndef MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_KERNEL_HPP
#define MODEL_MONOSTAR_MONOSTAR_HAMILTONIAN_KERNEL_HPP

#include<bfpt_common/hamiltonian_kernel.hpp>
#include<model_monostar/monostar_site_state.hpp>

#include<map>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_{af,fm}                         ##
// #######################################################################

namespace model_monostar {

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_af(double J_classical = 1.0, double J_quantum = 1.0); //TODO change name: add kernel -> hamiltonian_kernel

bfpt_common::HamiltonianKernel12<MonostarSiteState>
prepare_hamiltonian_kernel_12_fm(double J_classical = 1.0, double J_quantum = 1.0); //TODO change name: add hamiltonian -> hamiltonian_kernel

bfpt_common::HamiltonianKernel1<MonostarSiteState>
prepare_hamiltonian_kernel_1_af_fm(double B = 0.0);

} // end of namespace model_monostar

#endif
