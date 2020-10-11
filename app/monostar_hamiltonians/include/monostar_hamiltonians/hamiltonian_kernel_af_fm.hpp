#pragma once

#include<monostar_hamiltonians/hamiltonian_params_af_fm.hpp>

#include<monostar_system/monostar_site_state.hpp>

#include<bfpt_common/operator_kernel.hpp>

#include<kstate_trait/trait_site_state.hpp>

#include<map>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_{af,fm}                         ##
// #######################################################################

namespace monostar_hamiltonians {

bfpt_common::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_af(double J_classical = 1.0, double J_quantum = 1.0);

bfpt_common::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_af(const HamiltonianParamsAfFm&);

bfpt_common::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_fm(double J_classical = 1.0, double J_quantum = 1.0);

bfpt_common::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_fm(const HamiltonianParamsAfFm&);

bfpt_common::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_af_fm(double B = 0.0);

bfpt_common::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_af_fm(const HamiltonianParamsAfFm&);

} // end of namespace monostar_hamiltonians
