#pragma once

#include<monostar_hamiltonians/hamiltonian_params_fo.hpp>

#include<monostar_system/monostar_site_state.hpp>

#include<bfpt_common/operator_kernel.hpp>

#include<map>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_fo                              ##
// #######################################################################

namespace monostar_hamiltonians {

bfpt_common::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_fo(double Pzz_coef, double Pxz_coef, double Pxx_coef, double orbital_theta);

bfpt_common::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_fo(const HamiltonianParamsFo&, double orbital_theta);

bfpt_common::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_fo(double tau_z_coef, double tau_minus_coef, double orbital_theta);

bfpt_common::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_fo(const HamiltonianParamsFo&, double orbital_theta);

} // end of namespace monostar_hamiltonians
