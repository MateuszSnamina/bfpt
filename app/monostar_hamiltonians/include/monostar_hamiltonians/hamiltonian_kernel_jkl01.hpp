#pragma once

#include <monostar_hamiltonians/hamiltonian_params_jkl01.hpp>

#include <monostar_system/monostar_site_state.hpp>

#include <chainkernel/operator_kernel.hpp>

#include <kstate_trait/trait_site_state.hpp>  //TODO remove?

#include <map>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12,123,1234}_jkl01                  ##
// #######################################################################

namespace monostar_hamiltonians {

chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_jkl01();

chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1_jkl01(const HamiltonianParamsJkl01&);

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_jkl01();

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_12_jkl01(const HamiltonianParamsJkl01&);

chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_123_jkl01(double L, double L_1);

chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_123_jkl01(const HamiltonianParamsJkl01&);

chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1234_jkl01(double J, double J_0, double J_1, double K, double K_0, double K_1);

chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>
prepare_hamiltonian_kernel_1234_jkl01(const HamiltonianParamsJkl01&);

}  // end of namespace monostar_hamiltonians
