#pragma once

#include <so_hamiltonians/hamiltonian_params_af_fo.hpp>

#include <so_system/so_site_state.hpp>

#include <chainkernel/operator_kernel.hpp>

#include <map>

// #######################################################################
// ## prepare_hamiltonian_kernel_{1,12}_fo                              ##
// #######################################################################

namespace so_hamiltonians {

chainkernel::OperatorKernel12<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_12_fo(
        double ss_coef,
        double Pzz_coef, double Pxz_coef, double Pxx_coef,
        double ss_Pzz_coef, double ss_Pxz_coef, double ss_Pxx_coef,
        double orbital_theta);

chainkernel::OperatorKernel12<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_12_fo(const HamiltonianParamsAfFo&, double orbital_theta);

chainkernel::OperatorKernel1<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_1_fo(
        double tau_z_coef, double tau_minus_coef,
        double orbital_theta);

chainkernel::OperatorKernel1<so_system::SoSiteStateTrait>
prepare_hamiltonian_kernel_1_fo(const HamiltonianParamsAfFo&, double orbital_theta);

}  // end of namespace so_hamiltonians
