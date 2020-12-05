#pragma once

#include <monostar_app/interpreted_program_options.hpp>

#include <monostar_hamiltonians/hamiltonian_params_jkl01.hpp>
#include <monostar_hamiltonians/hamiltonian_params_af_fm_site_matrices.hpp>
#include <monostar_hamiltonians/hamiltonian_params_fo_site_matrices.hpp>
#include <monostar_hamiltonians/hamiltonian_params_agile_affo.hpp>

#include <monostar_hamiltonians/hamiltonian_reference_energies.hpp>
#include <chainkernel/operator_kernel.hpp>
#include <monostar_system/monostar_site_state.hpp>
#include <utility/named.hpp>

#include <armadillo>
#include <memory>
#include <vector>

namespace monostar_app {

chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>
get_prepare_hamiltonian_kernel_1(const InterpretedProgramOptions& interpreted_program_options);

chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>
get_prepare_hamiltonian_kernel_12(const InterpretedProgramOptions& interpreted_program_options);

chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>
get_prepare_hamiltonian_kernel_123(const InterpretedProgramOptions& interpreted_program_options);

chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>
get_prepare_hamiltonian_kernel_1234(const InterpretedProgramOptions& interpreted_program_options);

std::vector<utility::Named<arma::cx_mat22>>
get_one_site_metrices_for_average_calculation(const InterpretedProgramOptions& interpreted_program_options);

std::vector<utility::Named<arma::cx_mat44>>
get_two_site_metrices_for_average_calculation(const InterpretedProgramOptions& interpreted_program_options);

std::shared_ptr<monostar_hamiltonians::HamiltonianReferenceEnergies>
get_hamiltonian_reference_energies(const InterpretedProgramOptions& interpreted_program_options);

}  // end of namespace monostar_app

// #include <monostar_system/monostar_site_state.hpp> //TOODO remove?
