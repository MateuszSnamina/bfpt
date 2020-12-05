#pragma once

#include <monostar_app/interpreted_program_options.hpp>

#include <monostar_hamiltonians/hamiltonian_reference_energies.hpp>

#include <monostar_hamiltonians/hamiltonian_params_fo.hpp>
#include <monostar_hamiltonians/hamiltonian_params_affo.hpp>

#include <bfpt_common/do_common_recipie.hpp>

// #######################################################################
// ## print_{foo,bar,baz}                                               ##
// #######################################################################

namespace monostar_app {

void print_input_data(
    const InterpretedProgramOptions& interpreted_program_options);

void print_results_tree(
    const InterpretedProgramOptions& interpreted_program_options,
    const std::shared_ptr<monostar_hamiltonians::HamiltonianReferenceEnergies> reference_energies,
    const std::optional<bfpt_common::CommonRecipeReceipt>& gs_receipt_optional,
    const std::optional<std::vector<bfpt_common::CommonRecipeReceipt>>& es_receipts_optional,
    const unsigned n_sites);

void print_post_data(
    const InterpretedProgramOptions& interpreted_program_options,
    /*const std::unique_ptr<monostar_app::ReferenceEnergies> reference_energies,*/
    const std::optional<bfpt_common::CommonRecipeReceipt>& gs_receipt_optional,
    const std::optional<std::vector<bfpt_common::CommonRecipeReceipt>>& es_receipts_optional);

void print_theta_opt(
    const monostar_hamiltonians::HamiltonianParamsFo& hamiltonian_fo_params,
    const std::optional<double> user_defined_overrule);

void print_theta_opt(
    const monostar_hamiltonians::HamiltonianParamsAffo& hamiltonian_affo_params,
    const std::optional<double> user_defined_overrule,
    const double average_ss);

void print_decayed_agile_affo_hamiltonian_params(
    const monostar_hamiltonians::HamiltonianParamsAgileAffo& hamiltonian_params_agile_affo,
    const std::optional<double> user_defined_overrule_for_theta_opt,
    const double average_ss);

}  // end of namespace monostar_app
