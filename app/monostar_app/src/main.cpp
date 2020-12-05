#include <monostar_app/raw_program_options.hpp>
#include <monostar_app/interpreted_program_options.hpp>
#include <monostar_app/main_printing.hpp>
#include <monostar_app/get_model_bits.hpp>

#include <monostar_system/monostar_basis.hpp>
#include <monostar_system/monostar_kstate.hpp>
#include <monostar_system/monostar_site_state.hpp>

#include <chainkernel/operator_kernel.hpp>
#include <koperator_impl/kernel_driven_kstate_operator_matrix.hpp>
#include <kpopulator_impl/kernel_driven_kstate_basis_populator.hpp>

#include <bfpt_common/do_common_recipie.hpp>

#include <armadillo>

#include <iostream>
#include <iomanip>

#include <cassert>
#include <cstdlib>

// #######################################################################
// ## main - helpers                                                    ##
// #######################################################################

bfpt_common::CommonRecipeReceipt bfpt_gs(
    const chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>& hamiltonian_kernel_1,
    const chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>& hamiltonian_kernel_12,
    const chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>& hamiltonian_kernel_123,
    const chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>& hamiltonian_kernel_1234,
    const size_t n_sites, const unsigned max_pt_order,
    const bfpt_common::CommonRecipePrintFlags& print_flags,
    const std::vector<utility::Named<arma::cx_mat22>>& one_site_metrices_for_average_calculation,
    const std::vector<utility::Named<arma::cx_mat44>>& two_sites_metrices_for_average_calculation,
    const unsigned n_threads) {
    using KstateT = monostar_system::MonostarKstate;
    using KstateTraitT = monostar_system::MonostarKstateTrait;
    using KpopulatorT = kpopulator_impl::KernelDrivenKstateBasisPopulator<KstateTraitT>;
    using KpopulatorTraitT = kpopulator_trait::TraitKpopulator<KpopulatorT>;
    using KoperatorT = koperator_impl::KernelDrivenKstateOperatorMatrix<KstateTraitT>;
    using KoperatorTraitT = koperator_trait::TraitKoperator<KoperatorT>;
    using BasisT = monostar_system::MonostarKstateBasis;
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(monostar_system::classical_gs_kstate(n_sites)));
    const KpopulatorT kstate_populator{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12, hamiltonian_kernel_123, hamiltonian_kernel_1234};
    const KoperatorT kstate_hamiltonian{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12, hamiltonian_kernel_123, hamiltonian_kernel_1234};
    return bfpt_common::do_common_recipe<KstateTraitT, KpopulatorTraitT, KoperatorTraitT, 2u>(
               kstate_populator, kstate_hamiltonian,
               basis, max_pt_order,
               0,
               print_flags,
               one_site_metrices_for_average_calculation,
               two_sites_metrices_for_average_calculation,
               "[gs] ",
               n_threads)
        .unwrap();
}

bfpt_common::CommonRecipeReceipt bfpt_kn_es(
    const chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait>& hamiltonian_kernel_1,
    const chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait>& hamiltonian_kernel_12,
    const chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait>& hamiltonian_kernel_123,
    const chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait>& hamiltonian_kernel_1234,
    const size_t n_sites, const unsigned max_pt_order, const unsigned k_n,
    const bfpt_common::CommonRecipePrintFlags& print_flags,
    const std::vector<utility::Named<arma::cx_mat22>>& one_site_metrices_for_average_calculation,
    const std::vector<utility::Named<arma::cx_mat44>>& two_sites_metrices_for_average_calculation,
    const unsigned n_threads) {
    using KstateT = monostar_system::MonostarKstate;
    using KstateTraitT = monostar_system::MonostarKstateTrait;
    using KpopulatorT = kpopulator_impl::KernelDrivenKstateBasisPopulator<KstateTraitT>;
    using KpopulatorTraitT = kpopulator_trait::TraitKpopulator<KpopulatorT>;
    using KoperatorT = koperator_impl::KernelDrivenKstateOperatorMatrix<KstateTraitT>;
    using KoperatorTraitT = koperator_trait::TraitKoperator<KoperatorT>;
    using BasisT = monostar_system::MonostarKstateBasis;
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(monostar_system::classical_es_kstate(n_sites)));
    const KpopulatorT kstate_populator{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12, hamiltonian_kernel_123, hamiltonian_kernel_1234};
    const KoperatorT kstate_hamiltonian{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12, hamiltonian_kernel_123, hamiltonian_kernel_1234};
    return bfpt_common::do_common_recipe<KstateTraitT, KpopulatorTraitT, KoperatorTraitT, 2u>(
               kstate_populator, kstate_hamiltonian,
               basis, max_pt_order,
               k_n,
               print_flags,
               one_site_metrices_for_average_calculation,
               two_sites_metrices_for_average_calculation,
               "[es (" + std::to_string(k_n) + ")] ",
               n_threads)
        .unwrap();
}

// #######################################################################
// ## main                                                              ##
// #######################################################################

int main(int argc, char** argv) {
    using namespace monostar_app;
    try {
        // ******************************************************************
        const RawProgramOptions raw_program_options = grep_program_options(argc, argv);
        const InterpretedProgramOptions interpreted_program_options = interpret_program_options(raw_program_options);
        // ******************************************************************
        print_input_data(interpreted_program_options);
        // ******************************************************************
        if (interpreted_program_options.model_type == ModelType::FO) {
            print_theta_opt(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
        }
        // ******************************************************************
        const chainkernel::OperatorKernel1<monostar_system::MonostarSiteStateTrait> hamiltonian_kernel_1 =
            get_prepare_hamiltonian_kernel_1(interpreted_program_options);
        const chainkernel::OperatorKernel12<monostar_system::MonostarSiteStateTrait> hamiltonian_kernel_12 =
            get_prepare_hamiltonian_kernel_12(interpreted_program_options);
        const chainkernel::OperatorKernel123<monostar_system::MonostarSiteStateTrait> hamiltonian_kernel_123 =
            get_prepare_hamiltonian_kernel_123(interpreted_program_options);
        const chainkernel::OperatorKernel1234<monostar_system::MonostarSiteStateTrait> hamiltonian_kernel_1234 =
            get_prepare_hamiltonian_kernel_1234(interpreted_program_options);
        const std::vector<utility::Named<arma::cx_mat22>> one_site_metrices_for_average_calculation =
            get_one_site_metrices_for_average_calculation(interpreted_program_options);
        const std::vector<utility::Named<arma::cx_mat44>> two_sites_metrices_for_average_calculation =
            get_two_site_metrices_for_average_calculation(interpreted_program_options);
        const std::shared_ptr<monostar_hamiltonians::HamiltonianReferenceEnergies> reference_energies =
            get_hamiltonian_reference_energies(interpreted_program_options);
        // ******************************************************************
        const std::optional<bfpt_common::CommonRecipeReceipt> gs_receipt =
            [&]() -> std::optional<bfpt_common::CommonRecipeReceipt> {
            if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
                std::cout << "------------------------------------------" << std::endl;
                return bfpt_gs(
                    hamiltonian_kernel_1,
                    hamiltonian_kernel_12,
                    hamiltonian_kernel_123,
                    hamiltonian_kernel_1234,
                    interpreted_program_options.n_sites, interpreted_program_options.n_pt,
                    interpreted_program_options.print_flags,
                    one_site_metrices_for_average_calculation,
                    two_sites_metrices_for_average_calculation,
                    interpreted_program_options.n_threads);
                std::cout << "------------------------------------------" << std::endl;
            }
            return std::nullopt;
        }();
        // ******************************************************************
        const std::optional<std::vector<bfpt_common::CommonRecipeReceipt>> es_receipts =
            [&]() -> std::optional<std::vector<bfpt_common::CommonRecipeReceipt>> {
            std::cout << "------------------------------------------" << std::endl;
            if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
                const auto es_momentum_range_sapn = es_momentum_domain_variant_to_momentum_range_sapn(
                    interpreted_program_options.es_momentum_domain,
                    interpreted_program_options.n_sites);
                std::vector<bfpt_common::CommonRecipeReceipt> es_receipts_builder;
                for (unsigned k_n = es_momentum_range_sapn.first; k_n < es_momentum_range_sapn.second; k_n++) {
                    std::cout << "[PROGRESS] "
                              << "solving n_k: " << k_n << std::endl;
                    const auto es_energy = bfpt_kn_es(
                        hamiltonian_kernel_1,
                        hamiltonian_kernel_12,
                        hamiltonian_kernel_123,
                        hamiltonian_kernel_1234,
                        interpreted_program_options.n_sites, interpreted_program_options.n_pt, k_n,
                        interpreted_program_options.print_flags,
                        one_site_metrices_for_average_calculation,
                        two_sites_metrices_for_average_calculation,
                        interpreted_program_options.n_threads);
                    es_receipts_builder.push_back(es_energy);
                    std::cout << "------------------------------------------" << std::endl;
                }
                return es_receipts_builder;
            }
            return std::nullopt;
        }();
        // ******************************************************************
        print_results_tree(
            interpreted_program_options,
            reference_energies,
            gs_receipt,
            es_receipts,
            interpreted_program_options.n_sites);
        print_post_data(
            interpreted_program_options,
            /*reference_energies,*/
            gs_receipt,
            es_receipts);
        // ******************************************************************
    } catch (std::exception& e) {
        std::cerr << "[ERROR  ] Abnormal termination!" << std::endl;
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
