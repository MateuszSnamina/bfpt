#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/interpreted_program_options.hpp>
#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>
#include <model_monostar/hamiltonian_kernel_af_fm.hpp>
#include <model_monostar/hamiltonian_params_af_fm_site_matrices.hpp>
#include <model_monostar/hamiltonian_reference_energies_af_fm.hpp>
#include <model_monostar/hamiltonian_kernel_fo.hpp>
#include <model_monostar/hamiltonian_params_fo_site_matrices.hpp>
#include <model_monostar/hamiltonian_reference_energies_fo.hpp>
#include <model_monostar/get_orbital_theta.hpp>
#include <model_monostar/main_printing.hpp>

#include <bfpt_common/operator_kernel.hpp>
#include <bfpt_common/kernel_driven_kstate_basis_populator.hpp>
#include <bfpt_common/kernel_driven_kstate_operator_matrix.hpp>
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
        const bfpt_common::OperatorKernel1<model_monostar::MonostarSiteStateTrait>& hamiltonian_kernel_1,
        const bfpt_common::OperatorKernel12<model_monostar::MonostarSiteStateTrait>& hamiltonian_kernel_12,
        const size_t n_sites, const unsigned max_pt_order,
        const bfpt_common::CommonRecipePrintFlags& print_flags,
        const std::vector<utility::Named<arma::cx_mat22>>& one_site_metrices_for_average_calculation,
        const std::vector<utility::Named<arma::cx_mat44>>& two_sites_metrices_for_average_calculation,
        unsigned n_threads) {
    using KstateT = model_monostar::DynamicMonostarKstate;
    using KstateTraitT = model_monostar::DynamicMonostarKstateTrait;
    using BasisT = model_monostar::DynamicMonostarKstateBasis;
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(model_monostar::classical_gs_kstate(n_sites)));
    const bfpt_common::KernelDrivenKstateBasisPopulator<KstateTraitT> kstate_populator{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12};
    const bfpt_common::KernelDrivenKstateOperatorMatrix<KstateTraitT> kstate_hamiltonian{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12};
    return bfpt_common::do_common_recipe<KstateTraitT>(kstate_populator, kstate_hamiltonian,
                                                       basis, max_pt_order,
                                                       0,
                                                       print_flags,
                                                       one_site_metrices_for_average_calculation,
                                                       two_sites_metrices_for_average_calculation,
                                                       "[gs] ",
                                                       n_threads).unwrap();
}

bfpt_common::CommonRecipeReceipt bfpt_kn_es(
        const bfpt_common::OperatorKernel1<model_monostar::MonostarSiteStateTrait>& hamiltonian_kernel_1,
        const bfpt_common::OperatorKernel12<model_monostar::MonostarSiteStateTrait>& hamiltonian_kernel_12,
        const size_t n_sites, const unsigned max_pt_order, const unsigned k_n,
        const bfpt_common::CommonRecipePrintFlags& print_flags,
        const std::vector<utility::Named<arma::cx_mat22>>& one_site_metrices_for_average_calculation,
        const std::vector<utility::Named<arma::cx_mat44>>& two_sites_metrices_for_average_calculation,
        unsigned n_threads) {
    using KstateT = model_monostar::DynamicMonostarKstate;
    using KstateTraitT = model_monostar::DynamicMonostarKstateTrait;
    using BasisT = model_monostar::DynamicMonostarKstateBasis;
    BasisT basis{n_sites};
    basis.add_element(std::make_shared<KstateT>(model_monostar::classical_es_kstate(n_sites)));
    const bfpt_common::KernelDrivenKstateBasisPopulator<KstateTraitT> kstate_populator{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12};
    const bfpt_common::KernelDrivenKstateOperatorMatrix<KstateTraitT> kstate_hamiltonian{n_sites, hamiltonian_kernel_1, hamiltonian_kernel_12};
    return bfpt_common::do_common_recipe<KstateTraitT>(kstate_populator, kstate_hamiltonian,
                                                       basis, max_pt_order,
                                                       k_n,
                                                       print_flags,
                                                       one_site_metrices_for_average_calculation,
                                                       two_sites_metrices_for_average_calculation,
                                                       "[es (" + std::to_string(k_n) + ")] ",
                                                       n_threads).unwrap();
}

// #######################################################################
// ## main                                                              ##
// #######################################################################

int main(int argc, char** argv) {

    //    {
    //        arma::cx_mat22 P_z_in_zx_basis = OrbitalSiteMatrices::get_P_z_in_zx_basis();
    //        arma::cx_mat22 P_x_in_zx_basis = OrbitalSiteMatrices::get_P_x_in_zx_basis();
    //        arma::cx_mat22 P_plus_in_zx_basis = OrbitalSiteMatrices::get_P_plus_in_zx_basis();
    //        arma::cx_mat22 P_minus_in_zx_basis = OrbitalSiteMatrices::get_P_minus_in_zx_basis();
    //        arma::cx_mat22 tau_z_in_zx_basis = OrbitalSiteMatrices::get_tau_z_in_zx_basis();
    //        arma::cx_mat22 tau_x_in_zx_basis = OrbitalSiteMatrices::get_tau_x_in_zx_basis();
    //        arma::cx_mat22 tau_plus_in_zx_basis = OrbitalSiteMatrices::get_tau_plus_in_zx_basis();
    //        arma::cx_mat22 tau_minus_in_zx_basis = OrbitalSiteMatrices::get_tau_minus_in_zx_basis();

    //        arma::cx_mat22 beta_from_zx_to_ge = OrbitalSiteMatrices::get_beta_from_zx_to_ge(0.253585);

    //        arma::cx_mat22 P_z_in_ge_basis = OrbitalSiteMatrices::get_P_z_in_ge_basis(0.253585);
    //        arma::cx_mat22 P_x_in_ge_basis = OrbitalSiteMatrices::get_P_x_in_ge_basis(0.253585);
    //        arma::cx_mat22 P_plus_in_ge_basis = OrbitalSiteMatrices::get_P_plus_in_ge_basis(0.253585);
    //        arma::cx_mat22 P_minus_in_ge_basis = OrbitalSiteMatrices::get_P_minus_in_ge_basis(0.253585);
    //        arma::cx_mat22 tau_z_in_ge_basis = OrbitalSiteMatrices::get_tau_z_in_ge_basis(0.253585);
    //        arma::cx_mat22 tau_x_in_ge_basis = OrbitalSiteMatrices::get_tau_x_in_ge_basis(0.253585);
    //        arma::cx_mat22 tau_plus_in_ge_basis = OrbitalSiteMatrices::get_tau_plus_in_ge_basis(0.253585);
    //        arma::cx_mat22 tau_minus_in_ge_basis = OrbitalSiteMatrices::get_tau_minus_in_ge_basis(0.253585);
    //    }

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
        const auto hamiltonian_kernel_1 = [&interpreted_program_options]() {
            switch (interpreted_program_options.model_type) {
            case ModelType::AF:
            case ModelType::FM:
                return model_monostar::prepare_hamiltonian_kernel_1_af_fm(interpreted_program_options.hamiltonian_params_af_fm);
            case ModelType::FO:
            {
                const double orbital_theta_to_use = get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
                return model_monostar::prepare_hamiltonian_kernel_1_fo(interpreted_program_options.hamiltonian_params_fo, orbital_theta_to_use);
            }
            default:
                throw std::domain_error("Invalid model_type enum value.");
            }
        }();
        const auto hamiltonian_kernel_12 = [&interpreted_program_options]() {
            switch (interpreted_program_options.model_type) {
            case ModelType::AF:
                return model_monostar::prepare_hamiltonian_kernel_12_af(interpreted_program_options.hamiltonian_params_af_fm);
            case ModelType::FM:
                return model_monostar::prepare_hamiltonian_kernel_12_fm(interpreted_program_options.hamiltonian_params_af_fm);
            case ModelType::FO:
            {
                const double orbital_theta_to_use = get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
                return model_monostar::prepare_hamiltonian_kernel_12_fo(interpreted_program_options.hamiltonian_params_fo, orbital_theta_to_use);
            }
            default:
                throw std::domain_error("Invalid model_type enum value.");
            };
        }();
        const auto one_site_metrices_for_average_calculation =  [&interpreted_program_options]() {
            switch (interpreted_program_options.model_type) {
            case ModelType::AF:
            case ModelType::FM:
                return OneSiteSpinNamedMatrices::site_matrices_for_average_calculations_af_fm();
            case ModelType::FO:
            {
                const double orbital_theta_to_use = get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
                return OneSiteOrbitalNamedMatrices::matrices_for_average_calculations(orbital_theta_to_use);
            }
            default:
                throw std::domain_error("Invalid model_type enum value.");
            };
        }();
        const auto two_sites_metrices_for_average_calculation =  [&interpreted_program_options]() {
            switch (interpreted_program_options.model_type) {
            case ModelType::AF:
                return TwoSitesSpinNamedMatrices::matrices_for_average_calculations_af();
            case ModelType::FM:
                return TwoSitesSpinNamedMatrices::matrices_for_average_calculations_fm();
            case ModelType::FO:
            {
                const double orbital_theta_to_use = get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
                return TwoSitesOrbitalNamedMatrices::matrices_for_average_calculations(orbital_theta_to_use);
            }
            default:
                throw std::domain_error("Invalid model_type enum value.");
            };
        }();

        // ******************************************************************
        const std::shared_ptr<model_monostar::HamiltonianReferenceEnergies> reference_energies =
                [&interpreted_program_options]() {
            switch (interpreted_program_options.model_type) {
            case ModelType::AF:
                return std::dynamic_pointer_cast<model_monostar::HamiltonianReferenceEnergies>(
                            std::make_shared<model_monostar::HamiltonianReferenceEnergiesAf>(
                                interpreted_program_options.n_sites, interpreted_program_options.hamiltonian_params_af_fm));
            case ModelType::FM:
                return std::dynamic_pointer_cast<model_monostar::HamiltonianReferenceEnergies>(
                            std::make_shared<model_monostar::HamiltonianReferenceEnergiesFm>(
                                interpreted_program_options.n_sites, interpreted_program_options.hamiltonian_params_af_fm));
            case ModelType::FO:
            {
                const double orbital_theta_to_use = get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
                return std::dynamic_pointer_cast<model_monostar::HamiltonianReferenceEnergies>(
                            std::make_shared<model_monostar::HamiltonianReferenceEnergiesFo>(
                                interpreted_program_options.n_sites, interpreted_program_options.hamiltonian_params_fo, orbital_theta_to_use));
            }
            default:
                throw std::domain_error("Invalid model_type enum value.");
            }
        }();
        // ******************************************************************
        const std::optional<bfpt_common::CommonRecipeReceipt> gs_receipt =
                [&]() -> std::optional<bfpt_common::CommonRecipeReceipt> {
                if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
                std::cout << "------------------------------------------" << std::endl;
                return bfpt_gs(
                    hamiltonian_kernel_1,
                    hamiltonian_kernel_12,
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
            std::cout << "[PROGRESS] " << "solving n_k: " << k_n << std::endl;
            const auto es_energy = bfpt_kn_es(
                        hamiltonian_kernel_1,
                        hamiltonian_kernel_12,
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
