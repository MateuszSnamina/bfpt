#include <model_monostar/raw_program_options.hpp>
#include <model_monostar/interpreted_program_options.hpp>

#include <model_monostar/get_orbital_theta.hpp>
#include <model_monostar/hamiltonian_kernel_af_fm.hpp>
#include <model_monostar/hamiltonian_kernel_fo.hpp>
#include <model_monostar/hamiltonian_reference_energies_af_fm.hpp>
#include <model_monostar/hamiltonian_reference_energies_fo.hpp>
#include <model_monostar/monostar_basis.hpp>
#include <model_monostar/monostar_kstate.hpp>
#include <model_monostar/monostar_site_state.hpp>

#include <bfpt_common/operator_kernel.hpp>
#include <bfpt_common/kernel_driven_kstate_basis_populator.hpp>
#include <bfpt_common/kernel_driven_kstate_operator_matrix.hpp>
#include <bfpt_common/do_common_recipie.hpp>

#include <extensions/range_streamer.hpp>
#include <extensions/stream_fromat_stacker.hpp>

#include <armadillo>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/irange.hpp>

#include <iostream>

#include <cassert>
#include <cstdlib>

#include<model_monostar/hamiltonian_params_fo_site_matrices.hpp> //TODO remove, debug sake
#include<model_monostar/hamiltonian_params_af_fm_site_matrices.hpp> //TODO move to the right place
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
// ## print_foo                                                         ##
// #######################################################################

void print_input_data(const InterpretedProgramOptions& interpreted_program_options) {
    using namespace extension::boost::stream_pragma;
    // Stream RAII:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    // Print:
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_sites                            = " << interpreted_program_options.n_sites << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_pt                               = " << interpreted_program_options.n_pt << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] model_type                         = " << interpreted_program_options.model_type << std::endl;
    if (interpreted_program_options.model_type == ModelType::AF || interpreted_program_options.model_type == ModelType::FM){
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::J_classical     = " << interpreted_program_options.hamiltonian_params_af_fm.get_J_classical() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::J_quantum       = " << interpreted_program_options.hamiltonian_params_af_fm.get_J_quantum() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_af_fm::B               = " << interpreted_program_options.hamiltonian_params_af_fm.get_B() << std::endl;
    }
    if (interpreted_program_options.model_type == ModelType::FO) {
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::tau_z_coef         = " << interpreted_program_options.hamiltonian_params_fo.get_tau_z_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::tau_minus_coef     = " << interpreted_program_options.hamiltonian_params_fo.get_tau_minus_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pzz_coef           = " << interpreted_program_options.hamiltonian_params_fo.get_Pzz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pxz_coef           = " << interpreted_program_options.hamiltonian_params_fo.get_Pxz_coef() << std::endl;
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] hamiltonian_fo::Pxx_coef           = " << interpreted_program_options.hamiltonian_params_fo.get_Pxx_coef() << std::endl;
        if (interpreted_program_options.orbital_theta) {
            std::cout << "[INFO   ] [PROGRAM_OPTIONS] reference orbital theta            = " << *interpreted_program_options.orbital_theta << std::endl;
        } else {
            std::cout << "[INFO   ] [PROGRAM_OPTIONS] reference orbital theta            = " << "<auto: let the program choose the optimal value>" << std::endl;
        }
    }
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] run_type                           = " << interpreted_program_options.run_type << std::endl;
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        std::cout << "[INFO   ] [PROGRAM_OPTIONS] es_momentum_domain                 = " << interpreted_program_options.es_momentum_domain << std::endl;
    }
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::unpopulated_basis_flag      = " << interpreted_program_options.print_flags.print_unpopulated_basis_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::unpopulated_basis_size_flag = " << interpreted_program_options.print_flags.print_unpopulated_basis_size_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::populated_basis_flag        = " << interpreted_program_options.print_flags.print_populated_basis_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::populated_basis_size_flag   = " << interpreted_program_options.print_flags.print_populated_basis_size_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::hamiltonian_stats           = " << interpreted_program_options.print_flags.print_hamiltonian_stats << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::sp_hamiltonian_flag         = " << interpreted_program_options.print_flags.print_sp_hamiltonian_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::hamiltonian_flag            = " << interpreted_program_options.print_flags.print_hamiltonian_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::eigen_values_flag           = " << interpreted_program_options.print_flags.print_eigen_values_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::eigen_vectors_flag          = " << interpreted_program_options.print_flags.print_eigen_vectors_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::pretty_vectors_flag         = " << interpreted_program_options.print_flags.print_pretty_vectors_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::density_operator_flag       = " << interpreted_program_options.print_flags.print_density_operator_flag << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::pretty_min_max_n_kstates    = " << "[" << interpreted_program_options.print_flags.print_pretty_min_max_n_kstates.first << ":" << interpreted_program_options.print_flags.print_pretty_min_max_n_kstates.second << ")" << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] print::pretty_probability_treshold = " << interpreted_program_options.print_flags.print_pretty_probability_treshold << std::endl;
    std::cout << "[INFO   ] [PROGRAM_OPTIONS] n_threads                          = " << interpreted_program_options.n_threads << std::endl;
}


void print_results_tree(
        const InterpretedProgramOptions& interpreted_program_options,
        const std::shared_ptr<model_monostar::HamiltonianReferenceEnergies> reference_energies,
        const std::optional<bfpt_common::CommonRecipeReceipt>& gs_receipt_optional,
        const std::optional<std::vector<bfpt_common::CommonRecipeReceipt>>& es_receipts_optional) {
    // Helpers:
    const auto es_momentum_range_sapn = es_momentum_domain_variant_to_momentum_range_sapn(
                interpreted_program_options.es_momentum_domain,
                interpreted_program_options.n_sites);
    // [[expect]]:
    if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
        assert(gs_receipt_optional);
    }
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        assert(es_receipts_optional);
        assert(es_receipts_optional->size() == es_momentum_range_sapn.second - es_momentum_range_sapn.first);
    }
    // Print gs:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
        const auto gs_receipt = gs_receipt_optional.value();
        std::cout << " ├state: gs "  << std::endl;
        std::cout << " ││enery = " << gs_receipt.energy << std::endl;
        if (reference_energies) {
            if (const auto& gs_energy = reference_energies->get_gs_energy()){
                std::cout << " ││enery (reference) = " << *gs_energy << std::endl;
            }
        }
    }
    // Print es:
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        const auto es_receipts = es_receipts_optional.value();
        for (unsigned k_n =  es_momentum_range_sapn.first; k_n < es_momentum_range_sapn.second; k_n++) {
            const auto es_result = es_receipts[k_n];
            std::cout << " ├state: es [k_n = " << k_n << "]" << std::endl;
            std::cout << " ││abs. enery        = " << es_result.energy << std::endl;
            if (interpreted_program_options.run_type == RunType::EG) {
                const auto gs_receipt = gs_receipt_optional.value();
                std::cout << " ││exc. enery        = " << es_result.energy - gs_receipt.energy << std::endl;
            }
            if (reference_energies) {
                if (const auto& get_es_absolute_enery = reference_energies->get_es_absolute_enery(k_n)) {
                    std::cout << " ││abs. enery (ref.) = " << *get_es_absolute_enery << std::endl;
                }
                if (const auto& es_exciation_enery = reference_energies->get_es_exciation_enery(k_n)) {
                    std::cout << " ││exc. enery (ref.) = " << *es_exciation_enery << std::endl;
                }
            }
        }
    }
}

void print_post_data(
        const InterpretedProgramOptions& interpreted_program_options,
        /*const std::unique_ptr<model_monostar::ReferenceEnergies> reference_energies,*/
        const std::optional<bfpt_common::CommonRecipeReceipt>& gs_receipt_optional,
        const std::optional<std::vector<bfpt_common::CommonRecipeReceipt>>& es_receipts_optional) {
    // Helpers:
    const auto es_momentum_range_sapn = es_momentum_domain_variant_to_momentum_range_sapn(
                interpreted_program_options.es_momentum_domain,
                interpreted_program_options.n_sites);
    // [[expect]]:
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        assert(es_receipts_optional->size() == es_momentum_range_sapn.second - es_momentum_range_sapn.first);
    }
    // Stream RAII:
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    // Helpers:
    using  extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::operator<<;
    using namespace boost::adaptors;
    // Print gs:
    if (interpreted_program_options.run_type == RunType::G || interpreted_program_options.run_type == RunType::EG) {
        const auto gs_receipt = gs_receipt_optional.value();
        std::cout << "[RESULT] [POST] gs_energy: " << gs_receipt.energy << std::endl;
    }
    // Print es -- domain:
    if (interpreted_program_options.run_type == RunType::EG) {
        const auto& n_sites = interpreted_program_options.n_sites;
        const auto nk_to_k = [n_sites](int n_k)->double {return (2 * arma::datum::pi * n_k) / n_sites;};
        const auto domain = boost::irange(es_momentum_range_sapn.first, es_momentum_range_sapn.second) | transformed(nk_to_k);
        std::cout << "[RESULT] [POST] domain: " << (domain | RSS<double>().like_python_list()) << std::endl;
    }
    // Print es -- helpers:
    const auto receipt_to_energy = [](bfpt_common::CommonRecipeReceipt result)->double {return result.energy;};
    // Print es -- absolute energies:
    if (interpreted_program_options.run_type == RunType::E || interpreted_program_options.run_type == RunType::EG) {
        const auto es_receipt = es_receipts_optional.value();
        std::cout << "[RESULT] [POST] es_absolute_energies: " << (es_receipt | transformed(receipt_to_energy) | RSS<double>().like_python_list()) << std::endl;
    }
    // Print es -- excitation energies:
    if (interpreted_program_options.run_type == RunType::EG) {
        const auto gs_receipt = gs_receipt_optional.value();
        const auto es_receipts = es_receipts_optional.value();
        const auto gs_energy = gs_receipt.energy;
        const auto absolute_energy_into_excitation_energy = [gs_energy] (double es_energy) ->double {return es_energy - gs_energy;};
        const auto exciation_energies = es_receipts | transformed(receipt_to_energy) | transformed(absolute_energy_into_excitation_energy);
        std::cout << "[RESULT] [POST] es_excitation_energies: " << (exciation_energies | RSS<double>().like_python_list()) << std::endl;
    }
}

void print_theta_opt(const HamiltonianParamsFo& hamiltonian_fo_params, std::optional<double> user_defined_overrule) {
    using namespace extension::boost::stream_pragma;
    const extension::std::StreamFromatStacker stream_format_stacker(std::cout);
    using extension::boost::stream_pragma::RSS;
    using extension::boost::stream_pragma::operator|;
    using extension::boost::stream_pragma::operator<<;
    if (const auto & _ = hamiltonian_fo_params.get_theta_opt_analytical()) {
        std::cout << "[INFO   ] [THETA_OPT] optimal orbital theta (analytical) = " << (_.unwrap() | RSS<double>().like_python_set()) << std::endl;
    } else {
        std::cout << "[INFO   ] [THETA_OPT] optimal orbital theta (analytical) = " << "<no known analicycal solution solver>" << std::endl;
    }
    std::cout << "[INFO   ] [THETA_OPT] optimal orbital theta (numerical)  = " << (hamiltonian_fo_params.get_theta_opt_numerical()| RSS<double>().like_python_set()) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] optimal orbital theta              = " << (hamiltonian_fo_params.get_theta_opt() | RSS<double>().like_python_set() ) << std::endl;
    const double orbital_theta_to_use = get_orbital_theta(hamiltonian_fo_params, user_defined_overrule) ; //may thorw!
    std::cout << "[INFO   ] [THETA_OPT] used orbital theta                 = " << orbital_theta_to_use << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] cos, sin of used orbital theta     = " << std::cos(orbital_theta_to_use) << "," << std::sin(orbital_theta_to_use) << std::endl;
    std::cout << "[INFO   ] [THETA_OPT] cos², sin² of used orbital theta   = " << std::pow(std::cos(orbital_theta_to_use), 2) << ", " << std::pow(std::sin(orbital_theta_to_use), 2) << std::endl;
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
                return SpinSiteNamedMatrices::site_matrices_for_average_calculations_af_fm();
            case ModelType::FO:
            {
                const double orbital_theta_to_use = get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
                return OrbitalSiteNamedMatrices::site_matrices_for_average_calculations(orbital_theta_to_use);
            }
            default:
                throw std::domain_error("Invalid model_type enum value.");
            };
        }();
        const auto two_sites_metrices_for_average_calculation =  [&interpreted_program_options]() {
            switch (interpreted_program_options.model_type) {
            case ModelType::AF:
                return SpinTwoSiteNamedMatrices::two_site_matrices_for_average_calculations_af();
            case ModelType::FM:
                return SpinTwoSiteNamedMatrices::two_site_matrices_for_average_calculations_fm();
            case ModelType::FO:
            {
                const double orbital_theta_to_use = get_orbital_theta(interpreted_program_options.hamiltonian_params_fo, interpreted_program_options.orbital_theta);
                return OrbitalTwoSiteNamedMatrices::two_site_matrices_for_average_calculations(orbital_theta_to_use);
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
//        const bool print_density_operator_matrix = true;//TODO remove

// ******************************************************************
print_results_tree(
        interpreted_program_options,
        reference_energies,
        gs_receipt,
        es_receipts);
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
